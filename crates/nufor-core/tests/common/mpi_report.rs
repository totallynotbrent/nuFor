// each test binary compiles this module, so not every helper is used in
// every binary.
#![allow(dead_code)]

//! shared measurement recording for the mpi gates.
//!
//! a gate writes its measurements to the file named by `NF_MPI_REPORT`
//! when the environment names one. the case gate under `cases/mpi-parity`
//! sets that path per invocation, so the committed artefact comes from the
//! same process that performed the check rather than from a later
//! reconstruction.

use std::collections::BTreeMap;
use std::sync::Mutex;
use std::sync::OnceLock;

/// one measurement: a number, or a count of exact-match entries.
#[derive(Clone)]
pub enum Value {
    Number(f64),
    Integer(i64),
    Text(String),
}

impl Value {
    fn render(&self) -> String {
        match self {
            Value::Number(v) => format!("{v:?}"),
            Value::Integer(v) => v.to_string(),
            Value::Text(v) => format!("\"{}\"", v.replace('"', "'")),
        }
    }
}

static RECORDS: OnceLock<Mutex<BTreeMap<String, Value>>> = OnceLock::new();

fn records() -> &'static Mutex<BTreeMap<String, Value>> {
    RECORDS.get_or_init(|| Mutex::new(BTreeMap::new()))
}

/// record one measurement under a name.
pub fn record(name: &str, value: Value) {
    records().lock().unwrap().insert(name.to_string(), value);
}

/// record a floating-point measurement.
pub fn number(name: &str, value: f64) {
    record(name, Value::Number(value));
}

/// record an integer measurement.
pub fn integer(name: &str, value: i64) {
    record(name, Value::Integer(value));
}

/// record a text measurement.
pub fn text(name: &str, value: &str) {
    record(name, Value::Text(value.to_string()));
}

/// the rank count this process runs under.
pub fn world_size() -> i64 {
    std::env::var("OMPI_COMM_WORLD_SIZE")
        .ok()
        .and_then(|s| s.parse().ok())
        .or_else(|| std::env::var("PMI_SIZE").ok().and_then(|s| s.parse().ok()))
        .unwrap_or(1)
}

/// hold every rank at one collective point.
///
/// the launcher tears the job down as soon as the first rank exits, so a
/// rank that returns from a lane early takes rank 0 down with it while
/// rank 0 is still walking its report to disk. the fence makes the other
/// ranks wait until rank 0 has published its file.
pub fn fence() {
    if world_size() > 1 {
        let mut token = [0.0f64];
        let _ = nufor_core::ffi_mpi::mpi_bcast_field(&mut token);
    }
}

/// write the recorded measurements to the file the environment names.
///
/// the write is atomic: the file lands under a temporary name and is
/// renamed into place, so a reader never sees a half-written artefact.
pub fn flush() {
    // every rank shares the report path, so only rank 0 writes it.
    let rank: i64 = std::env::var("OMPI_COMM_WORLD_RANK")
        .ok()
        .and_then(|s| s.parse().ok())
        .or_else(|| std::env::var("PMI_RANK").ok().and_then(|s| s.parse().ok()))
        .unwrap_or(0);
    if rank != 0 {
        return;
    }
    let path = match std::env::var_os("NF_MPI_REPORT") {
        Some(p) => std::path::PathBuf::from(p),
        None => return,
    };
    let mut body = String::new();
    body.push_str("{\n");
    let entries: Vec<(String, Value)> = records()
        .lock()
        .unwrap()
        .iter()
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect();
    for (idx, (k, v)) in entries.iter().enumerate() {
        body.push_str(&format!("  \"{}\": {}", k, v.render()));
        if idx + 1 < entries.len() {
            body.push(',');
        }
        body.push('\n');
    }
    body.push_str("}\n");

    if let Some(parent) = path.parent() {
        let _ = std::fs::create_dir_all(parent);
    }
    // the report lands in one write: the launcher can tear the job down
    // the moment the first rank exits, so the fewer syscalls sit between
    // the measurement and the file, the less chance the job dies with
    // the report still unwritten.
    if let Err(e) = std::fs::write(&path, body.as_bytes()) {
        eprintln!("mpi report: write {} failed: {e}", path.display());
    }
}

//! round-trip and compatibility checks for the binary restart format.

use std::path::PathBuf;

use nufor_core::{read_restart, write_restart, ConservedState, Error, RESTART_VERSION};

fn state() -> ConservedState {
    ConservedState {
        rho: vec![1.0, 1.1, 0.9, 1.0],
        m: vec![0.5, 0.6, 0.4, 0.5],
        e: vec![2.5, 2.6, 2.4, 2.5],
    }
}

fn tmp(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("nufor_rst_{tag}_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn restart_round_trips_exactly_and_reports_version() {
    let dir = tmp("round");
    let path = dir.join("run.rst");
    write_restart(&path, &state(), 1.4, 0.22, 1337).unwrap();

    let got = read_restart(&path).unwrap();
    assert_eq!(got.state, state());
    assert_eq!(got.gamma, 1.4);
    assert_eq!(got.time, 0.22);
    assert_eq!(got.step, 1337);
    assert_eq!(got.version, RESTART_VERSION);
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn old_legacy_layout_is_still_readable_as_version_zero() {
    let dir = tmp("legacy");
    let path = dir.join("legacy.rst");
    // the pre-versioned format: "NUFR1\0" magic, then header, then arrays.
    let mut bytes = Vec::new();
    bytes.extend_from_slice(b"NUFR1\0");
    bytes.extend_from_slice(&4u64.to_le_bytes());
    bytes.extend_from_slice(&1.4f64.to_le_bytes());
    bytes.extend_from_slice(&0.22f64.to_le_bytes());
    bytes.extend_from_slice(&1337u64.to_le_bytes());
    for v in [&state().rho, &state().m, &state().e] {
        for x in v {
            bytes.extend_from_slice(&x.to_le_bytes());
        }
    }
    std::fs::write(&path, &bytes).unwrap();

    let got = read_restart(&path).unwrap();
    assert_eq!(got.version, 0);
    assert_eq!(got.state, state());
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn future_schema_version_is_rejected_as_incompatible() {
    let dir = tmp("future");
    let path = dir.join("future.rst");
    let mut bytes = Vec::new();
    bytes.extend_from_slice(b"NUFR");
    bytes.push(RESTART_VERSION + 1); // a schema this build does not yet understand.
    bytes.extend_from_slice(&4u64.to_le_bytes());
    bytes.extend_from_slice(&1.4f64.to_le_bytes());
    bytes.extend_from_slice(&0.0f64.to_le_bytes());
    bytes.extend_from_slice(&0u64.to_le_bytes());
    bytes.extend_from_slice(&[0u8; 4 * 8]);
    bytes.extend_from_slice(&[0u8; 4 * 8]);
    bytes.extend_from_slice(&[0u8; 4 * 8]);
    std::fs::write(&path, &bytes).unwrap();

    assert_eq!(read_restart(&path), Err(Error::IncompatibleVersion));
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn wrong_magic_or_short_or_oversized_files_are_rejected() {
    let dir = tmp("bad");
    // not a restart at all.
    let bad = dir.join("bad.rst");
    std::fs::write(&bad, b"NOTARESTARTFILE").unwrap();
    assert!(read_restart(&bad).is_err());
    // versioned magic but truncated before the header.
    let short = dir.join("short.rst");
    std::fs::write(&short, b"NUFR\x01").unwrap();
    assert!(read_restart(&short).is_err());
    // valid round-trip size, printed, then one stray byte appended.
    let path = dir.join("tampered.rst");
    write_restart(&path, &state(), 1.4, 0.22, 1337).unwrap();
    {
        let mut bytes = std::fs::read(&path).unwrap();
        bytes.push(0);
        std::fs::write(&path, &bytes).unwrap();
    }
    assert!(read_restart(&path).is_err());
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn non_finite_physics_in_a_well_shaped_file_is_rejected() {
    let dir = tmp("nan");
    let path = dir.join("nan.rst");
    let mut bytes = Vec::new();
    bytes.extend_from_slice(b"NUFR");
    bytes.push(RESTART_VERSION);
    bytes.extend_from_slice(&4u64.to_le_bytes());
    bytes.extend_from_slice(&f64::NAN.to_le_bytes());
    bytes.extend_from_slice(&0.0f64.to_le_bytes());
    bytes.extend_from_slice(&0u64.to_le_bytes());
    bytes.extend_from_slice(&[0u8; 4 * 8]);
    bytes.extend_from_slice(&[0u8; 4 * 8]);
    bytes.extend_from_slice(&[0u8; 4 * 8]);
    std::fs::write(&path, &bytes).unwrap();

    assert!(read_restart(&path).is_err());
    std::fs::remove_dir_all(&dir).unwrap();
}

//! smoke tests that drive the built `nufor` binary end to end.

use std::path::PathBuf;
use std::process::{Command, Output};

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_nufor")
}

fn run(args: &[&str], cwd: &PathBuf) -> Output {
    Command::new(bin())
        .current_dir(cwd)
        .args(args)
        .output()
        .unwrap()
}

fn stdout(o: &Output) -> String {
    String::from_utf8_lossy(&o.stdout).into_owned()
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("nufor_cli_{tag}_{}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    d
}

#[test]
fn version_prints_a_version() {
    let d = tmpdir("v");
    let o = run(&["version"], &d);
    assert!(stdout(&o).contains("nufor"));
    std::fs::remove_dir_all(&d).unwrap();
}

#[test]
fn mesh_reports_the_grid() {
    let d = tmpdir("m");
    let o = run(&["mesh", "8", "0", "8"], &d);
    assert!(stdout(&o).contains("cells: 8"));
    assert!(stdout(&o).contains("dx: 1"));
    std::fs::remove_dir_all(&d).unwrap();
}

#[test]
fn run_writes_a_restart_and_inspect_reads_it() {
    let d = tmpdir("r");
    let rst: PathBuf = d.join("run.rst");
    let run_out = run(&["run", "100", "0.2", "sod"], &d);
    assert!(
        run_out.stderr.is_empty(),
        "run stderr: {}",
        String::from_utf8_lossy(&run_out.stderr)
    );
    // the run above writes no restart unless the path is given; run again with it.
    let run2 = run(&["run", "100", "0.2", "sod", rst.to_str().unwrap()], &d);
    assert!(run2.status.success());
    let o = run(&["inspect", rst.to_str().unwrap()], &d);
    assert!(stdout(&o).contains("cells: 100"));
    std::fs::remove_dir_all(&d).unwrap();
}

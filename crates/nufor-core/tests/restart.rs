//! round-trip and validation checks for the binary restart format.

use std::path::PathBuf;

use nufor_core::{read_restart, write_restart, ConservedState};

#[test]
fn restart_round_trips_exactly() {
    let state = ConservedState {
        rho: vec![1.0, 1.1, 0.9, 1.0],
        m: vec![0.5, 0.6, 0.4, 0.5],
        e: vec![2.5, 2.6, 2.4, 2.5],
    };
    let dir = std::env::temp_dir().join(format!("nufor_rst_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path: PathBuf = dir.join("run.rst");
    write_restart(&path, &state, 1.4, 0.22, 1337).unwrap();

    let got = read_restart(&path).unwrap();
    assert_eq!(got.state, state);
    assert_eq!(got.gamma, 1.4);
    assert_eq!(got.time, 0.22);
    assert_eq!(got.step, 1337);
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn restart_rejects_a_corrupt_or_wrong_file() {
    let dir = std::env::temp_dir().join(format!("nufor_rst_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    // wrong magic.
    let bad = dir.join("bad.rst");
    std::fs::write(&bad, b"NOTARESTARTFILE").unwrap();
    assert!(read_restart(&bad).is_err());
    // truncated header.
    let short = dir.join("short.rst");
    std::fs::write(&short, b"NUFR1\0\x02\x00").unwrap();
    assert!(read_restart(&short).is_err());
    std::fs::remove_dir_all(&dir).unwrap();
}

#!/usr/bin/env bash
# the regression gate: format, clippy, the whole test suite, the mpi
# parity gates, and the repo hygiene check must all pass. exits
# non-zero if any group fails. run from the repo root.
#
# every group runs even when an earlier one fails, so one invocation
# reports every problem in it. ci runs the same groups in
# .github/workflows/ci.yml.
set -u

export PATH="$HOME/.cargo/bin:$PATH"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$here"

fails=0

echo "== fmt --check =="
cargo fmt --check >/dev/null 2>&1 || { echo "FAIL: formatting"; fails=$((fails+1)); }

echo "== clippy -D warnings =="
cargo clippy --workspace --all-targets -- -D warnings >/tmp/reg_clippy 2>&1 || { echo "FAIL: clippy"; fails=$((fails+1)); }

echo "== test --workspace =="
# one run, one log: the suite is the slowest group here and running it
# twice doubled the wall time for no extra coverage.
test_log=/tmp/reg_test
if cargo test --workspace >"$test_log" 2>&1; then
  suites=$(grep -cE "test result: ok" "$test_log")
  passed=$(grep -oE "ok\. [0-9]+ passed" "$test_log" | awk '{s += $2} END {print s + 0}')
  ignored=$(grep -oE "[0-9]+ ignored" "$test_log" | awk '{s += $1} END {print s + 0}')
  echo "  $suites suites, $passed tests passed, $ignored ignored"
else
  echo "FAIL: tests"
  grep -E 'FAILED|error\[|panicked' "$test_log" | head -5
  fails=$((fails+1))
  suites=0
fi

echo "== mpi parity (one gate per process) =="
# each mpi gate runs in its own process, the way the ci runs them. the
# gates share one mpi communicator, and a process that runs them back to
# back carries the previous gate's collective state into the next one:
# the sa gate then lands on the east domain edge with a drift of 4.8e-5
# instead of its own 8.2e-14. one gate per process keeps every gate's
# communicator state its own.
mpi_bin=$(ls -t target/debug/deps/mpi_parity-* 2>/dev/null | grep -v '\.d$' | head -1)
if [ -z "$mpi_bin" ]; then
  cargo build --tests >/dev/null 2>&1
  mpi_bin=$(ls -t target/debug/deps/mpi_parity-* 2>/dev/null | grep -v '\.d$' | head -1)
fi
if [ -n "$mpi_bin" ]; then
  for gate in mpi_np1_matches_rust_sod2d mpi_np1_matches_rust_uniform_eqair \
              mpi_np1_visc_matches_rust mpi_np1_matches_rust_capsule \
              mpi_np1_sa_matches_rust; do
    if ! "$mpi_bin" --exact "$gate" --nocapture --ignored >/tmp/reg_mpi 2>&1; then
      echo "FAIL: mpi gate $gate"
      grep -E "worst rel diff|panicked" /tmp/reg_mpi | head -3
      fails=$((fails+1))
    fi
  done
  echo "  mpi gates: 5 run in their own process"
fi

echo "== case gates =="
# the shared runner owns the gate list, so the local gate and the ci
# gate run the same cases. a gate whose committed artifacts are in the
# repository must run here; a gate that needs a marched field reports a
# skip when the field is not on disk.
if ! bash tools/case_gates.sh; then
  echo "FAIL: case gates"
  fails=$((fails+1))
fi

echo "== repo hygiene =="
if git status --porcelain | grep -q .; then echo "FAIL: dirty tree"; fails=$((fails+1)); fi
# the private-path rule lives in tools/repo_hygiene.py so ci runs the
# same check. docs/ and docs-site/ are published product, not private
# scratch, so they are allowed to be tracked.
if ! python3 tools/repo_hygiene.py; then echo "FAIL: repo hygiene"; fails=$((fails+1)); fi

if [ "$fails" -eq 0 ]; then
  echo "regression gate: PASS ($suites suites green)"
  exit 0
fi
echo "regression gate: FAIL ($fails check(s) failing)"
exit 1
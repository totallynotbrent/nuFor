#!/usr/bin/env bash
# the regression gate: format, clippy, and the whole test suite must pass.
# exits non-zero if any check fails. run from the repo root.
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
suites=$(cargo test --workspace 2>&1 | grep -cE "test result: ok")
failed=$(cargo test --workspace 2>&1 | grep -cE "FAILED|error\[")
echo "  passing suites: $suites, failing/errors: $failed"
if [ "$suites" -eq 0 ] || [ "$failed" -ne 0 ]; then
  echo "FAIL: tests"; fails=$((fails+1))
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

echo "== repo hygiene =="
if git status --porcelain | grep -q .; then echo "FAIL: dirty tree"; fails=$((fails+1)); fi
for p in plans/ goals/ nuFor_README.md docs/research docs-site; do
  if git ls-files | grep -q "^$p"; then echo "FAIL: tracked private path $p"; fails=$((fails+1)); fi
done

if [ "$fails" -eq 0 ]; then
  echo "regression gate: PASS ($suites suites green)"
  exit 0
fi
echo "regression gate: FAIL ($fails check(s) failing)"
exit 1
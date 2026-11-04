#!/usr/bin/env bash
# rebuild every published figure of the chemistry page from the repo's
# own data.
#
# two figure families live here. the verification figures come from the
# crate's example, which compares the generated mechanism's net rates of
# progress against the checked-in cantera fixture. the ignition figures
# come from the marched parcel's own trace against the cantera oracle
# beside the case.
#
# both are reproductions of the measurement rather than drawings of it,
# and the plotter needs the system's libraries, which do not live in an
# agent's virtual environment, so it runs with the system interpreter and
# a cleaned PYTHONPATH.
#
# usage:
#   bash regen-figures.sh [out-dir]
set -eu

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
OUT=${1:-"$REPO/docs/img"}
WORK="${NUFOR_FIGURE_WORK:-$HOME/.hermes/cache/scratch/chemistry-figures}"
CASE="$HERE/gri30-stoich-ignition"
mkdir -p "$OUT" "$WORK"

echo "the chemistry page's figures, rebuilt from $REPO"
echo "output: $OUT"

# the mechanism's own numbers, against the cantera fixture in the crate
CSV="$WORK/verification.csv"
echo "generating the verification csv"
( cd "$REPO" && cargo run --quiet --release -p nufor-core \
    --example gri30_verification ) > "$CSV"
env -u PYTHONPATH /usr/bin/python3 "$HERE/figures.py" "$CSV" "$OUT"

# the parcel's own trace, against the oracle beside the case
if [ ! -f "$CASE/ignition-trace.csv" ]; then
  echo "missing $CASE/ignition-trace.csv; march it first:" >&2
  echo "  nufor ignite $CASE/case.toml" >&2
  exit 1
fi
env -u PYTHONPATH /usr/bin/python3 "$HERE/ignition_figures.py" "$CASE" "$OUT"

echo
echo "done"

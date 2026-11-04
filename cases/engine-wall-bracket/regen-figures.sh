#!/usr/bin/env bash
# rebuild every published figure of the engine wall bracket from one
# campaign's data.
#
# each script reads the campaign's marched fields and progress records,
# so the pictures are reproductions of the measurement rather than
# drawings of it. the gate's verdict and the figures come from the same
# code path, so a picture cannot quietly disagree with the gate.
#
# usage:
#   bash regen-figures.sh <runs root> [out-dir]
#
# the figures themselves need the system's plotting libraries, which do
# not live in an agent's virtual environment, so they are invoked with
# the system interpreter and a cleaned PYTHONPATH. the instruments that
# judge the data need python3.12, because Cantera does not load under
# the older interpreter on this host.
set -eu

RUNS=${1:?usage: regen-figures.sh <runs root> [out-dir]}
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
OUT=${2:-"$REPO/docs/img"}
mkdir -p "$OUT"

echo "the bracket's figures, rebuilt from $RUNS"
echo "output: $OUT"

# the four long-standing pictures of the wall's behaviour
env -u PYTHONPATH /usr/bin/python3 "$HERE/figures.py" "$RUNS" "$OUT"

# the settlement criterion and its window
env -u PYTHONPATH /usr/bin/python3 "$HERE/settlement_figures.py" "$RUNS" "$OUT"

# the wall's effect on the gas, and the chamber's common mode
env -u PYTHONPATH /usr/bin/python3 "$HERE/effect_figures.py" "$RUNS" "$OUT"

echo
echo "done"

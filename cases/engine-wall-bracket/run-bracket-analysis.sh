#!/usr/bin/env bash
# run the whole bracket's instruments over a campaign's data and print
# one report.
#
# the three instruments share the campaign's period, so running them
# together is what lets their numbers be compared. the gate is the
# verdict, the analyser reports the observables, and the near-wall
# instrument asks the gas whether it agrees with the wall's claim.
#
# usage:
#   bash run-bracket-analysis.sh <runs root>
set -u

RUNS=${1:?usage: run-bracket-analysis.sh <runs root>}
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"

echo "======================================================================"
echo "the engine wall bracket"
echo "runs root: $RUNS"
echo "======================================================================"

echo
echo "----------------------------------------------------------------------"
echo "the gate's verdict"
echo "----------------------------------------------------------------------"
( cd "$HERE" && NUFOR_BRACKET_RUNS="$RUNS" python3.12 gates.py )
GATE_RC=$?

echo
echo "----------------------------------------------------------------------"
echo "the analyser's numbers"
echo "----------------------------------------------------------------------"
( cd "$HERE" && python3.12 analyse.py "$RUNS" --json "$RUNS/bracket-results.json" )
ANALYSE_RC=$?

echo
echo "----------------------------------------------------------------------"
echo "the near-wall gas, against the adiabatic leg"
echo "----------------------------------------------------------------------"
( cd "$HERE" && python3.12 near_wall_response.py "$RUNS" )
NEAR_RC=$?

echo
echo "======================================================================"
echo "the gate exited $GATE_RC, the analyser $ANALYSE_RC, the near-wall $NEAR_RC"
echo "a non-zero gate exit is the bracket's verdict, not a malfunction."
echo "======================================================================"
exit $GATE_RC

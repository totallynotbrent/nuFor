#!/usr/bin/env bash
# the case gate suite: one list of case gates, run the same way by ci
# and by the local regression script.
#
# each case carries a mode:
#   committed  the gate reads artifacts that live in the repository, so
#              it must run and pass everywhere. a skip here is a failure.
#   field      the gate reads the marched field next to the case file.
#              that field is a local artifact, so the gate runs where the
#              field exists and reports a skip where it does not.
#
# run from the repo root:
#   bash tools/case_gates.sh
#
# exit status is non-zero when any gate failed.
set -u

here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$here"

case_gates=(
  "cases/cylinder-m20-park5 committed"
  "cases/cylinder-stagnation-fay committed"
  "cases/cylinder-m20-park5-visc committed"
  "cases/chamber-throat committed"
  "cases/nozzle-cdv field"
  "cases/aeroshell-station field"
)

fails=0
for entry in "${case_gates[@]}"; do
  case_dir="${entry%% *}"
  mode="${entry##* }"
  if [ ! -f "$case_dir/gates.py" ]; then
    echo "FAIL $case_dir: no gates.py"
    fails=$((fails + 1))
    continue
  fi
  out=$(cd "$case_dir" && python3 gates.py 2>&1)
  rc=$?
  if [ $rc -eq 0 ]; then
    echo "PASS $case_dir"
    continue
  fi
  if printf '%s' "$out" | grep -qiE "no marched field|not marched|run the case|MISS"; then
    if [ "$mode" = "field" ]; then
      echo "SKIP $case_dir: no marched field here"
    else
      echo "FAIL $case_dir: the gate needs artifacts that should be in the repository"
      printf '%s\n' "$out" | tail -6
      fails=$((fails + 1))
    fi
    continue
  fi
  echo "FAIL $case_dir"
  printf '%s\n' "$out" | tail -8
  fails=$((fails + 1))
done

exit $fails

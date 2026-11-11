#!/usr/bin/env bash
# install the linux toolchain a job needs, then prove it is usable.
#
# the github runners sometimes lose dns for the ubuntu archive halfway
# through an install. retrying covers the transient case. when the retry
# still fails, the job only dies if a tool is genuinely missing: the
# runner images already carry some of these, and a step that insists on
# reinstalling what is already there turns a network wobble into a
# failed build.
#
# usage:
#   .github/scripts/linux_toolchain.sh <package> [<package>...]
#
# the packages are installed with apt. the tools are then checked by
# running them.
set -uo pipefail

attempts="${APT_ATTEMPTS:-5}"
delay="${APT_DELAY_SECONDS:-15}"

packages=("$@")
if [ "${#packages[@]}" -eq 0 ]; then
  echo "usage: $0 <package> [<package>...]" >&2
  exit 2
fi

install_once() {
  sudo apt-get update && sudo apt-get install -y --no-install-recommends "${packages[@]}"
}

installed=1
for i in $(seq 1 "$attempts"); do
  if install_once; then
    installed=0
    break
  fi
  echo "apt attempt $i of $attempts failed, retrying in ${delay}s"
  sleep "$delay"
done
if [ "$installed" -ne 0 ]; then
  echo "the archive stayed unreachable after $attempts attempts"
fi

# the tools a job actually runs. each check is the cheapest proof that the
# tool is present and usable, not just that a package unpacked.
require() {
  local label="$1"
  shift
  if "$@" > /dev/null 2>&1; then
    echo "OK $label"
    return 0
  fi
  echo "MISSING $label: $* did not run"
  MISSING_TOOLS="$MISSING_TOOLS $label"
  return 1
}

missing=""
mf() {
  MISSING_TOOLS=""
  require "gfortran" gfortran --version || missing="$missing gfortran"
  require "cmake" cmake --version || missing="$missing cmake"
  require "numpy" python3 -c "import numpy" || missing="$missing numpy"
  require "hdf5 headers" h5cc -showconfig || missing="$missing hdf5"

  if printf '%s\n' "${packages[@]}" | grep -qx "openmpi-bin"; then
    require "mpirun" mpirun --version || missing="$missing mpirun"
    require "mpicc" mpicc --version || missing="$missing mpicc"
  fi

  if [ -n "$missing" ]; then
    echo "these tools are unavailable:$missing"
    exit 1
  fi
  echo "toolchain ready"
}

mf

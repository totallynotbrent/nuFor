#!/usr/bin/env bash
# prepare the linux toolchain a job needs, then prove it is usable.
#
# the github runners sometimes lose dns for the ubuntu archive, and an
# apt call against an unreachable mirror hangs rather than failing
# cleanly. the runner image already carries much of what a job needs, so
# the check comes first: when the tools run, the job goes ahead without
# touching apt at all. only a genuinely missing tool leads to an
# install, and that install is bounded, so a dead mirror fails the job
# in minutes instead of hanging it for an hour.
#
# usage:
#   .github/scripts/linux_toolchain.sh <package> [<package>...]
#
# the arguments are apt packages. an openmpi package also demands
# mpirun and mpicc.
set -uo pipefail

attempts="${APT_ATTEMPTS:-3}"
delay="${APT_DELAY_SECONDS:-10}"

packages=("$@")
if [ "${#packages[@]}" -eq 0 ]; then
  echo "usage: $0 <package> [<package>...]" >&2
  exit 2
fi

# the tools a job actually runs, each paired with the cheapest proof that
# it is present and usable rather than merely unpacked.
TOOLS=(
  "gfortran:gfortran --version"
  "cmake:cmake --version"
  "numpy:numpy_check"
  "hdf5:h5cc -showconfig"
)
OPENMPI_PACKAGES="libopenmpi-dev openmpi-bin"

wants_openmpi=0
for pkg in "${packages[@]}"; do
  for omp in $OPENMPI_PACKAGES; do
    [ "$pkg" = "$omp" ] && wants_openmpi=1
  done
done
if [ "$wants_openmpi" -eq 1 ]; then
  TOOLS+=("mpirun:mpirun --version" "mpicc:mpicc --version")
fi

runs() {
  case "$1" in
    numpy_check)
      python3 -c "import numpy" > /dev/null 2>&1
      ;;
    *)
      $1 > /dev/null 2>&1
      ;;
  esac
}

# the names of the tools that do not run right now.
missing_tools() {
  local entry name cmd
  for entry in "${TOOLS[@]}"; do
    name="${entry%%:*}"
    cmd="${entry#*:}"
    runs "$cmd" || echo "$name"
  done
}

missing="$(missing_tools)"
if [ -z "$missing" ]; then
  echo "the toolchain is already in place:"
  for entry in "${TOOLS[@]}"; do
    echo "  OK ${entry%%:*}"
  done
  exit 0
fi

echo "missing before the install:${missing}"

install_once() {
  # the acquire timeouts keep a dead mirror from stalling the job: apt
  # gives up on an unreachable host instead of waiting it out.
  sudo apt-get \
    -o Acquire::http::Timeout=20 \
    -o Acquire::https::Timeout=20 \
    -o Acquire::Retries=1 \
    update &&
    sudo apt-get \
      -o Acquire::http::Timeout=20 \
      -o Acquire::https::Timeout=20 \
      -o Acquire::Retries=1 \
      install -y --no-install-recommends "${packages[@]}"
}

for i in $(seq 1 "$attempts"); do
  if install_once; then
    break
  fi
  if [ "$i" -lt "$attempts" ]; then
    echo "apt attempt $i of $attempts failed, retrying in ${delay}s"
    sleep "$delay"
  fi
done

missing="$(missing_tools)"
if [ -n "$missing" ]; then
  echo "these tools are unavailable:${missing}"
  exit 1
fi

echo "the toolchain is ready:"
for entry in "${TOOLS[@]}"; do
  echo "  OK ${entry%%:*}"
done

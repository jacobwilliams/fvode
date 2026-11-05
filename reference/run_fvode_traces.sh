#!/usr/bin/env bash
# Build and run the FVODE step-trace drivers (reference/fortran/trace_*.f90),
# each of which writes reference/out/fvode_<name>.csv.
# Extra arguments are passed to fpm (e.g. --flag "-ffp-contract=off").
# Run through pixi:  pixi run -e reference fvode-traces
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE/.."
mkdir -p reference/out

for f in reference/fortran/trace_*.f90; do
    name="$(basename "${f%.f90}")"
    case "$name" in trace_common|*_driver) continue ;; esac
    fpm run --example "$name" --profile release "$@" > /dev/null
    echo "wrote reference/out/fvode_${name#trace_}.csv"
done

#!/usr/bin/env bash
# Build and run the reference C programs in reference/ against the pinned
# SUNDIALS build, writing their output to reference/out/.
# Run through pixi:  pixi run -e reference ref-run
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
SRC="$HERE/sundials"
INS="$HERE/install"
OUT="$HERE/out"
mkdir -p "$OUT"

CC="${CC:-cc}"
CFLAGS="-O2 -ffp-contract=off -I$INS/include -I$SRC/src/cvode -I$SRC/src/sundials"
LIBDIR="$INS/lib"; [ -d "$INS/lib64" ] && LIBDIR="$INS/lib64"
LIBS="-L$LIBDIR -lsundials_cvode -lsundials_nvecserial -lsundials_sunmatrixdense \
 -lsundials_sunmatrixband -lsundials_sunlinsoldense -lsundials_sunlinsolband \
 -lsundials_sunlinsolspgmr -lsundials_sunlinsolspfgmr -lsundials_sunlinsolspbcgs \
 -lsundials_sunlinsolsptfqmr -lsundials_sunlinsolpcg -lsundials_sunnonlinsolnewton \
 -lsundials_sunnonlinsolfixedpoint -lsundials_core -lm"

# 1. constants table (header values from C, .c-local #defines parsed by python)
"$CC" $CFLAGS "$HERE/constants.c" -o "$OUT/constants"
{
    "$OUT/constants"
    python -I "$HERE/parse_defines.py" "$SRC"
} | LC_ALL=C sort > "$ROOT/test/data/sundials_constants.txt"
echo "wrote test/data/sundials_constants.txt"

# 2. step-trace drivers (trace_*.c), one CSV trace per driver
for f in "$HERE"/trace_*.c; do
    [ -e "$f" ] || continue
    exe="$OUT/$(basename "${f%.c}")"
    "$CC" $CFLAGS "$f" -o "$exe" $LIBS
    "$exe" > "$exe.csv"
    echo "wrote ${exe#$ROOT/}.csv"
done

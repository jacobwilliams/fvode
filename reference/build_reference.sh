#!/usr/bin/env bash
# Fetch and build the pinned reference SUNDIALS (C) used to verify FVODE.
# Run through pixi:  pixi run -e reference ref-build
set -euo pipefail

TAG=v7.9.0
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/sundials"
BLD="$HERE/build"
INS="$HERE/install"

if [ ! -d "$SRC" ]; then
    git clone --depth 1 --branch "$TAG" https://github.com/LLNL/sundials.git "$SRC"
fi

if [ -f "$INS/.built-$TAG" ]; then
    echo "reference SUNDIALS $TAG already built in $INS"
    exit 0
fi

# -ffp-contract=off: the reference does not fuse multiply-adds, so an FVODE
# build with the same flag can be compared exactly (pixi task compare-exact).
# Normal FVODE builds are optimized and compared with a tolerance.
cmake -S "$SRC" -B "$BLD" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_FLAGS="-ffp-contract=off" \
    -DCMAKE_INSTALL_PREFIX="$INS" \
    -DBUILD_SHARED_LIBS=OFF \
    -DBUILD_STATIC_LIBS=ON \
    -DBUILD_ARKODE=OFF -DBUILD_CVODES=OFF -DBUILD_IDA=OFF -DBUILD_IDAS=OFF \
    -DBUILD_KINSOL=OFF -DBUILD_CVODE=ON \
    -DEXAMPLES_ENABLE_C=ON -DEXAMPLES_INSTALL=OFF \
    -DSUNDIALS_PRECISION=double -DSUNDIALS_INDEX_SIZE=64
cmake --build "$BLD"
cmake --install "$BLD"
touch "$INS/.built-$TAG"

#!/bin/bash
# Cross-build yafu for aarch64 on an x86 Linux/WSL host, to run under
# qemu-aarch64-static for correctness checks of the NEON kernels.
#
# Needs: gcc-aarch64-linux-gnu, qemu-user-static, rsync, and a static
# aarch64 GMP, e.g. from a GMP source tree:
#   ./configure --host=aarch64-linux-gnu --prefix=$HOME/neon/gmp-arm64 --disable-shared
#   make -j && make install
#
# Usage: tools/neon_cross_build.sh [make args]     (e.g. FORCE_GENERIC=1)
# Env:   DST       build directory     (default ~/neon/yafu)
#        GMP_ARM   aarch64 GMP prefix  (default ~/neon/gmp-arm64)
#
# The source tree is synced into DST (without config.mk, which holds the
# host's x86 settings, and without host objects), so the checkout itself
# is never touched.
set -e
SRC=$(cd "$(dirname "$0")/.." && pwd)
DST=${DST:-$HOME/neon/yafu}
GMP_ARM=${GMP_ARM:-$HOME/neon/gmp-arm64}

mkdir -p "$DST"
rsync -a --delete \
    --exclude .git --exclude .claude --exclude config.mk \
    --exclude '*.exe' --exclude '*.a' --exclude '*.o' --exclude '*.d' \
    --exclude '*.qo' --exclude '*.no' --exclude '*.obj' \
    --exclude build.vc22 --exclude 'siqs_tune_out*' --exclude test1 \
    --exclude __pycache__ --exclude bin --exclude yafu --exclude '*.log' \
    --exclude '*.dat' \
    "$SRC/" "$DST/"
cd "$DST"

# the Makefile calls plain ar/ranlib; point them at the cross binutils
mkdir -p "$DST/.crossbin"
ln -sf "$(command -v aarch64-linux-gnu-ar)" "$DST/.crossbin/ar"
ln -sf "$(command -v aarch64-linux-gnu-ranlib)" "$DST/.crossbin/ranlib"
export PATH="$DST/.crossbin:$PATH"

rm -f lib*.a
make -j"$(nproc)" yafu CC=aarch64-linux-gnu-gcc GMP_PREFIX="$GMP_ARM" \
    ECM_INCDIR=/nonexistent ECM= USER_CFLAGS="-mcpu=neoverse-n1" \
    USER_LDFLAGS="-static" "$@" > build.log 2>&1 \
    || { grep -E "error|Error" build.log | head -40; exit 1; }
echo "warnings: $(grep -c 'warning:' build.log)"
file yafu | cut -c1-100

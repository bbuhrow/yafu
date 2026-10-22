#!/bin/bash
# Build yafu once per ISA level, each in its own copy of the source tree
# (objects are built in-tree, so ISA levels can't share one tree).
#
# usage: tools/build_matrix.sh [tag ...]
#   tags: avx512 avx2 sse41 sse2 generic   (default: all)
#
# Each tag builds into $BUILD_ROOT/<tag>/src and copies the binary to
# $BUILD_ROOT/<tag>/yafu.  Copies are refreshed with rsync, so rebuilds are
# incremental.  Extra make arguments can be passed in $MAKE_ARGS, e.g.
#   MAKE_ARGS="USER_CFLAGS=-DDIRECT_SCATTER_MODE=2" tools/build_matrix.sh avx512
#
# Environment:
#   BUILD_ROOT  where the copies go          (default ~/yafu_matrix)
#   GMP_PREFIX  GMP install prefix, absolute (default <src>/../gmp-install/6.2.1)
#   ECM_PREFIX  GMP-ECM install prefix, absolute; ECM is linked only if set
#               or if <src>/../ecm-install/7.0.6-wsl exists
#   JOBS        parallel make jobs           (default nproc)

set -e

SRC=$(cd "$(dirname "$0")/.." && pwd)
BUILD_ROOT=${BUILD_ROOT:-$HOME/yafu_matrix}
GMP_PREFIX=${GMP_PREFIX:-$(cd "$SRC/.." && pwd)/gmp-install/6.2.1}
if [ -z "$ECM_PREFIX" ] && [ -d "$SRC/../ecm-install/7.0.6-wsl" ]; then
    ECM_PREFIX=$(cd "$SRC/../ecm-install/7.0.6-wsl" && pwd)
fi
JOBS=${JOBS:-$(nproc)}

TAGS="$*"
[ -z "$TAGS" ] && TAGS="avx512 avx2 sse41 sse2 generic"

isa_flags() {
    # every level sets all the USE_* flags explicitly so that config.mk's
    # choices (e.g. USE_AVX512IFMA := 1) are overridden
    local off="USE_AVX512IFMA=0 USE_AVX512PF=0 USE_NATIVE="
    case "$1" in
        avx512)  echo "$off USE_AVX512=1" ;;
        avx2)    echo "$off USE_AVX512=0 USE_AVX2=1 USE_BMI2=1" ;;
        sse41)   echo "$off USE_AVX512=0 USE_AVX2=0 USE_BMI2=0 USE_SSE41=1" ;;
        sse2)    echo "$off USE_AVX512=0 USE_AVX2=0 USE_BMI2=0 USE_SSE41=0" ;;
        generic) echo "$off USE_AVX512=0 USE_AVX2=0 USE_BMI2=0 USE_SSE41=0 FORCE_GENERIC=1" ;;
        *) echo "unknown tag $1" >&2; return 1 ;;
    esac
}

for tag in $TAGS; do
    flags=$(isa_flags "$tag")
    dir="$BUILD_ROOT/$tag"
    mkdir -p "$dir/src"
    rsync -a --delete \
        --exclude='.git' --exclude='*.o' --exclude='*.no' --exclude='*.qo' \
        --exclude='*.a' --exclude='.deps' --exclude='/yafu' --exclude='/yafu.exe' \
        --exclude='/yafu-x64.exe' --exclude='/siqs_demo.exe' --exclude='/test1' \
        --exclude='siqs.dat' --exclude='*.log' --exclude='/build.vc22' \
        "$SRC/" "$dir/src/"
    ecm="ECM="
    [ -n "$ECM_PREFIX" ] && ecm="ECM=1 ECM_PREFIX=$ECM_PREFIX"
    echo "=== $tag: make yafu $flags $MAKE_ARGS"
    if make -C "$dir/src" yafu -j"$JOBS" $flags GMP_PREFIX="$GMP_PREFIX" $ecm \
        $MAKE_ARGS > "$dir/build.log" 2>&1; then
        cp "$dir/src/yafu" "$dir/yafu"
        echo "    ok -> $dir/yafu"
    else
        echo "    FAILED, see $dir/build.log"
        grep -E "error|undefined reference" "$dir/build.log" | head -20
    fi
done

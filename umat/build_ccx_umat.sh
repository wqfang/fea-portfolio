#!/usr/bin/env bash
# Build CalculiX (ccx) from source with an Abaqus-format UMAT linked in.
#
# Mirrors the Homebrew formula in costerwi/homebrew-calculix (same sources, patches
# and flags), but replaces ccx's empty umat.f template with your routine and
# compiles it with -DCCX_COROT so the UMAT returns the corotated Cauchy stress that
# ccx's ABAQUSNL wrapper expects. In a deck, select it with a material name that
# starts with ABAQUSNL (finite strain) and *USER MATERIAL, CONSTANTS=...
#
# Usage:
#   ./build_ccx_umat.sh                          # src/umat_neohooke.f, ccx 2.23
#   ./build_ccx_umat.sh --umat src/other.f       # a different UMAT
#   ./build_ccx_umat.sh --ccx-tarball ~/Downloads/ccx_2.23.src.tar.bz2 \
#                       --spooles-tarball ~/Downloads/spooles.2.2.tgz   # offline
#
# Needs (macOS): brew install gcc arpack pkg-config
# Output: bin/ccx_<version>_umat
set -euo pipefail

VER=2.23
HERE=$(cd "$(dirname "$0")" && pwd)
UMAT="$HERE/src/umat_neohooke.f"
CCX_TGZ=""
SPOOLES_TGZ=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --umat) UMAT=$(cd "$(dirname "$2")" && pwd)/$(basename "$2"); shift 2 ;;
    --version) VER=$2; shift 2 ;;
    --ccx-tarball) CCX_TGZ=$2; shift 2 ;;
    --spooles-tarball) SPOOLES_TGZ=$2; shift 2 ;;
    -h|--help) sed -n '2,19p' "$0"; exit 0 ;;
    *) echo "unknown option $1" >&2; exit 1 ;;
  esac
done

BUILD="$HERE/build/ccx_$VER"
JOBS=$( (sysctl -n hw.ncpu 2>/dev/null || nproc) )
FC=$(command -v gfortran || ls /opt/homebrew/bin/gfortran-* 2>/dev/null | sort -V | tail -1 || true)
[[ -n "$FC" ]] || { echo "gfortran not found: brew install gcc" >&2; exit 1; }
[[ -f "$UMAT" ]] || { echo "UMAT file not found: $UMAT" >&2; exit 1; }
echo "ccx $VER | UMAT $UMAT | gfortran $FC | $JOBS jobs"

mkdir -p "$BUILD"
cd "$BUILD"

# 1. Sources (downloaded once, reused on later builds)
if [[ -z "$CCX_TGZ" ]]; then
  CCX_TGZ="$BUILD/ccx_$VER.src.tar.bz2"
  [[ -f "$CCX_TGZ" ]] || curl -fL -o "$CCX_TGZ" "https://www.dhondt.de/ccx_$VER.src.tar.bz2"
fi
if [[ -z "$SPOOLES_TGZ" ]]; then
  SPOOLES_TGZ="$BUILD/spooles.2.2.tgz"
  [[ -f "$SPOOLES_TGZ" ]] || curl -fL -o "$SPOOLES_TGZ" "https://www.netlib.org/linalg/spooles/spooles.2.2.tgz"
fi

rm -rf CalculiX
tar xjf "$CCX_TGZ"                       # -> CalculiX/ccx_<ver>/src
SRC="$BUILD/CalculiX/ccx_$VER/src"
mkdir -p CalculiX/spooles
SPOOLES_LIST=$(tar tzf "$SPOOLES_TGZ")       # (no pipe into grep -q: SIGPIPE + pipefail)
if grep -qx "Make.inc" <<<"$SPOOLES_LIST"; then
  tar xzf "$SPOOLES_TGZ" -C CalculiX/spooles
else                                     # tarball with a top-level folder
  tar xzf "$SPOOLES_TGZ" -C CalculiX/spooles --strip-components=1
fi

# 2. SPOOLES sparse solver (same patches as the Homebrew formula), static library
SP="$BUILD/CalculiX/spooles"
[[ -f "$SP/Make.inc" ]] || { echo "unexpected SPOOLES tarball layout" >&2; exit 1; }
CC_BIN=$(command -v cc)
perl -pi -e "s#/usr/lang-4.0/bin/cc#$CC_BIN#" "$SP/Make.inc"
perl -pi -e 's/drawTree\.c/tree.c/' "$SP/Tree/src/makeGlobalLib"
perl -pi -e 's/IVinit\(nfront, NULL\)/IVinit(nfront, 0)/' "$SP/ETree/src/transform.c"
make -s -C "$SP" lib >/dev/null 2>&1
make -s -C "$SP/MT/src" makeLib >/dev/null 2>&1
echo "SPOOLES built"

# 3. ccx patches from the formula: drop the bundled-library dependency, fix a void return
perl -pi -e "s/^(ccx_$VER: \\\$\\(OCCXMAIN\\) ccx_$VER\\.a)\\s+\\\$\\(LIBS\\)/\$1/" "$SRC/Makefile"
perl -0pi -e 's/(\*iprfnp=iprfn;\*konrfnp=konrfn;\*ratiorfnp=ratiorfn;\s*)return NULL;/$1return;/' "$SRC/readnewmesh.c"

# 4. Put the UMAT in place of the empty template
cp "$UMAT" "$SRC/umat.f"
printf '      IMPLICIT REAL*8(A-H,O-Z)\n' > "$SRC/ABA_PARAM.INC"

# 5. Build ccx
if [[ "$(uname)" == "Darwin" ]]; then
  ARPACK=$(pkg-config --libs arpack)
  LIBS="../../spooles/spooles.a $ARPACK -framework Accelerate"
else
  ARPACK=$(pkg-config --libs arpack 2>/dev/null || echo "-larpack")
  LIBS="../../spooles/spooles.a $ARPACK -llapack -lblas -lpthread"
fi
make -s -j"$JOBS" -C "$SRC" "ccx_$VER" \
  CC=cc FC="$FC" \
  CFLAGS="-O2 -w -I../../spooles -DARCH=Linux -DSPOOLES -DARPACK -DMATRIXSTORAGE -DUSE_MT=1" \
  FFLAGS="-O2 -w -fopenmp -cpp -DCCX_COROT -fallow-argument-mismatch" \
  LIBS="$LIBS" >"$BUILD/make.log" 2>&1 || { tail -30 "$BUILD/make.log"; exit 1; }

mkdir -p "$HERE/bin"
cp "$SRC/ccx_$VER" "$HERE/bin/ccx_${VER}_umat"
echo "Built $HERE/bin/ccx_${VER}_umat"

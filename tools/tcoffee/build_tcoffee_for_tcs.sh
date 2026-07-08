#!/usr/bin/env bash
#
# build_tcoffee_for_tcs.sh
# -----------------------------------------------------------------------------
# Build a T-Coffee `t_coffee` binary that works with FunVIP's optional TCS
# (Transitive Consistency Score) alignment-validation step.
#
# Why this script exists
# ----------------------
# The prebuilt (e.g. bioconda) t-coffee indexes an internal table by the raw OS
# process id, but the table is sized by a COMPILE-TIME constant MAX_N_PID=260000.
# Modern Linux kernels default kernel.pid_max=4194304, so t-coffee's own PID
# overflows that table on essentially every run -> out-of-bounds access ->
# SIGSEGV -> its signal handler re-runs the faulting instruction forever, growing
# the heap until all system memory is gone. It is a crash-loop (not a real leak)
# and it fires on ANY invocation, even `t_coffee -version`. The MAX_N_PID_4_TCOFFEE
# environment variable does NOT fully fix it: even where the runtime check honours
# it, the table is still allocated with the compile-time MAX_N_PID and overflows.
#
# The only correct fix is to recompile t-coffee with MAX_N_PID >= kernel.pid_max.
# This script downloads the T-Coffee source, applies exactly that one-line change
# (tools/tcoffee/max_n_pid.patch), builds only the `t_coffee` binary, verifies it,
# and installs it.
#
# FunVIP does NOT ship a t-coffee binary: t-coffee is GPL but its authors add a
# "non-commercial pipelines only" caveat (see README.md / the printed notice), so
# redistributing the binary is avoided. You build it yourself with this script and
# thereby accept T-Coffee's license directly.
#
# Usage
# -----
#   bash build_tcoffee_for_tcs.sh [--prefix DIR] [--max-n-pid N] [--jobs N] [--keep]
#
#   --prefix DIR    install t_coffee into DIR/bin (default: $CONDA_PREFIX, else
#                   ~/.local). DIR/bin must be on your PATH for FunVIP to find it.
#   --max-n-pid N   override the compiled limit (default: kernel.pid_max + margin,
#                   floored at 4194305). Increase only if your kernel.pid_max is
#                   unusually large.
#   --jobs N        parallel compile jobs (default: all cores).
#   --keep          keep the build directory (default: removed on success).
# -----------------------------------------------------------------------------
set -euo pipefail

# --- pinned upstream source (matches bioconda t-coffee 13.46.2) ---------------
TC_VERSION="13.46.2.7c9e712d"
TC_URL="https://s3.eu-central-1.amazonaws.com/tcoffee-packages/Archives/T-COFFEE_distribution_Version_${TC_VERSION}.tar.gz"
TC_SHA256="84f9b4076767d39dec6619c5eb91c9538a7c58c68a3731a92ebbf2e1f914296f"

# --- defaults -----------------------------------------------------------------
PREFIX="${CONDA_PREFIX:-$HOME/.local}"
JOBS="$( (command -v nproc >/dev/null && nproc) || echo 4)"
KEEP=0
# kernel.pid_max + margin, floored at PID_MAX_LIMIT+1 (4194304+1) so the table can
# hold any 64-bit-Linux PID even if this build host's pid_max is lower than a
# future run host's.
_PIDMAX="$(cat /proc/sys/kernel/pid_max 2>/dev/null || echo 4194304)"
MAX_N_PID=$(( _PIDMAX + 16 )); [ "$MAX_N_PID" -lt 4194305 ] && MAX_N_PID=4194305

while [ $# -gt 0 ]; do
  case "$1" in
    --prefix)    PREFIX="$2"; shift 2 ;;
    --max-n-pid) MAX_N_PID="$2"; shift 2 ;;
    --jobs)      JOBS="$2"; shift 2 ;;
    --keep)      KEEP=1; shift ;;
    -h|--help)   sed -n '2,45p' "$0"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

say() { printf '\n\033[1m==> %s\033[0m\n' "$*"; }
die() { printf '\nERROR: %s\n' "$*" >&2; exit 1; }

# --- prerequisites ------------------------------------------------------------
say "Checking build tools"
MISSING=""
for t in make g++ sed awk tar; do command -v "$t" >/dev/null || MISSING="$MISSING $t"; done
command -v curl >/dev/null || command -v wget >/dev/null || MISSING="$MISSING curl-or-wget"
[ -n "$MISSING" ] && die "missing required tools:$MISSING (install e.g. 'conda install -c conda-forge gxx make', or a build-essential package)"
echo "  ok (g++ $(g++ -dumpversion), make, $( command -v curl >/dev/null && echo curl || echo wget ))"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/tcoffee_build.XXXXXX")"
cleanup() { [ "$KEEP" -eq 1 ] || rm -rf "$WORK"; }
trap cleanup EXIT

# --- download + verify --------------------------------------------------------
say "Downloading T-Coffee ${TC_VERSION} source"
TARBALL="$WORK/tcoffee_src.tar.gz"
if command -v curl >/dev/null; then curl -fSL -o "$TARBALL" "$TC_URL"
else wget -O "$TARBALL" "$TC_URL"; fi

say "Verifying checksum"
if command -v sha256sum >/dev/null; then
  echo "${TC_SHA256}  ${TARBALL}" | sha256sum -c - || die "checksum mismatch (download corrupt or upstream changed)"
elif command -v shasum >/dev/null; then
  echo "${TC_SHA256}  ${TARBALL}" | shasum -a 256 -c - || die "checksum mismatch"
else
  echo "  (no sha256 tool found; skipping verification)"
fi

tar -xzf "$TARBALL" -C "$WORK"
SRCROOT="$(find "$WORK" -maxdepth 1 -type d -name 'T-COFFEE_distribution_*' | head -1)"
[ -d "$SRCROOT/t_coffee_source" ] || die "unexpected source layout under $SRCROOT"

# --- patch MAX_N_PID ----------------------------------------------------------
say "Patching MAX_N_PID 260000 -> ${MAX_N_PID} (kernel.pid_max=${_PIDMAX})"
DEFS="$SRCROOT/t_coffee_source/coffee_defines.h"
grep -qE '^#define MAX_N_PID[[:space:]]+260000' "$DEFS" \
  || die "MAX_N_PID definition not found/changed upstream in $DEFS; update this script"
sed -i.orig -E "s/^#define MAX_N_PID[[:space:]]+260000/#define MAX_N_PID       ${MAX_N_PID}/" "$DEFS"
grep -E '^#define MAX_N_PID' "$DEFS" | sed 's/^/  /'

# --- build only the t_coffee binary (modern-g++ compatible flags) -------------
say "Building t_coffee (-j${JOBS}); this takes ~1 minute"
( cd "$SRCROOT/t_coffee_source"
  make t_coffee -j"$JOBS" CC=g++ \
    CFLAGS="-O3 -fpermissive -fsigned-char -Wno-write-strings -Wno-register -Wno-return-type -Dregister=" )
BIN="$SRCROOT/t_coffee_source/t_coffee"
[ -x "$BIN" ] || die "build did not produce a t_coffee binary"

# --- self-test: must NOT crash-loop and MUST produce a score ------------------
say "Self-test: running TCS on a tiny alignment"
T="$WORK/selftest"; mkdir -p "$T"
cat > "$T/aln.fasta" <<'FASTA'
>s1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>s2
ACGTACGTACGTACGTTCGTACGTACGTACGTACGTACGT
>s3
ACGTACGTACGTACGTACGTACGTACGAACGTACGTACGT
FASTA
( cd "$T"
  # 8 GB address-space guard so a still-broken build cannot take down the host
  ( ulimit -v 8388608 2>/dev/null || true
    timeout 120 "$BIN" -infile aln.fasta -method fast_pair -type DNA -evaluate \
      -output score_ascii -outfile out.tcs -quiet </dev/null >tc.log 2>&1 ) || true
  if grep -q "MAX_N_PID exceded" tc.log; then
    die "self-test still hit MAX_N_PID -- increase --max-n-pid above your kernel.pid_max"
  fi
  [ -s out.tcs ] || { echo "----- t_coffee output -----"; tail -20 tc.log; die "self-test produced no .tcs score file"; }
  echo "  ok: produced a valid TCS score file"
)

# --- install ------------------------------------------------------------------
say "Installing to ${PREFIX}/bin/t_coffee"
mkdir -p "$PREFIX/bin"
install -m 0755 "$BIN" "$PREFIX/bin/t_coffee"

cat <<EOF

Done. Installed: ${PREFIX}/bin/t_coffee  (MAX_N_PID=${MAX_N_PID})

  - Ensure ${PREFIX}/bin is on your PATH so FunVIP finds it:
      export PATH="${PREFIX}/bin:\$PATH"
  - FunVIP auto-enables TCS when t_coffee is on PATH (disable with --notcs).
  - Quick check:  t_coffee -version   (should print a version, not hang)

$(sed -n '/^T-COFFEE LICENSE NOTICE/,/^END NOTICE/p' "$SCRIPT_DIR/README.md" 2>/dev/null || true)
EOF

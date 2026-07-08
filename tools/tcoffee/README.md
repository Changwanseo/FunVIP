# Building T-Coffee for FunVIP's optional TCS step

FunVIP can validate its alignments with **TCS** (Transitive Consistency Score,
from T-Coffee). TCS is **optional** and off unless a `t_coffee` binary is on your
PATH. FunVIP itself works fully without it; only enable TCS if you specifically
want that alignment-quality metric (for example to reproduce a published TCS
analysis).

FunVIP does **not** bundle a `t_coffee` binary. The prebuilt (bioconda) t-coffee
is broken on modern Linux, and T-Coffee's license discourages redistribution
(see the notice below). This directory instead ships a **build recipe** so you can
compile a working `t_coffee` yourself.

## Why a prebuilt t-coffee breaks (and eats all your RAM)

The stock t-coffee indexes an internal table by the **raw OS process id**, but the
table is sized by a **compile-time constant `MAX_N_PID = 260000`**. Modern Linux
kernels default `kernel.pid_max = 4194304`, so t-coffee's own PID overflows the
table on essentially every run:

```
out-of-bounds table access  ->  SIGSEGV  ->  t-coffee's signal handler re-runs the
faulting instruction forever  ->  heap grows without bound  ->  all RAM consumed
```

It is a **crash-loop, not a memory leak**, and it triggers on *any* invocation,
even `t_coffee -version`. The `MAX_N_PID_4_TCOFFEE` environment variable does not
fully fix it: even when the runtime check honours it, the table is still allocated
with the compile-time `MAX_N_PID` and overflows anyway. The only correct fix is to
**recompile with `MAX_N_PID >= kernel.pid_max`**, which is exactly what
[`max_n_pid.patch`](max_n_pid.patch) does (a one-line change to `coffee_defines.h`).

## Usage

```bash
bash tools/tcoffee/build_tcoffee_for_tcs.sh            # installs into $CONDA_PREFIX/bin
# or
bash tools/tcoffee/build_tcoffee_for_tcs.sh --prefix ~/tcoffee   # -> ~/tcoffee/bin
```

Requirements: a C++ compiler (`g++`), `make`, and `curl` or `wget`. In a conda
env: `conda install -c conda-forge gxx make`.

The script downloads the pinned T-Coffee source (checksum-verified), applies the
`MAX_N_PID` patch, builds only the `t_coffee` binary, **self-tests it** (runs TCS
on a tiny alignment and fails if it still crashes), and installs it. Make sure the
install `bin` directory is on your PATH; FunVIP then enables TCS automatically
(disable with `--notcs`).

If your kernel's `pid_max` is unusually large, pass `--max-n-pid N` with a value
above it.

## Alternatives (no compiling)

- Run FunVIP **without TCS** (the default when t-coffee is absent) — nothing to do.
- On a host where you have root, lowering `kernel.pid_max` below 260000 also avoids
  the crash, but it is system-wide and affects everything else, so it is not
  recommended on shared machines.

---

T-COFFEE LICENSE NOTICE
  T-Coffee is developed by Cedric Notredame (CNRS) and distributed under the GNU
  General Public License, with the authors' added condition that it may be
  incorporated only into NON-COMMERCIAL pipelines; for commercial use, contact the
  T-Coffee authors (https://tcoffee.org). By building and using t_coffee via this
  script you obtain and accept T-Coffee directly under its own license -- FunVIP
  neither relicenses nor redistributes the T-Coffee binary. See:
  https://tcoffee.readthedocs.io/en/latest/tcoffee_license.html
  Please also cite the T-Coffee / TCS papers when you use this functionality.
END NOTICE

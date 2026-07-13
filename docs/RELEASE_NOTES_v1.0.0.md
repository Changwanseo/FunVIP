# FunVIP 1.0.0

First stable release. FunVIP moved from ete3 to ete4, became substantially faster,
was hardened across the board, and now installs natively on Windows with the same
one-line `pip install` as Linux and macOS.

## Highlights

- **New ete4 tree engine.** The whole tree layer was ported from ete3 to ete4.
  Supported Python is now 3.9-3.13.
- **Native Windows support.** ete4 has no Windows wheel on PyPI, so FunVIP bundles
  prebuilt ete4 wheels and installs the right one on first run (building from a
  patched source only if no wheel matches). Windows now uses the same install steps
  as the other platforms.
- **About 13x faster tree interpretation**, plus faster clustering, dataset
  building, and concatenation.
- **Conserved-gene (5.8S) crash fixed** and a **result-affecting cluster-cutoff bug
  fixed** (an explicit cutoff used to be discarded).
- **TCS no longer exhausts system memory.** The optional t-coffee step is now safely
  contained and skipped when unavailable.
- **Test suite + CI, reproducible installs** (`environment.yml`, `Dockerfile`), and
  a **self-contained HTML report**.

## Install

Linux (external tools via conda, FunVIP + ete4 via pip):
```
conda create -n FunVIP python=3.12
conda activate FunVIP
conda config --add channels conda-forge
conda install -c bioconda raxml iqtree "modeltest-ng==0.1.7" mmseqs2 "blast>=2.12" mafft trimal gblocks fasttree
pip install FunVIP
FunVIP --test terrei --email <your email>
```

Windows is the same recipe (`conda create` / `activate` / `pip install FunVIP` /
`--test`); FunVIP installs the bundled ete4 automatically on first run. See
`tutorial/installation.md` for macOS and the one-file `environment.yml` recipe.

## Upgrading from 0.5.x

0.5.x is built on **ete3** and 1.0 on **ete4** (different packages, not two versions
of one). Do **not** `pip install --upgrade` across this jump; rebuild the conda
environment from scratch instead. Your input data, databases, and result folders are
untouched. Full steps are in `tutorial/installation.md`.

## Changes

### Performance
- Tree interpretation about 13x faster (memoized tree search, per-tree set-cache in
  `decide_type`, numpy `calculate_zero`, within-cluster `diff_min` instead of an
  O(n^2) distance matrix); removed an O(depth^2) copy in `seperate_clade`.
- Vectorized group clustering; `generate_dataset` pre-indexes sequences by group.
- Pool workers share the sequence universe via fork copy-on-write instead of
  per-task pickling.
- `concatenate` frees per-gene tables and dictionary-encodes string columns.

### Fixed
- Conserved-gene (5.8S) interpretation crash (silent RecursionError / pickle
  failure) now interprets correctly.
- Cluster cutoff given explicitly (preset or CLI) is now applied instead of being
  discarded.
- TCS (t-coffee) no longer hangs and consumes all memory: FunVIP never executes
  t-coffee just to detect it, caps and times out the TCS run, and skips it cleanly
  on failure. See `tools/tcoffee/` for a recipe to build a working t-coffee.
- Many correctness and robustness fixes: external-tool exit codes are checked,
  option/preset validation bugs, several always-wrong guards, a list-mutation bug,
  `--all` handling, `homogenize`, and GenMine executable resolution.

### Added
- Typed exception hierarchy with a top-level handler and a logging overhaul.
- Method-aware external-tool preflight and a more robust version probe.
- Unit-test suite and CI (pip matrix on Python 3.10-3.13 plus a conda install smoke
  test) and a marked terrei end-to-end integration test.
- `environment.yml` and a `Dockerfile`.
- Self-contained HTML report and per-group statistics.
- Parameter reference (`docs/parameters.md`), refreshed README and installation
  docs, and build recipes under `tools/` (working t-coffee for TCS; ete4 Windows
  wheels).

### Changed (behavior)
- Cluster e-value default is now `1e-4`.
- Group-assignment ties are broken by top-bitscore majority.
- Missing per-gene bitscores in the concatenated table are filled by projecting onto
  the fitted regression line.
- Removed confirmed-dead / broken code paths.

## Known issues

- TCS is optional and Linux-only. The prebuilt bioconda t-coffee is broken on modern
  kernels (see `tools/tcoffee/`); FunVIP skips TCS automatically when t-coffee is
  absent.
- On Windows, MAFFT's bundled shell prints a harmless "could not find /tmp" warning;
  it does not affect results.

## Citation

If you use FunVIP, please cite the FunVIP publication (see the repository README).

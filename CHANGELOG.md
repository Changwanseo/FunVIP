# Changelog

All notable changes to FunVIP are documented here. This project adheres to
[Semantic Versioning](https://semver.org/).

## [1.0.1] - 2026-09-10

### Fixed
- **Crash while cleaning up empty alignments.** A sequence left with no residues
  after trimming was dropped from its dataset with a warning that looked its id up
  in `dict_hash_id`, a dictionary that is never populated, so the run died with
  `KeyError: 'HS<n>HE'` instead of continuing. The warning now reads the id from the
  removed sequence itself.
- **UnicodeEncodeError when decoding hashes on a non-UTF-8 Windows locale.**
  `hasher.decode` read and wrote with the platform default encoding, so on a cp949
  (Korean) console any non-ASCII character in a sequence name (an accented author
  name, for instance) aborted the decode of trees, SVGs and alignments. Both ends
  are now explicitly UTF-8.

## [1.0.0] - 2026-07-09

First stable release. The pipeline was ported from ete3 to ete4, made
substantially faster, hardened across the board, and given native Windows support.

Upgrading from 0.5.x: 0.5.x is built on ete3 and 1.0 on ete4 (different packages),
so rebuild the conda environment instead of `pip install --upgrade`. See
`tutorial/installation.md` (Migrating from 0.5.x to 1.0).

### Major
- **ete4 engine.** Migrated the whole tree layer from ete3 to ete4. Supported
  Python is now 3.9-3.13.
- **Native Windows support.** ete4 has no Windows wheel on PyPI; FunVIP now bundles
  prebuilt ete4 wheels (`funvip/_vendor/ete4_wheels`) and installs the matching one
  on first run, falling back to building ete4 from source (patched for Windows) when
  no wheel matches. The Windows install is now the same `pip install FunVIP` as
  Linux/macOS.

### Performance
- Tree interpretation roughly 13x faster (memoized tree search, per-tree set-cache
  in `decide_type`, numpy `calculate_zero`, within-cluster `diff_min` instead of an
  O(n^2) distance matrix).
- `seperate_clade` no longer re-deep-copies the resolved comb at every level
  (dropped an O(depth^2) copy).
- Vectorized group clustering; `generate_dataset` pre-indexes sequences by group
  (was O(groups x genes x N)).
- Pool workers share the sequence universe via fork copy-on-write instead of
  per-task pickling.
- `concatenate` frees per-gene search tables after building the concatenated table
  and dictionary-encodes its string columns.

### Fixed
- **5.8S / conserved-gene interpretation crash.** The conserved-gene tree step used
  to fail silently (RecursionError / pickle failure); it now interprets correctly.
- **Result-affecting preset bug.** An explicitly given cluster cutoff was discarded,
  so runs used the wrong cutoff; presets and CLI cutoffs are now applied correctly.
- **TCS (t-coffee) no longer eats all system memory.** The bundled t-coffee
  crash-loops on modern kernels (PID-indexed static array vs. large
  `kernel.pid_max`); FunVIP no longer executes t-coffee just to detect it, caps and
  times out the TCS run, and skips TCS cleanly on failure. See `tools/tcoffee/` for
  a recipe to build a working t-coffee.
- Numerous correctness and robustness fixes: external-tool exit codes are checked,
  option/preset validation bugs, several always-wrong guards, a list-mutation bug,
  `--all` handling, `homogenize`, and the GenMine executable resolution.

### Added
- Typed exception hierarchy with a top-level handler, plus a logging overhaul across
  the core modules.
- Method-aware external-tool preflight (checks only the tools the chosen methods
  need) and a more robust `Version` probe.
- Unit-test suite and CI (pip matrix on Python 3.10-3.13 + a conda install smoke
  test), and a marked terrei end-to-end integration test.
- `environment.yml` and a `Dockerfile` for reproducible installs.
- Self-contained HTML report and per-group statistics.
- Parameter reference (`docs/parameters.md`), refreshed README and installation
  docs, and build recipes under `tools/` (working t-coffee for TCS; ete4 Windows
  wheels).

### Changed (behavior)
- Cluster e-value default is now `1e-4` (a single clean preset key).
- Group-assignment ties are broken by top-bitscore majority.
- Missing per-gene bitscores in the concatenated table are filled by projecting onto
  the fitted regression line.
- Removed confirmed-dead / broken code paths.

[1.0.1]: https://github.com/Changwanseo/FunVIP/releases/tag/v1.0.1
[1.0.0]: https://github.com/Changwanseo/FunVIP/releases/tag/v1.0.0

# Changelog

All notable changes to FunVIP are documented here. This project adheres to
[Semantic Versioning](https://semver.org/).

## [1.0.2] - 2026-09-30

Faster outgroup selection and tree drawing on large databases, earlier input
validation, and a fix for runs at the default verbosity on Windows and macOS.
Analysis results are unchanged.

### Performance
- **Outgroup selection groups the search table once.** Preparing the outgroup tasks
  built a new `groupby` over the whole concatenated search table for every group,
  which is O(groups x rows) and ran in the main process before the worker pool
  started. On a 121,133-sequence, 5,266-genus ITS database this took about 2 h with no
  log output. The table is now grouped once. At 1,500 groups and 300,000 search hits,
  preparing the tasks went from 34.8 s to 0.15 s with the same task list.
- **Single-gene groups interpret and draw each tree once.** With one gene, the
  concatenated tree is a copy of the gene tree, but it was interpreted and rendered a
  second time. When the tree, the trimmed alignment, the partition and the query and
  outgroup lists of both datasets are identical, the concatenated tree now reuses the
  gene tree's interpretation, SVG and report rows. Otherwise it is interpreted as
  before. On the 121,133-sequence ITS database about half of the visualize step was
  this duplicate work.

### Fixed
- **Crash on Windows and macOS at the default verbosity.** Outgroup selection and tree
  interpretation handed the sequence records to their worker processes through a
  module-level store that only fork-started workers inherit. Windows and macOS start
  workers by spawning, so every run below `--verbose 3` stopped at the outgroup step
  with `KeyError: 'list_FI'`. Each worker now receives the records when it starts.
  The bundled test presets use `--verbose 3`, which is why `--test` runs passed.
- **Free text in a gene column no longer reaches GenMine.** A cell was treated as an
  accession whenever an accession appeared anywhere in it, and the whole cell was
  sent to GenMine, which failed on free text and accession lists. A cell is now
  downloaded only when it is a single token. Other cells containing an accession stop
  the run at input validation, reporting the table, line and column. FASTA-formatted
  cells (starting with `>`) are read as sequences even when the header holds an
  accession.
- **IDs that collide after non-ASCII replacement stop at input validation.**
  `LÖ21-04` and `LO21-04`, or IDs that differ only in a Unicode hyphen, became the
  same ID and were merged into one record, or failed later as a colliding genus or
  datatype. Each such pair is now reported with both original strings, across DB and
  query tables. An ID repeated verbatim is still merged with a warning, as before.
- IDs without any letter or digit (`--`, `...`, a lone dash) are rejected as empty
  IDs.

### Documentation
- The `--outgroupoffset` help and `docs/parameters.md` now state that the value also
  discards every search hit at or below it before clustering, besides setting the
  bitscore gap between ingroup and outgroup.

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

[1.0.2]: https://github.com/Changwanseo/FunVIP/releases/tag/v1.0.2
[1.0.1]: https://github.com/Changwanseo/FunVIP/releases/tag/v1.0.1
[1.0.0]: https://github.com/Changwanseo/FunVIP/releases/tag/v1.0.0

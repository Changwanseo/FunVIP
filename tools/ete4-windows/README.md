# Windows wheels for ete4

FunVIP depends on **ete4**, which upstream does **not** ship for Windows:

- PyPI has only an **sdist** (source) for ete4 — so `pip install` tries to compile,
  and that compile fails on Windows.
- conda-forge builds ete4 only for **linux-64 / osx-64 / osx-arm64** (no win-64);
  bioconda has no ete4 at all.

Because ete4 is a hard FunVIP dependency, this means FunVIP currently **cannot be
installed with pip or conda on Windows** — only WSL/Docker work. This directory
produces a prebuilt Windows wheel so ordinary Windows users can install FunVIP
without a compiler or a terminal-heavy workflow.

## The fix

The only change needed is the small path-separator fix from
[etetoolkit/ete PR #783](https://github.com/etetoolkit/ete/pull/783), captured here
as [`ete4-windows.patch`](ete4-windows.patch) (three lines):

- `setup.py`: derive Cython module names with `os.path.sep` instead of a hardcoded
  `'/'` (Windows uses `\`), so the extensions are named/placed correctly.
- `ete4/config.py`: use `os.path.expanduser('~')` instead of `os.environ['HOME']`
  (Windows has no `HOME`).

The patch is a no-op on Linux/macOS (verified: the patched source still builds a
normal Linux wheel), so it is safe to apply unconditionally.

## Build the wheels (GitHub Actions)

1. Push this repo to GitHub (the workflow lives in
   `.github/workflows/build-ete4-windows-wheels.yml`).
2. Actions tab → **build-ete4-windows-wheels** → **Run workflow** (optionally set the
   ete4 version; default `4.4.0`).
3. When it finishes, download the **`ete4-<version>-windows-wheels`** artifact. It
   contains `ete4-<version>-cp3XX-cp3XX-win_amd64.whl` for CPython 3.10–3.13.

The workflow downloads the pinned ete4 sdist, applies the patch, builds with
`cibuildwheel` on `windows-latest`, and runs an import test that loads the compiled
extension (so a green run means the wheel actually *works* on Windows, not just that
it compiled). To make a first "does it even compile?" run faster, narrow
`CIBW_BUILD` in the workflow to a single version, e.g. `cp312-win_amd64`.

## Using the wheels

Host the wheels where users can reach them (e.g. attach them to a **GitHub Release**),
then a Windows install becomes, for example:

```
pip install ete4 --only-binary :all: --find-links https://github.com/<you>/<repo>/releases/download/<tag>/
pip install funvip
```

(or wrap that in a one-line `.bat` so users do not type anything). FunVIP itself
already bundles the Windows binaries of the other external tools under
`funvip/external/`, so ete4 is the last missing piece for a Windows install.

## Even cleaner: upstream the patch

The best outcome is to get PR #783 merged and ask the ete4 maintainers to publish
official `win_amd64` wheels. Then you redistribute nothing and Windows users install
ete4 straight from PyPI. Hosting your own wheels (below) is a fine stopgap until then.

---

ETE4 LICENSE NOTICE
  ete4 is distributed under the GNU General Public License v3 or later
  (GPL-3.0-or-later) -- a standard GPL with no extra restrictions. You may build and
  redistribute these wheels, provided you (1) offer the corresponding source for the
  exact version built, (2) ship this patch and note that the source is modified, and
  (3) keep ete4's LICENSE and copyright notices (the wheel already carries them). No
  additional restrictions may be added. See: https://github.com/etetoolkit/ete
END NOTICE

# Releasing FunVIP

Maintainer checklist for cutting a release (e.g. 1.0.0). Run from a clean checkout
of `development` (or the release branch).

## 1. Bundle the ete4 Windows wheels
ete4 has no Windows wheel on PyPI, so FunVIP ships them under
`funvip/_vendor/ete4_wheels/` and installs the right one on first Windows run.

- [ ] Run the **build-ete4-windows-wheels** GitHub Action (Actions tab -> Run
      workflow). It builds `cp310`-`cp313` `win_amd64` wheels and import-tests them.
- [ ] Download the artifact and place all four `.whl` files in
      `funvip/_vendor/ete4_wheels/`, then commit them.
- [ ] Confirm the wheel version satisfies the ete4 pin in `pyproject.toml`
      (`>=4.4.0,<4.5.0`).

(The `cp312` wheel is already vendored and validated on a real Windows machine;
the Action produces the remaining versions on a clean runner.)

## 2. Version and changelog
- [ ] Bump `version` in `pyproject.toml`. It is the single source of truth
      (`FunVIP --version` reads it via `importlib.metadata`).
- [ ] Add/finish the release entry in `CHANGELOG.md` and set its date.

## 3. Verify (do not skip)
- [ ] `pytest` passes (unit suite).
- [ ] Linux end-to-end: `FUNVIP_TEST_EMAIL=<you> pytest -m integration` (or
      `FunVIP --test terrei --email <you>` in a temp dir) produces a result CSV.
- [ ] Windows end-to-end: in a fresh Python 3.12 conda env, `pip install .` then
      `FunVIP --test terrei --email <you>` -- confirm it prints
      "installing bundled ete4" (uses the wheel, no compile) and completes.
- [ ] macOS: run the recipe once, or mark Intel Mac "experimental" in the docs.
- [ ] At least one real, dataset-scale run (the intended workload) completes.

## 4. Build and publish
- [ ] Clean build: `python -m build` (produces sdist + a `py3-none-any` wheel;
      confirm `funvip/_vendor/ete4_wheels/*.whl`, `funvip/external/*`, `funvip/data`,
      and `funvip/preset` are inside the built artifacts).
- [ ] `twine check dist/*`.
- [ ] Upload: `twine upload dist/*` (test first on TestPyPI if unsure). Note: the
      PyPI `FunVIP` was last published at 0.5.x (ete3); 1.0.0 is the ete4 line.
- [ ] Do NOT put a direct git/URL ete4 dependency in `pyproject.toml`; PyPI rejects
      those. The Windows wheels ship inside the package, not as a URL dependency.

## 5. Tag and announce
- [ ] Tag: `git tag v1.0.0 && git push origin v1.0.0`.
- [ ] Create the GitHub release from the tag; paste the CHANGELOG entry.
- [ ] Merge `development` into `main` if that is the release branch.

## Known issues to mention in the release notes
- TCS is optional and Linux-only; the prebuilt bioconda t-coffee is broken (see
  `tools/tcoffee/`). FunVIP skips TCS automatically when t-coffee is absent.
- On Windows, MAFFT's bundled bash prints a harmless "could not find /tmp" warning;
  it does not affect results.

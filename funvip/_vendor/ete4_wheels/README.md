# Bundled ete4 Windows wheels

FunVIP installs the matching wheel from this directory on the **first Windows run**
(see `_ensure_ete4` in `funvip/main.py`), because ete4 has no Windows wheel on PyPI
and cannot compile there. On Linux/macOS this directory is unused (ete4 installs
normally from PyPI).

## What goes here

`ete4-<version>-cp3XX-cp3XX-win_amd64.whl` for each supported CPython (3.10–3.13),
e.g. `ete4-4.4.0-cp312-cp312-win_amd64.whl`.

## How to produce them

Run the **build-ete4-windows-wheels** GitHub Actions workflow (see
`tools/ete4-windows/`), download the `ete4-<version>-windows-wheels` artifact, and
drop the `.whl` files here, then commit them. They are declared as package data in
`pyproject.toml`, so they ship inside the FunVIP wheel/sdist.

The wheel version must satisfy FunVIP's ete4 pin (`>=4.4.0,<4.5.0`).

## License

These wheels are ete4, distributed under **GPL-3.0-or-later** (see `LICENSE` in this
directory). FunVIP bundles them unmodified except for the small Windows build fix in
`tools/ete4-windows/ete4-windows.patch`; the corresponding source is ete4's public
release plus that patch. FunVIP (also GPL-3.0) may redistribute them under the GPL.

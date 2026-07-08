## Installation

FunVIP needs two things: the Python package (installed with `pip`, which pulls in
ete4 and the other Python dependencies) and a set of external command-line tools
(installed with `conda` from the bioconda / conda-forge channels). ete4 installs
cleanly from PyPI on **Linux and macOS**; on **Windows** use conda / WSL / Docker
(see the Windows section), because ete4's Cython build does not install with a
bare `pip install` there.

The bundled test run `FunVIP --test terrei --email <your email>` at the end of
each recipe checks the installation. `--test` is case-insensitive, so `terrei`
and `Terrei` both work.

<br><br/>

### Recommended: one conda environment file
From the repository root (or after downloading `environment.yml`):
```
conda env create -f environment.yml     # or: mamba env create -f environment.yml
conda activate funvip
FunVIP --test terrei --email <your email>
```
This installs the external tools and FunVIP (from PyPI) together. To use the
development version, install from source (below) instead.

<br><br/>

### Linux
1. ```conda create -n FunVIP python=3.12```
2. ```conda activate FunVIP```
3. ```conda config --add channels conda-forge```
4. ```conda install -c bioconda raxml iqtree "modeltest-ng==0.1.7" mmseqs2 "blast>=2.12" mafft trimal gblocks fasttree "t-coffee>=13"```
5. ```pip install FunVIP```
6. run ```FunVIP --test terrei --email <your email>``` to check installation

* t-coffee (used only for the optional TCS alignment-validation step) can be
  omitted; FunVIP detects that it is missing and skips TCS automatically. Omit it
  if you hit memory problems.
* For an Intel Mac this recipe may also work, but it is untested. Feedback welcome.

<br><br/>

### Apple Silicon Mac

#### Native installation (faster, Gblocks unavailable)
1. ```conda create -n FunVIP python=3.12```
2. ```conda activate FunVIP```
3. ```conda install conda-forge::pyqt bioconda::fasttree bioconda::raxml bioconda::iqtree bioconda::mmseqs2 bioconda::blast conda-forge::mafft bioconda::trimal```
4. ```pip install FunVIP```
5. run ```FunVIP --test terrei --email <your email>``` to check installation

#### Rosetta installation (slower, Gblocks available)
1. ```softwareupdate --install-rosetta```
2. ```CONDA_SUBDIR=osx-64 conda create -n FunVIP python=3.12```
3. ```conda activate FunVIP```
4. ```conda config --env --set subdir osx-64```
5. ```conda install pyqt openssl ca-certificates```
6. ```conda install -c bioconda raxml iqtree "mmseqs2<=16" "blast>=2.12" mafft trimal gblocks```
7. ```CONDA_SUBDIR=osx-arm64 conda install -c bioconda fasttree```
8. ```pip install FunVIP```
9. run ```FunVIP --test terrei --email <your email>``` to check installation

- Native installation is recommended; Gblocks can be substituted with trimAl.

<br><br/>

### Windows
A bare `pip install FunVIP` does not work on Windows because ete4's Cython build
fails to install. Use one of:

- **Docker (simplest):** from the repository root,
  ```
  docker build -t funvip .
  docker run --rm -v "%cd%:/data" funvip --test terrei --email <your email> --outdir /data/out
  ```
- **WSL2:** open an Ubuntu (WSL) shell and follow the Linux recipe above.

<br><br/>

### Installation from source (for developers)
1. ```git clone https://github.com/Changwanseo/FunVIP.git```
2. ```cd FunVIP```
3. ```conda create -n FunVIP python=3.12```
4. ```conda activate FunVIP```
5. ```conda config --add channels conda-forge```
6. ```conda install -c bioconda raxml iqtree "modeltest-ng==0.1.7" mmseqs2 "blast>=2.12" mafft trimal gblocks fasttree "t-coffee>=13"```
7. ```pip install -e ".[test]"```   (editable install + test dependencies; use ```pip install ./``` for a plain install)
8. run ```pytest``` for the unit tests, and ```FunVIP --test terrei --email <your email>``` for an end-to-end check

<br><br/>

### Upgrade FunVIP
```pip install FunVIP --upgrade```

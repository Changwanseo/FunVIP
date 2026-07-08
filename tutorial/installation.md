## Installation

FunVIP needs two things: the Python package (installed with `pip`, which pulls in
ete4 and the other Python dependencies) and a set of external command-line tools
(installed with `conda` from the bioconda / conda-forge channels). ete4 installs
from PyPI on **Linux and macOS**; **Windows has no PyPI ete4**, so FunVIP bundles a
prebuilt ete4 for Windows and installs it automatically on the first run (the
external tools are bundled for Windows too). The Windows steps are therefore the
same `pip install` as the other platforms.

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
4. ```conda install -c bioconda raxml iqtree "modeltest-ng==0.1.7" mmseqs2 "blast>=2.12" mafft trimal gblocks fasttree```
5. ```pip install FunVIP```
6. run ```FunVIP --test terrei --email <your email>``` to check installation

* **TCS (optional): do not `conda install t-coffee`.** t-coffee is deliberately
  left out of the command above. It is needed only for the optional TCS
  alignment-validation step, and the prebuilt bioconda t-coffee crash-loops and can
  consume all system memory on modern kernels (those with a large
  `kernel.pid_max`). FunVIP detects when t-coffee is absent and skips TCS
  automatically, so most users need to do nothing. If you specifically need TCS,
  build a working t-coffee with ```tools/tcoffee/build_tcoffee_for_tcs.sh```
  (see ```tools/tcoffee/README.md```).
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
Same as Linux, with nothing extra to install: the external tools are bundled for
Windows (in `funvip/external/`), and FunVIP installs a prebuilt **ete4** for you on
the first run (ete4 has no PyPI wheel for Windows). Use Python 3.10-3.13.
1. ```conda create -n FunVIP python=3.12```
2. ```conda activate FunVIP```
3. ```pip install FunVIP```
4. run ```FunVIP --test terrei --email <your email>``` to check installation

The first run prints `installing bundled ete4 ...` once and then continues. WSL or
Docker also work if you prefer: open an Ubuntu (WSL) shell and follow the Linux
recipe, or build the Docker image (`docker build -t funvip .`).

<br><br/>

### Installation from source (for developers)
1. ```git clone https://github.com/Changwanseo/FunVIP.git```
2. ```cd FunVIP```
3. ```conda create -n FunVIP python=3.12```
4. ```conda activate FunVIP```
5. ```conda config --add channels conda-forge```
6. ```conda install -c bioconda raxml iqtree "modeltest-ng==0.1.7" mmseqs2 "blast>=2.12" mafft trimal gblocks fasttree```   (t-coffee omitted on purpose; see the TCS note in the Linux section)
7. ```pip install -e ".[test]"```   (editable install + test dependencies; use ```pip install ./``` for a plain install)
8. run ```pytest``` for the unit tests, and ```FunVIP --test terrei --email <your email>``` for an end-to-end check

<br><br/>

### Upgrade FunVIP
```pip install FunVIP --upgrade```

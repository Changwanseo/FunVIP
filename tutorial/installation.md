## Installation

FunVIP needs two things: the Python package (installed with `pip`, which pulls in
ete4 and the other Python dependencies) and a set of external command-line tools
(installed with `conda` from the bioconda / conda-forge channels).

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

* TCS (optional alignment validation) is left out on purpose: the bioconda
  t-coffee can cause memory problems, and FunVIP skips TCS automatically when it is
  absent. To use TCS, build t-coffee with ```tools/tcoffee/build_tcoffee_for_tcs.sh```
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
1. ```conda create -n FunVIP python=3.12```
2. ```conda activate FunVIP```
3. ```pip install FunVIP```
4. run ```FunVIP --test terrei --email <your email>``` to check installation

* Docker and WSL2 also work (build the image with `docker build -t funvip .`, or
  open an Ubuntu WSL shell and follow the Linux recipe).

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

### Upgrade FunVIP (within the 1.x series)
Both sides use ete4, so a plain pip upgrade is fine:
```pip install FunVIP --upgrade```

### Migrating from 0.5.x to 1.0

**Do not `pip install --upgrade` across this jump.** FunVIP 0.5.x is built on
**ete3**; 1.0 switched to **ete4**. These are different packages, not two versions
of one, so a pip upgrade leaves the old `ete3` (and other stale 0.5.x dependencies)
behind alongside the new ones, which can conflict. Rebuild the conda environment
from scratch instead:

```
# 1. leave and delete the old environment (use your existing env's name)
conda deactivate
conda env remove -n FunVIP

# 2. recreate it fresh, pick the recipe for your platform above, e.g. Linux:
conda create -n FunVIP python=3.12
conda activate FunVIP
conda config --add channels conda-forge
conda install -c bioconda raxml iqtree "modeltest-ng==0.1.7" mmseqs2 "blast>=2.12" mafft trimal gblocks fasttree
pip install FunVIP
FunVIP --test terrei --email <your email>
```

Only the Python environment is rebuilt: your input data, databases, and result
folders are untouched. If you installed with the one-file recipe, the equivalent is
`conda env remove -n funvip` followed by `conda env create -f environment.yml`.

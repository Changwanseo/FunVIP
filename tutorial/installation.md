## Installation
### Windows
1. ```conda create -n FunVIP python=3.12```
2. ```conda activate FunVIP```
3. ```pip install FunVIP```
4. run ```FunVIP --test Terrei --email [your email] ``` to check installation
<br><br/>
## Linux
1. ```conda create -n FunVIP python=3.11```
2. ```conda activate FunVIP```
3. ```pip install FunVIP```
4. ```conda config --add channels conda-forge```
5. ```conda install -c bioconda raxml iqtree "modeltest-ng==0.1.7" mmseqs2 "blast>=2.12" mafft trimal gblocks fasttree "t-coffee>=13"```
6. run ```FunVIP --test Terrei --email [your email] ``` to check installation
* For intel mac system, this method may work, but we couldn't test it because we don't have any intel mac device. We're looking for feedbacks in intel mac
* If you have memory leakage problems, install without t-coffee 
<br><br/>

### Apple Silicon Mac
1. ```softwareupdate --install-rosetta```
2. ```CONDA_SUBDIR=osx-64 conda create -n FunVIP python=3.12```
3. ```conda activate FunVIP```
4. ```conda config --env --set subdir osx-64```
5. ```conda install pyqt```
6. ```pip install FunVIP```
7. ```conda install -c bioconda raxml iqtree "mmseqs2<=16" "blast>=2.12" mafft trimal gblocks fasttree```
8. run ```FunVIP --test Terrei --email [your email] ``` to check installation
<br><br/>
### Installation from source (For developers and core users)
* this is for developmental steps
1. ```git clone https://github.com/Changwanseo/FunVIP.git```
2. Move to ```~/FunVIP```
3. ```conda create -n FunVIP python=3.12```
4. ```conda activate FunVIP```
5. ```pip install ./```
6. run ```FunVIP --test Terrei --email [your email]``` to check installation
<br><br/>
### Upgrade FunVIP
``` pip install FunVIP --upgrade ```

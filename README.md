[![DOI](https://zenodo.org/badge/588465720.svg)](https://zenodo.org/doi/10.5281/zenodo.10714946)

#### This is Beta release. Bug reports are welcomed

## Scheduling
### Beta release part 1 (2023 Feburary ~ As paper published, ver 0.3)
- Will be tested by our lab memebers to fix bugs and advance features

### Beta release part 2 (As paper published ~ When pipeline gets stabled, ver 0.4)
- Will be tested by peer taxonomists

### Stable release (ver 1.0)


# FunVIP
"Fun"gal "V"alidation & "I"dentification "P"ipeline

An automatic tree-based sequence identification and validation pipeline for fungal species

- Automatic tree-based identification
- Works with multigene
- Data validation algorithm implemented


## See [tutorial](https://github.com/Changwanseo/FunVIP/blob/main/tutorial/tutorial.md) for step by step tutorial
## See [documentation](https://github.com/Changwanseo/FunVIP/blob/main/Documentation.md) for advanced usage


## Requirements
- Conda environment (See [https://www.anaconda.com/products/individual](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) for how to install conda environment)

<!--
## Installation with conda (May not work with Linux or Mac)
1. ```conda create -n FunVO{ python=3.10```
2. ```conda activate FunVIP```
3. ```conda install -c cwseo FunVIP```
4. run ```FunVIP --test Terrei --email [your email] ``` to check installation
If this one fails, use next one
-->

## Installation
### Windows
1. Install visual c++ [here](https://visualstudio.microsoft.com/visual-cpp-build-tools/)
2. ```conda create -n FunVIP python>=3.8```
3. ```conda activate FunVIP```
4. ```pip install FunVIP```
5. run ```FunVIP --test Terrei --email [your email] ``` to check installation

* For upgrade use this command
``` pip install FunVIP --upgrade ```

### Linux
1. ```conda create -n FunVIP python>=3.8```
2. ```conda activate FunVIP```
3. ```pip install FunVIP```
4. ```conda config --add channels conda-forge```
5. ```conda install -c bioconda raxml iqtree "modeltest-ng>=0.1.7" mmseqs2 "blast>=2.12" mafft trimal gblocks fasttree```
6. run ```FunVIP --test Terrei --email [your email] ``` to check installation


* For intel mac system, this method probably work, but we couldn't test it because we don't have any intel mac device. We're looking for feedbacks in intel mac

### Apple Silicon Mac
1. ```CONDA_SUBDIR=osx-64 conda create -n FunVIP python>=3.8```
2. ```conda activate FunVIP```
3. ```conda config --env --set subdir osx-64```
4. ```conda install pyqt```
5. ```pip install FunVIP```
6. ```conda install -c bioconda raxml iqtree mmseqs2 "blast>=2.12" mafft trimal gblocks fasttree```
7. run ```FunVIP --test Terrei --email [your email] ``` to check installation

### Installation from source (For developers and core users)
* this is for developmental steps
1. ```git clone https://github.com/Changwanseo/FunVIP.git```
2. Move to ```~/FunVIP```
3. ```conda create -n FunVIP python=3.10```
4. ```conda activate FunVIP```
5. ```pip install ./```
6. run ```FunVIP --test Terrei --email [your email]``` to check installation


## Usage
```FunVIP --db {Your database file} --query {Your query file} --email {Your email} --gene {Your genes} --preset {fast or accurate}```

### Example
```FunVIP --db Penicillium.xlsx --query Query.xlsx --email {Your email} --gene ITS BenA RPB2 CaM --preset fast```


\* See documentation for detailed usage



<!--### GUI mode (\*Currently under development)
1. Go to ~/FunID-dev
2. ```streamlit run FunID_GUI.py```
* GUI run is on experimental
* If you want to edit GUI options, edit ```Option_manager.xlsx``` and variables in ```FunID_GUI.py```

### Server mode (\* Currently under development)-->



## How to make database?
![Fig 2 Database and command configuration of FunID (ver2) ](https://github.com/Changwanseo/FunVIP/assets/64393882/9ba71eb9-91e9-4c0b-ac60-b9b7be993694)




[See example database here](https://github.com/Changwanseo/FunVIP/blob/main/funvip/test_dataset/penicillium/DB/DB_Penicillium.xlsx)


<!--## 
## What query formats can be used?
#### Query formats can be either 
fasta (```.fa```, ```.fna```, ```.fas```, ```.fasta```, ```.txt```) or
tabular (```.xlsx```, ```.csv```,  ```.parquet```, ```.ftr```) form

- fasta form : Do not use ambiguous accessions in your fasta name. For example, accessions "A1234" and "A123" can be confused in pipeline. Section and genus name of the sequences will be automatically assigned according to your database. So if you want to fix it, use tabular form
- tabular form : your table should include ```ID```, and ```{gene names}``` (highly recommended for multigene analysis)-->

<!--## Tips for method selection
* SEARCH_METHOD : blast is faster for smaller dataset, while mmseqs are faster in huge dataset, but consumes a lot of memory
* ALIGNMENT_METHOD : currently mafft is only available.
* TRIMMING_METHOD : use trimal or gblocks, in your favor. gblocks usally cuts more, but can be differ by advanced option. Use none if you have enough time and resource for calculation
* MODEL_METHOD : model method is currently not working good enough please wait
* TREE_METHOD : fasttree is fastest, but least accurate (However, still a lot accurate than NJ tree). It is treated that iqtree is faster but slightly less accurate than raxml, but iqtree requires at least 1000 bootstrap. So in case of speed, raxml could be a little bit faster when low bootstrap selected-->

## Results
* ```Section Assignment.xlsx``` : Your clustering result is here. You can find which of your sequences are clustered to which section 
* ```Identification_result.xlsx``` : Your final identification result. Shows how your sequences were assigned to species level through tree-based identification
* ```report.xlsx``` : overall statistics about the tree. If your find taxon ends with numbers, these taxon are found to be paraphyletic, so should be checked
* ```/Tree/{section}_{gene}.svg``` : Final collapsed tree in svg format. Can be edited in vector graphics programs, or in powerpoint (by ungroup)
* ```/Tree/{section}_{gene}_original.svg ``` : Uncollapsed tree for inspection

## How does FunVIP work?
![figure1 - ver4](https://github.com/Changwanseo/FunID/assets/64393882/6a366d32-6aaf-4d0c-8102-8c7dd5fda4c2)




## License
GPL 3.0

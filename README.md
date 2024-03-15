[![DOI](https://zenodo.org/badge/588465720.svg)](https://zenodo.org/doi/10.5281/zenodo.10714946)

#### This is Beta release. Unstable. We're eagaring bug reports

## Scheduling
### Beta release part 1 (2023 Feburary ~ As paper published, ver 0.3)
- ~~Will operate normally most of the case~~ We tried our best, but it's still very buggy. Please report bugs for us
- Will be tested by our lab memebers to fix bugs and advance features

### Beta release part 2 (As paper published ~ When pipeline gets stabled, ver 0.4)
- Will be tested by peer taxonomists

### Stable release (ver 1.0)




# FunID
Fungal Identification Pipeline

A automatic tree-based sequence identification pipeline for fungal species

- Automatic tree-based identification
- Works with multigene
- Data validation algorithm implemented


# See [tutorial](https://github.com/Changwanseo/FunID/blob/main/tutorial/tutorial.md) for detailed usage


## Requirements
- Conda environment (See https://www.anaconda.com/products/individual to install)

<!--
## Installation with conda (May not work with Linux or Mac)
1. ```conda create -n FunID python=3.10```
2. ```conda activate FunID```
3. ```conda install -c cwseo FunID```
4. run ```FunID --test Penicillium ``` to check installation
If this one fails, use next one
-->

## Windows
1. Install visual c++ [here](https://visualstudio.microsoft.com/ko/visual-cpp-build-tools/)
2. ```conda create -n FunID python=3.9```
3. ```conda activate FunID```
4. ```pip install FunID```
5. run ```FunID --test Penicillium ``` to check installation

* For upgrade use this command
``` pip install FunID --upgrade ```

## Linux
1. ```conda create -n FunID python=3.9```
2. ```conda activate FunID```
3. ```pip install FunID```
4. ```conda install -c bioconda raxml iqtree modeltest-ng mmseqs2 "blast>=2.12" mafft trimal gblocks fasttree```
5. run ```FunID --test Penicillium ``` to check installation


* For intel mac system, this method probably work, but we couldn't test it because we don't have any intel mac device. We're looking for feedbacks in intel mac

## Apple Silicon Mac
1. ```CONDA_SUBDIR=osx-64 conda create -n FunID python=3.10```
2. ```conda activate FunID```
3. ```conda config --env --set subdir osx-64```
4. ```conda install pyqt```
5. ```pip install FunID```
6. ```conda install -c bioconda raxml iqtree mmseqs2 "blast>=2.12" mafft trimal gblocks fasttree```
7. run ```FunID --test Penicillium ``` to check installation

## Installation from source (For developers)
* this is for developmental steps
1. ```git clone https://github.com/Changwanseo/FunID.git```
2. Move to ```~/FunID```
3. ```conda create -n FunID python=3.10```
4. ```conda activate FunID```
5. ```pip install ./```
6. run ```FunID --test Penicillium ``` to check installation


## Usage
```FunID --db {Your database file} --query {Your query file} --email {Your email} --gene {Your genes} --preset {fast or accurate}```

### Example
```FunID --db Penicillium.xlsx --query Query.xlsx --email {Your email} --gene ITS BenA RPB2 CaM --preset fast```


\* See documentation for detailed usage



<!--### GUI mode (\*Currently under development)
1. Go to ~/FunID-dev
2. ```streamlit run FunID_GUI.py```
* GUI run is on experimental
* If you want to edit GUI options, edit ```Option_manager.xlsx``` and variables in ```FunID_GUI.py```

### Server mode (\* Currently under development)-->



## How to make database?
  ![Fig 2 FunID interface, usage, and input example](https://github.com/Changwanseo/FunID/assets/64393882/78e5dd58-a656-4395-bda9-0ee671411e2c)




[See example database here](https://github.com/Changwanseo/FunID/blob/main/funid/test_dataset/penicillium/DB/DB_Penicillium.xlsx)


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

## How does FunID work?
![figure1 - ver4](https://github.com/Changwanseo/FunID/assets/64393882/6a366d32-6aaf-4d0c-8102-8c7dd5fda4c2)




## License
GPL 3.0

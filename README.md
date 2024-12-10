
# FunVIP [![DOI](https://zenodo.org/badge/588465720.svg)](https://zenodo.org/doi/10.5281/zenodo.10714946)
### **Fun**gal **V**alidation & **I**dentification **P**ipeline
#### An automatic tree-based sequence identification and validation pipeline for fungal (or maybe other) species

- Automatic tree-based identification
- Works with multigene
- Data validation algorithm implemented

![figure1 - ver17A](https://github.com/user-attachments/assets/22a50a62-14e8-41a7-87a0-8f5a1f9c3f62)

This is Beta release. Bug reports are welcomed
<br><br/>

## Tutorial
* [Part 1 - Getting started!](https://github.com/Changwanseo/FunVIP/blob/main/tutorial/tutorial.md)
* [Part 2 - Preparing database and query](https://github.com/Changwanseo/FunVIP/blob/main/tutorial/tutorial2.md)
* [Advanced tips](https://github.com/Changwanseo/FunVIP/blob/main/tutorial/advanced.md)
<br><br/>
## Documentation
* See [Documentation](https://github.com/Changwanseo/FunVIP/blob/main/Documentation.md) for advanced usage !
<br><br/>
## Requirements
- Conda environment

\* See [https://www.anaconda.com/products/individual](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) for how to install conda environment
<br><br/>
## Installation
* [Windows](https://github.com/Changwanseo/FunVIP/blob/main/tutorial/installation.md##Windows)
* [Mac - apple silicon](https://github.com/Changwanseo/FunVIP/blob/main/tutorial/installation.md##Apple )
* [Linux](https://github.com/Changwanseo/FunVIP/blob/main/tutorial/installation.md##Linux)
* [from source](https://github.com/Changwanseo/FunVIP/blob/main/tutorial/installation.md##Installation)
<br><br/>
## Usage
```FunVIP --db {Your database file} --query {Your query file} --email {Your email} --gene {Your genes} --preset {fast or accurate}```
<br><br/>
### Example
```FunVIP --db Penicillium.xlsx --query Query.xlsx --email {Your email} --thread 8 --gene ITS BenA RPB2 CaM --preset fast```

\* See documentation for detailed usage
<br><br/>






## How to make database?
![figure1 - ver17B](https://github.com/user-attachments/assets/0043e6f6-7470-4c2b-bc96-b51f41c43ee4)




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

## Scheduling
1. Beta release part 1 (2023 Feburary ~ As paper published, ver 0.3)
Will be tested by our lab memebers to fix bugs and advance features
2. Beta release part 2 (As paper published ~ When pipeline gets stabled, ver 0.4)
Will be tested by peer taxonomists
3. Stable release (ver 1.0)

## License
[GPL 3.0](https://github.com/Changwanseo/FunVIP/blob/main/LICENSE)


<!--
## Installation with conda (May not work with Linux or Mac)
1. ```conda create -n FunVO{ python=3.10```
2. ```conda activate FunVIP```
3. ```conda install -c cwseo FunVIP```
4. run ```FunVIP --test Terrei --email [your email] ``` to check installation
If this one fails, use next one
-->
<!--### GUI mode (\*Currently under development)
1. Go to ~/FunID-dev
2. ```streamlit run FunID_GUI.py```
* GUI run is on experimental
* If you want to edit GUI options, edit ```Option_manager.xlsx``` and variables in ```FunID_GUI.py```
### Server mode (\* Currently under development)-->

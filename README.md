### This is developmental repository.

# FunID
Fungal Identification Pipeline

A automatic tree-based sequence identification pipeline for fungal species

- Precisely identify sequences from fungal strains / OTU
- Automatic database validation


![Workflow](https://user-images.githubusercontent.com/64393882/165916028-48f86e26-76a2-4e98-a066-d3372bb6ba61.png)


## Requirements
- Conda environment (See https://www.anaconda.com/products/individual to install)

## Installation
1. git clone https://github.com/Changwanseo/FunID-dev.git
2. Move to ~/FunID-dev
3. conda create -n FunID python=3.9
4. conda activate FunID
5. pip install ./

## Usage
```FunID -d {Your database file} -q {Your query file} -e {Your email}```

\* See documentation for detailed options 



<!--### GUI mode (\*Currently under development)
1. Go to ~/FunID-dev
2. ```streamlit run FunID_GUI.py```
* GUI run is on experimental
* If you want to edit GUI options, edit ```Option_manager.xlsx``` and variables in ```FunID_GUI.py```

### Server mode (\* Currently under development)-->



## How to make database?
Database should be tabular files, ```.xlsx, .csv, .parquet or .ftr``` files without unicodes (unicodes will be automatically edited, but not recomended)
Use parquet or feather datatype if your database is really big.
There are essential columns that should be included database
- ```ID``` : the numbers or symbols that were displayed in reports and figures. It can be NCBI accession, but not necessarily to be
- ```Genus``` : genus of the species
- ```Species``` : species epithet of the species. We recommend not to use 'sp.' only, because it can confused with multiple sp.s over clades. Please add numbers (like sp. 1) or other expressions (like aff. amilaria, tmpspecies1)
- ```{gene names}``` : each of the sequences should be added in {gene names} columns. Old database may condtain ```seq``` column instead of ```{gene names}```, which cannot be applied in multigene mode. {gene names} used in database should be recognized by "GENE" in ```Options.config```  


## What query formats can be used?
Query formats can be fasta (```.fa```, ```.fna```, ```.fas```, ```.fasta```, ```.txt```)
or tabular form (```.xlsx```, ```.csv```,  ```.parquet```, ```.ftr```) (use parquet or feather files for large db)

- fasta form : It is important to use not ambiguous accessions in your fasta name. For example, accessions "A1234" and "A123" can be confused in pipeline. Section and genus name of the sequences will be automatically assigned according to your database. So if you want to fix it, use tabular form
- tabular form : your table should include ```Accession```, and ```{gene names}```. ```Genus```, ```Species```, ```Section``` are optional


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

## License
Effective license will be added in the seperate file. This is an abstract.

0. Please wait for a while (may be by 2023 Feburary) for advanced (such as re-distribution) usage. We are working on finalizing stage
1. For softwares in /Bin/External_Programs, each of the software follows their own license 
2. In non-commercial use, free to use it and redistribute without edit
3. You may edit for non-commercial use, but should not redistribute without permission
4. Contact me with email for commercial use


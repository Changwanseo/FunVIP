
# FunID Tutorial

This is basic tutorial for how to run FunID

  

## Prerequisites

* [Conda installation](https://www.anaconda.com/products/individual)

* Database file

* Query file

  

### 1. Install FunID by following instructions

* [Windows](https://github.com/Changwanseo/FunID/#Windows)

* [Linux](https://github.com/Changwanseo/FunID/#Linux)

* [Mac - Apple Silicon](https://github.com/Changwanseo/FunID/#Apple-Silicon-Mac)

  

### 2. Prepare database and query file

Here, we'll going to use *Aspergillus* section *Terrei* dataset from the paper ([link will be added after publication]())

  

Download database and query file by clicking

* [Database](https://github.com/Changwanseo/FunID/blob/main/tutorial/FunID_Aspergillus_db.xlsx)

* [Query](https://github.com/Changwanseo/FunID/blob/main/tutorial/FunID_Aspergillus_query.xlsx)

  
  

### 3. Run FunID

Turn on your FunID conda environment

*  ```conda activate FunID```

  

Run FunID with single command.

*  ```FunID --db FunID_Aspergillus_db.xlsx --query FunID_Aspergillus_query.xlsx --email <your email> --gene ITS BenA CaM RPB2 --preset fast --level section```

  
  

### 4. See results

Open the result directory. The basic result directory will be named by timestamp when you run the program.
ex) 20240118-101519

* See ```<TIMESTAMP>.result.csv``` file to see overall results
* See ```./07_Tree/<TIMESTAMP>_Terrei_concatenated.svg``` to see phylogenetic tree




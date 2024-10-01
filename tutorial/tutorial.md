
# FunVIP Tutorial
* [Part 1 - Getting started](https://github.com/Changwanseo/FunVIP/edit/main/tutorial/tutorial.md) - this page
* [Part 2 - Preparing database and query]()
* [Part 3 - How to use FunVIP for final publication]()

This is basic tutorial for how to run FunVIP
In this tutorial, we will run FunVIP with already prepared input files.


## Prerequisites
### Conda environment should be installed to follow the tutorial 
* [Conda installation](https://www.anaconda.com/products/individual)


### 1. Install FunVIP by following instructions, depending on your os system

* [Windows](https://github.com/Changwanseo/FunVIP/#Windows)

* [Linux](https://github.com/Changwanseo/FunVIP/#Linux)

* [Mac - Apple Silicon](https://github.com/Changwanseo/FunVIP/#Apple-Silicon-Mac)
  

### 2. Prepare database and query file

Here, we'll going to use *Aspergillus* section *Terrei* dataset from the paper ([link will be added after publication]())

Download database and query file by clicking

* [Database](https://github.com/Changwanseo/FunVIP/tree/main/funvip/test_dataset/terrei/DB/FunVIP_Aspergillus_db.xlsx)

* [Query](https://github.com/Changwanseo/FunVIP/tree/main/funvip/test_dataset/terrei/Query/FunVIP_Aspergillus_query.xlsx)


If your directory looks like this you are going well

![image](https://github.com/user-attachments/assets/a6b65405-4828-4bf6-9589-1ffb3e7b4ae7)



### 3. Run FunVIP

Turn on your FunVIP conda environment

*  ```conda activate FunVIP```

![image](https://github.com/user-attachments/assets/4a28393c-3afe-45c0-a41e-38d7a2ed5380)

A. Your environment should be "FunVIP"

B. Your current directory should be the place where your database and query file exists
 
   

Now you are ready to go. Run FunVIP with single command.

*  ```FunVIP --db FunVIP_Aspergillus_db.xlsx --query FunID_Aspergillus_query.xlsx --email <your email> --gene ITS BenA CaM RPB2 --preset fast --level section```

  

### 4. See results

Open the result directory. The basic result directory will be named by timestamp when you run the program.
ex) 20240118-101519

* See ```<TIMESTAMP>.result.csv``` file to see overall results
* See ```./07_Tree/<TIMESTAMP>_Terrei_concatenated.svg``` to see phylogenetic tree




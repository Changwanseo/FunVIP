
# FunVIP Tutorial
* [Part 1 - Getting started](https://github.com/Changwanseo/FunVIP/edit/main/tutorial/tutorial.md)
* [Part 2 - Preparing database and query]() - this page
* [Part 3 - How to use FunVIP for final publication]()

This is basic tutorial for how to run FunVIP
In this tutorial, we will prepare own Sanghuangporus database and identify unknown Sanghuangporus sequneces from GenBank

### 1. Generate Sanghuangporus database from [Zhou et al. 2021](https://link.springer.com/article/10.1186/s43008-021-00059-x)
Open the following link about [work of Zhou et al. 2021](https://link.springer.com/article/10.1186/s43008-021-00059-x)

Copy the table from [Table 1](https://link.springer.com/article/10.1186/s43008-021-00059-x/tables/1) like this

![image](https://github.com/user-attachments/assets/513d58df-1acd-4a3b-a81a-62ae76f98522)

...and paste to your table managing program, such as microsoft excel

![image](https://github.com/user-attachments/assets/73c0fbb8-58b9-4687-a280-96661b80b8a2)

Before editing, save this file to Sanghuangporus_db.xlsx

Now, we sould edit the file in the format of FunVIP database, following this format

![image](https://github.com/user-attachments/assets/4d4b0d5a-af27-4fda-afa6-6e140d052eb4)

As long as you can fit it into the given format, you can proceed in any way. 

First, we will make "ID" column. We will change the column "Voucher No." to "ID" so that FunVIP can understand

![image](https://github.com/user-attachments/assets/7a4300fb-5858-4234-963d-635690bb671b)

Next, we will make "Genus" and "Species" column. But before that, we have to unmerge cells

You can use split by text to easily seperate genus name and specie name

Now we are almost there. Change the column name "GenBank No." to "ITS"

Also, we can see unnecessary subscripts in GenBank accession. Remove them.



### 2. Generate Sanghuangporus query from GenBank













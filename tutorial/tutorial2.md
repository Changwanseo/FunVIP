
# FunVIP Tutorial
* [Part 1 - Getting started](https://github.com/Changwanseo/FunVIP/blob/main/tutorial/tutorial.md)
* [Part 2 - Preparing database and query](https://github.com/Changwanseo/FunVIP/blob/main/tutorial/tutorial2.md) - this page
* [Part 3 - How to use FunVIP for final publication]()

This is basic tutorial for how to run FunVIP
In this tutorial, we will prepare own Sanghuangporus database and identify unknown _Sanghuangporus_ sequneces from GenBank

### 1. Generate _Sanghuangporus_ ITS database from [Zhou et al. 2021](https://link.springer.com/article/10.1186/s43008-021-00059-x)
Open the following link about [work of Zhou et al. 2021](https://link.springer.com/article/10.1186/s43008-021-00059-x)

Copy the table from [Table 1](https://link.springer.com/article/10.1186/s43008-021-00059-x/tables/1) like this

![image](https://github.com/user-attachments/assets/513d58df-1acd-4a3b-a81a-62ae76f98522)

...and paste to your spreadsheet managing software, such as Microsoft Excel

![image](https://github.com/user-attachments/assets/73c0fbb8-58b9-4687-a280-96661b80b8a2)

Before editing, save this file to Sanghuangporus_db.xlsx

Now, we sould edit the file in the format of FunVIP database, following this format

![image](https://github.com/user-attachments/assets/4d4b0d5a-af27-4fda-afa6-6e140d052eb4)

As long as you can fit it into the given format, you can proceed in any way. 



First, we will make "Genus" and "Species" column. But before that, we have to unmerge cells

![image](https://github.com/user-attachments/assets/ca684cfd-31dd-4c5e-8816-7e1a87e116ef)
![image](https://github.com/user-attachments/assets/a070fdf7-cbb2-49a6-9705-a0ed7000fd27)

You can use split by text to easily seperate genus name and specie name
![image](https://github.com/user-attachments/assets/8572cbc1-5dca-4fff-a642-2b78254356b1)
![image](https://github.com/user-attachments/assets/1a3d44a9-3697-4722-a38e-348568f0016f)
![image](https://github.com/user-attachments/assets/f42c0532-c25e-4b84-89e9-73b3b7641f65)
If needed, do some manual curations
![image](https://github.com/user-attachments/assets/3db67448-10a9-4222-bf28-ed247e20b2c4)
![image](https://github.com/user-attachments/assets/f1436f66-b89e-4bd1-902c-7812ccb65ba7)
Don't forget to change column names to "Genus" and "Species"
![image](https://github.com/user-attachments/assets/da8f8a22-db12-4db1-9bc1-6fb4a7fd7bdc)

Now we are almost there. As FunVIP is tool for multi-gene phylogeny, we have to designate which gene does the sequence means.
Change the column name "GenBank No." to "ITS"
![image](https://github.com/user-attachments/assets/ef2cad61-d710-425a-ab34-4891b89029cd)

Also, there are unnecessary subscripts in GenBank accession. Remove them.

![image](https://github.com/user-attachments/assets/8bf1c5ba-2a20-4656-8d07-8867bba0c0a6)

Last, we will make "ID" column. We will change the column "Voucher No." to "ID" so that FunVIP can understand
![image](https://github.com/user-attachments/assets/7a4300fb-5858-4234-963d-635690bb671b)

(Optional) you can remove unnecessary columns for clear look. It doesn't matter if you leave them or not

If you got table like this, now you have a good database for FunVIP!

![image](https://github.com/user-attachments/assets/9dc7d583-4ba2-431f-b4d9-b0c23d2da00c)




### 2. Generate _Sanghuangporus_ query from GenBank

In the case of genus _Sanghaungporus_, the database only includes single genetic marker, ITS

So in this case, FunVIP can accept FASTA file as query input.


![image](https://github.com/user-attachments/assets/9f88fb26-0c2e-4417-ad18-0517fde70c72)












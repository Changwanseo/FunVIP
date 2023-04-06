import os
import sys
import shutil
import re
import subprocess
import json
import pandas as pd
import numpy as np
import logging

from copy import deepcopy

from unidecode import unidecode
from Bio import SeqIO
from funid.src.tool import (
    initialize_path,
    get_genus_species,
    get_id,
    manage_unicode,
    mkdir,
)
from pathlib import Path
from funid.src import save
from funid.src.logics import isnewicklegal, isuniquecolumn, isvalidcolor
from funid.src.hasher import decode, newick_legal, hash_funinfo_list


# newick
NEWICK_ILLEGAL = (
    "(",
    '"',
    "[",
    ":",
    ";",
    "/",
    "[",
    "]",
    "{",
    "}",
    "(",
    ")",
    ",",
    "]",
    "+",
    '"',
    ")",
    " ",
)


# default funinfo class
class Funinfo:
    def __init__(self):
        self.original_id = ""  # original id, can be newick illegal
        self.id = ""  #  newick illegal characters removed
        self.hash = ""  # hash : HSXXHE
        self.description = ""  # full description from fasta
        self.ori_genus = ""  # original genus
        self.genus = ""  # final genus
        self.ori_species = ""  # original species
        self.bygene_species = {"concatenated": ""}  # species name designated by gene
        self.final_species = ""  # final species designated by concatenated analysis
        # species identifier if multiple branches with same species exists - ambiguous in result
        self.species_identifier = 0
        self.source = ""
        self.datatype = ""  # DB or Query
        self.group = ""  # original taxonomic group
        self.adjusted_group = ""  # adjusted taxonomic group by group clustering
        self.seq = {}
        self.unclassified_seq = []
        self.color = None  # color for highlighting in phylogenetic tree
        self.flat = []  # list of flat genes

    def update_seqrecord(self, seq, gene=None):
        self.description = seq.description
        self.genus, self.ori_species = get_genus_species(seq.description)

        if gene in self.seq:
            logging.error(
                f"More than 1 sequence for {gene} found for {self.id} during update_seqrecord"
            )
            raise Exception
        elif gene is None:
            self.unclassified_seq.append(str(seq.seq.ungap("-")))
        else:
            self.seq[gene] = str(seq.seq.ungap("-"))

        self.bygene_species[gene] = self.ori_species

    def update_seq(self, gene, seq):  # get input as Entrez seqrecord! Important!
        if gene in self.seq:
            if (
                self.seq[gene] != seq
            ):  # if more than 1 sequence per gene gets in, and if they are different
                logging.error(
                    f"More than 1 sequence for {gene} found for {self.id} during update_seq"
                )
                raise Exception
            else:
                pass
        else:
            self.seq[gene] = seq

        self.bygene_species[gene] = self.ori_species

        # Update concatenated
        self.bygene_species["concatenated"] = self.ori_species

    def update_description(self, description):
        self.description = description

    def update_genus(self, genus):
        # Try to solve illegal unicode characters
        if pd.isnull(genus):
            genus = ""

        # Genus with space causes error while mafft
        genus = genus.strip().replace(" ", "_")
        genus = manage_unicode(genus, column="Genus")

        # Check ambiguity
        if self.genus != "" and self.genus != genus:
            logging.error(
                f"Colliding genus info found for {funinfo}, {self.genus} and {genus}"
            )
            raise Exception

        # Update original if should
        if self.ori_genus == "":
            self.ori_genus = genus

        # Update genus
        self.genus = genus

    def update_ori_species(self, species):
        # Try to solve illegal unicode characters
        if pd.isnull(species):
            species = ""
        species = species.strip()
        species = manage_unicode(species, column="Species")

        # Check ambiguity
        if self.ori_species != "" and self.ori_species != species:
            logging.error(
                f"Colliding species info found for {self.original_id}, {self.ori_species} and {species}"
            )
            raise Exception

        # Update original if should
        if self.ori_species == "":
            self.ori_species = species

    def update_species(self, gene, species):
        self.bygene_species[gene] = species

    def update_group(self, group):
        # Try to solve illegal unicode characters
        if pd.isnull(group):
            group = ""

        # Group with space causes error while mafft
        group = group.strip().replace(" ", "_")
        group = manage_unicode(group, column="Group")

        # Check ambiguity
        if self.group != "" and self.group != group:
            logging.error(
                f"Colliding group info found for {funinfo}, {self.group} and {group}"
            )
            raise Exception

        # Update group
        self.group = group

    def update_color(self, color):
        if pd.isnull(color):
            color = None
            self.color = color
        else:
            color = manage_unicode(str(color).strip(), column="Color")
            if isvalidcolor(color) is True:
                self.color = color
            else:
                logging.error(
                    f"Color {color} does not seems to be valid svg color nor hex code"
                )
                raise Exception

        logging.debug(f"Updated color {color}")

    def update_datatype(self, datatype):
        # Available datatypes : db, query
        if not (datatype in ("db", "query", "outgroup")):
            logging.error(f"{datatype} is not available datatype")
            raise Exception

        # Check ambiguity
        if self.datatype != "" and self.datatype != datatype:
            logging.error(
                f"DEVELOPMENTAL ERROR : Colliding datatype found for {funinfo}, {self.datatype} and {datatype}"
            )
            raise Exception

        self.datatype = datatype

    def update_id(self, id_, regexs=None):
        if not regexs == None:
            id_ = get_id(id_, tuple(regexs))

        # if cannot find id by regex
        id_ = str(id_)
        self.original_id = id_
        id_ = newick_legal(id_)
        self.id = id_

    def update_hash(self, n):
        self.hash = f"HS{n}HE"

    def __repr__(self):
        return f"FI: {self.id}"

    def __hash__(self):
        return hash((self.original_id, self.hash, self.description))

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        return (
            self.original_id == other.original_id
            and self.hash == other.hash
            and self.description == other.description
        )


# getting data input from fasta file
def input_fasta(path, opt, fasta_list, funinfo_dict, datatype):
    # initialize path to use function "get_genus_species"
    initialize_path(path)

    # Fasta files only
    for file in fasta_list:
        # Copy input files to designation
        if datatype == "query":
            shutil.copy(file, path.out_query)
        elif datatype == "db":
            shutil.copy(file, path.out_db)

        tmp_list = []

        full = ".".join(file.split(".")[:-1])  # full file path
        name = ".".join(file.split("/")[-1].split(".")[:-1])  # name

        logging.info(f"{file}: Fasta file")

        try:
            fasta_list = list(SeqIO.parse(file, "fasta"))
            for seq in fasta_list:
                if not opt.regex == None:
                    id_ = get_id(seq.description, tuple(opt.regex))
                else:
                    id_ = seq.description

                id_ = newick_legal(id_)

                if id_ in funinfo_dict:
                    funinfo_dict[id_].update_seqrecord(seq)
                    funinfo_dict[id_].update_datatype(datatype)
                    funinfo_dict[id_].update_group("")
                    if get_genus_species(seq.description)[0] != "":
                        funinfo_dict[id_].update_genus(
                            get_genus_species(seq.description)[0]
                        )

                    if get_genus_species(seq.description)[1] != "":
                        funinfo_dict[id_].update_ori_species(
                            get_genus_species(seq.description)[1]
                        )

                else:
                    newinfo = Funinfo()
                    newinfo.update_seqrecord(seq)
                    newinfo.update_datatype(datatype)
                    newinfo.update_group("")  # because group not designated yet
                    # id by regex match
                    newinfo.update_id(seq.description, regexs=opt.regex)
                    if get_genus_species(seq.description)[0] != "":
                        newinfo.update_genus(get_genus_species(seq.description)[0])

                    if get_genus_species(seq.description)[1] != "":
                        newinfo.update_ori_species(
                            get_genus_species(seq.description)[1]
                        )

                    funinfo_dict[id_] = deepcopy(newinfo)

        except:
            logging.warning(f"{file} does not seems to be valid fasta file skipping")

    return funinfo_dict


# getting datafile from excel or tabular file
def input_table(funinfo_dict, path, opt, table_list, datatype):
    string_error = 0

    initialize_path(path)  # this one is ugly
    df_list = []

    # extensionto filetype translation
    dict_extension = {
        ".csv": "csv",
        ".tsv": "csv",
        ".xlsx": "excel",
        ".xls": "excel",
        ".parquet": "parquet",
        ".ftr": "feather",
        ".feather": "feather",
    }

    # Running table by table operations
    for table in table_list:
        # Read each of the table by each of the extensions
        flag_read_table = 0
        for extension in dict_extension:
            if table.endswith(extension):
                try:
                    if dict_extension[extension] == "csv":
                        df = pd.read_csv(table)
                        flag_read_table = 1
                    elif dict_extension[extension] == "excel":
                        df = pd.read_excel(table)
                        flag_read_table = 1
                    elif dict_extension[extension] == "parquet":
                        df = pd.read_parquet(table, engine="pyarrow")
                        flag_read_table = 1
                    elif dict_extension[extension] == "feather":
                        df = pd.read_feather(table, use_threads=True)
                        flag_read_table = 1
                except:
                    logging.error(
                        f"Table {table} cannot be read as {dict_extension[extension]} file. Please check files, extensions and seperators"
                    )
                    raise Exception

        if flag_read_table == 0:
            logging.error(
                f"Table {table} cannot recognized as either csv, xlsx, feather or parquet. Please check if extensions endswith .csv, .tsv, .xlsx, .parquet, .ftr or .feather"
            )
            raise Exception

        # Lower case column names
        df.columns = df.columns.str.lower()
        df.columns = df.columns.str.strip()
        df_list.append(df)

        # Clean up columns
        # Check if "id" column exists and unique
        flag_id = isuniquecolumn(
            list_column=list(df.columns), column=("accession", "id"), table_name=table
        )

        # If old accession column, change to id
        if flag_id == "accession":
            df.rename(columns={"accession": "id"}, inplace=True)

        # Check if "genus" column exists and unique
        # Column "genus" is mandatory in db, and optional in query
        check_none = True if datatype == "db" else False
        flag_genus = isuniquecolumn(
            list_column=list(df.columns),
            column=tuple(("genus",)),
            table_name=table,
            check_none=check_none,
        )

        # Check if "species" column exists and unique
        # Column "species" is mandatory in db, and optional in query
        check_none = True if datatype == "db" else False
        flag_species = isuniquecolumn(
            list_column=list(df.columns),
            column=tuple(("species",)),
            table_name=table,
            check_none=check_none,
        )

        # Check "opt.level" column exists and unique
        # Column "opt.level" is mandatory in db, and optional inquery
        check_none = True if datatype == "db" else False
        flag_level = isuniquecolumn(
            list_column=list(df.columns),
            column=tuple((opt.level,)),
            table_name=table,
            check_none=check_none,
        )

        # Check column color
        # color column is optional
        flag_color = isuniquecolumn(
            list_column=list(df.columns),
            column=tuple(("color",)),
            table_name=table,
            check_none=False,
        )

        # Sequence column operations, download sequences with GenMine
        if 1:  # If download on/off option added, change this part
            download_dict = {}  # for downloaded sequences
            download_set = set()
            # 1 letter + 5 digit regex should be last, because they overlap with 2 letter + 6 digit ids
            regex_genbank = r"(([A-Z]{1}[0-9]{5})(\.[0-9]{1}){0,1})|(([A-Z]{2}[\_]{0,1}[0-9]{6}){1}([\.][0-9]){0,1})"

            # if gene name were not designated by user, use seq
            opt.gene = list(set([gene.lower().strip() for gene in opt.gene]))

            # find all NCBI accessions in seq
            for gene in opt.gene:
                if isuniquecolumn(
                    list_column=df.columns,
                    column=tuple((gene,)),
                    table_name=table,
                    check_none=False,
                ):
                    for n, _ in enumerate(df[gene]):
                        if not (pd.isna(df[gene][n])):
                            if re.search(regex_genbank, df[gene][n]):
                                # remove unexpected indents with strip
                                download_set.add(df[gene][n].strip())

            # if NCBI accessions detected in sequence part, download it
            if len(download_set) > 0:
                logging.info(
                    f"Running GenMine to download {len(download_set)} sequences from GenBank"
                )
                logging.info(download_set)

                # Write GenMine input file
                with open(f"{path.GenMine}/Accessions.txt", "w") as fg:
                    for acc in download_set:
                        fg.write(f"{acc.strip()}\n")

                # Run GenMine
                cmd = f"GenMine -c {path.GenMine}/Accessions.txt -o {path.GenMine} -e {opt.email}"
                return_code = subprocess.call(cmd, shell=True)

                if return_code != 0:
                    logging.error(f"GenMine failed with return_code: {return_code}")
                    logging.error(f"This is usually ncbi server connection error")
                    logging.error(
                        f"Check your network problems, or replace all accession numbers in your db file with actual sequences to run locally"
                    )
                    raise Exception

                GenMine_df_list = [
                    file
                    for file in os.listdir(path.GenMine)
                    if file.endswith("_transformed.xlsx")
                ]

                if len(GenMine_df_list) == 1:
                    download_df = pd.read_excel(f"{path.GenMine}/{GenMine_df_list[0]}")

                    # Generate download_dict (I think this can be done with pandas operation, but a bit tricky. Will be done later)
                    for n, acc in enumerate(download_df["acc"]):
                        download_dict[acc.strip()] = download_df["seq"][n]

                    # replace accession to sequence downloaded
                    for n, _ in enumerate(df["id"]):
                        for gene in opt.gene:
                            if gene in df.columns:
                                if not (pd.isna(df[gene][n])):
                                    if df[gene][n].strip() in download_dict:
                                        df[gene][n] = download_dict[
                                            df[gene][n].strip().split(".")[0]
                                        ]
                                    # Removed download fails
                                    elif (
                                        not (df[gene][n].strip() in download_dict)
                                        and df[gene][n].strip() in download_set
                                    ):
                                        df[gene][n] = ""

                    # Remove GenMine results to prevent collision with next set
                    for file in os.listdir(path.GenMine):
                        if "transformed.xlsx" in file:
                            os.remove(f"{path.GenMine}/{file}")

                elif len(GenMine_df_list) == 0:
                    logging.warning(
                        f"None of the GenMine results were succesfully parsed"
                    )
                else:
                    logging.error(
                        f"DEVELOPMENTAL ERROR: Multiple GenMine result colliding!"
                    )
                    raise Exception

        # Generate funinfo by each row
        for n, acc in enumerate(df["id"]):
            # Check if each of the ids are unique
            # Remove non-unicode first
            new_acc = True
            df["id"][n] = manage_unicode(str(df["id"][n]), column="ID/Accession", row=n)
            # Generate funinfo for each id
            if df["id"][n] in funinfo_dict:
                newinfo = funinfo_dict[df["id"][n]]
                new_acc = False
                logging.warning(f"Duplicate id {df['id'][n]} found!")
            else:
                funinfo_dict[df["id"][n]] = Funinfo()
                newinfo = funinfo_dict[df["id"][n]]
                newinfo.update_id(df["id"][n])

            # if flag_genus is true, try to parse genus
            if not (flag_genus is None or flag_genus is False):
                newinfo.update_genus(df["genus"][n])

            # if flag_species is true, try to parse species
            if not (flag_species is None or flag_species is False):
                newinfo.update_ori_species(df["species"][n])

            # if flag_level is true, try to parse the optimal taxonomic group
            if not (flag_level is None or flag_level is False):
                newinfo.update_group(df[flag_level][n])

            # if flag_color is true, try to parse color for taxon
            if not (flag_color is None or flag_color is False):
                newinfo.update_color(df[flag_color][n])

            # update datatype
            newinfo.update_datatype(datatype)

            # parse each of the genes
            for gene in opt.gene:
                seq_error = 0
                if gene in df.columns:
                    if not (pd.isna(df[gene][n])) or str(df[gene][n]).strip() == "":
                        # skip blank sequences
                        if df[gene][n].startswith(">"):
                            # remove fasta header
                            seq_string = "".join(df[gene][n].split("\n")[1:])
                        else:
                            seq_string = df[gene][n]

                        # adjust seq_string
                        seq_string = seq_string.replace("\n", "").replace(" ", "")

                        # Finding if sequence contains error
                        seq_error_cnt = 0
                        seq_error_list = []
                        for x in seq_string:  # x is every character of sequence
                            if not x.lower() in "acgtryswkmbdhvn-.":
                                seq_error_cnt += 1
                                seq_error_list.append(x)

                        if seq_error_cnt > 0:
                            logging.warning(
                                f"Illegal DNA character {seq_error_list} found in {gene} of DB {df['id'][n]}"
                            )
                        elif seq_error_cnt == 0:
                            # remove gaps for preventing BLAST error
                            newinfo.update_seq(
                                gene, seq_string.replace("-", "").replace(".", "")
                            )

        # After successfully parsed this table, save if
        save.save_df(
            df,
            f"{path.out_db}/Saved_{'.'.join(table.split('.')[:-1])}.{opt.matrixformat}",
            fmt=opt.matrixformat,
        )

    # make it to list at last
    # list_funinfo = [funinfo_dict[x] for x in funinfo_dict]

    return funinfo_dict


def db_input(funinfo_dict, opt, path) -> list:
    # Get DB input
    logging.info(f"Input DB list: {opt.db}")
    # db_namelist = [str(os.path.basename(db)) for db in opt.db]

    funinfo_dict = input_table(
        funinfo_dict=funinfo_dict, path=path, opt=opt, table_list=opt.db, datatype="db"
    )

    """
    #### MOVE THIS TO TABLE INPUT
    # do it after save_db option enabled
    for n, db in enumerate(db_namelist):
        save.save_df(
            df_list[n],
            f"{path.out_db}/Saved_{'.'.join(db.split('.')[:-1])}.xlsx",
            fmt=opt.matrixformat,
        )
    """

    # validate dataset
    # if only one group exists, outgroup cannot work
    group_set = set(funinfo_dict[key].group for key in funinfo_dict)
    group_set.discard("")
    if len(group_set) <= 1:
        logging.error(
            f"Only 1 group soley detected : {group_set}. Please add outgroup sequences"
        )
        raise Exception

    # check if minimum outgroup number count exceeds minimum group
    group_cnt_dict = {x: 0 for x in group_set}
    for key in funinfo_dict:
        if type(funinfo_dict[key].group) is str:
            if funinfo_dict[key].group in group_set:
                group_cnt_dict[funinfo_dict[key].group] += 1

    if any(group_cnt_dict[x] < opt.maxoutgroup for x in group_cnt_dict):
        logging.warning(
            f"Sequences in database of some group has lower number than MINIMUM_OUTGROUP_COUNT. It may cause error when outgroup selection, or may select not most appropriate outgroup to group. Please lower number of MINIMUM_OUTGROUP_COUNT in option or add more sequences to these groups"
        )

    return funinfo_dict


def query_input(funinfo_dict, opt, path):
    query_fasta = [
        file
        for file in opt.query
        if any(file.endswith(x) for x in (".fa", ".fna", ".fas", ".fasta", ".txt"))
    ]
    query_table = [
        file
        for file in opt.query
        if any(
            file.endswith(x)
            for x in (".csv", ".tsv", ".xlsx", ".ftr", ".feather", ".parquet")
        )
    ]

    funinfo_dict = input_table(
        funinfo_dict=funinfo_dict,
        path=path,
        opt=opt,
        table_list=query_table,
        datatype="query",
    )
    funinfo_dict = input_fasta(
        path=path,
        opt=opt,
        fasta_list=query_fasta,
        funinfo_dict=funinfo_dict,
        datatype="query",
    )

    # Save query initially in raw format, because gene has not been assigned yet
    for file in query_table:
        shutil.copy(f"{file}", f"{path.out_query}")

    for file in query_fasta:
        shutil.copy(f"{file}", f"{path.out_query}")

    logging.info(
        f"Total {len([funinfo_dict[key].datatype =='query' for key in funinfo_dict.keys()])} sequences parsed from query"
    )

    return funinfo_dict


# combined db and query input
def data_input(V, R, opt, path):
    funinfo_dict = {}

    # get database input
    funinfo_dict = db_input(funinfo_dict=funinfo_dict, opt=opt, path=path)

    # get query input
    funinfo_dict = query_input(funinfo_dict=funinfo_dict, opt=opt, path=path)

    # combine all data
    V.list_FI = [funinfo_dict[key] for key in funinfo_dict]

    # hashing data for safety in tree analysis
    V.list_FI = hash_funinfo_list(V.list_FI)

    # make hash dict
    for FI in V.list_FI:
        V.dict_hash_FI[FI.hash] = FI

    return V, R, opt


# This function updates initial input query files to saved files with matrices and downloaded sequences
def update_queryfile(V, path, opt):
    # Dictionary to generate dataframe
    out_dict = {"ID": []}

    for gene in opt.gene:
        out_dict[gene] = []

    for FI in V.list_FI:
        if FI.datatype == "query":
            out_dict["ID"].append(FI.original_id)
            for gene in opt.gene:
                if gene in FI.seq:
                    out_dict[gene].append(FI.seq[gene])
                else:
                    out_dict[gene].append("")

    df_out = pd.DataFrame(out_dict)
    save.save_df(
        df_out, f"{path.out_query}/Saved_query.{opt.matrixformat}", fmt=opt.matrixformat
    )

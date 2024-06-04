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
from funvip.src.tool import (
    initialize_path,
    get_genus_species,
    get_id,
    manage_unicode,
    mkdir,
)
from pathlib import Path
from time import sleep

from funvip.src import save
from funvip.src.logics import isnewicklegal, isuniquecolumn, isvalidcolor
from funvip.src.hasher import decode, newick_legal, hash_funinfo_list


## newick illegal characters
# fmt: off
NEWICK_ILLEGAL = ("(",'"',"[",":",";","/","[","]","{","}","(",")",",","]","+",'"',")"," ",)
# fmt: on


# Default funinfo class
# Abbreviated as FI in most of the codes
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
        self.flat = []  # list of flat species in concatenated tree
        self.issues = set()  # list of issues to this FI
        self.source = set()  # source which FI was from

    def update_source(self, source):
        self.source.add(source)

    def update_seqrecord(self, seq, gene=None):
        error = None
        self.description = seq.description
        self.genus, self.ori_species = get_genus_species(seq.description)

        if gene in self.seq:
            error = f"In {self.source}, More than 1 sequence for {gene} found for {self.id} during update_seqrecord"
        elif gene is None:
            self.unclassified_seq.append(str(seq.seq.ungap("-")))
        else:
            self.seq[gene] = str(seq.seq.ungap("-"))

        self.bygene_species[gene] = self.ori_species

        return error

    def update_seq(self, gene, seq):  # get input as Entrez seqrecord! Important!
        error = None
        if gene in self.seq:
            if (
                self.seq[gene] != seq
            ):  # if more than 1 sequence per gene gets in, and if they are different
                error = f"In {self.source}, More than 1 sequence for {gene} found for {self.id} during update_seq"
            else:
                pass
        else:
            self.seq[gene] = seq

        self.bygene_species[gene] = self.ori_species

        # Update concatenated
        self.bygene_species["concatenated"] = self.ori_species

        return error

    def update_description(self, description):
        self.description = description

    def update_genus(self, genus):
        error = None
        # Try to solve illegal unicode characters
        if pd.isnull(genus):
            genus = ""

        # Genus with space causes error while mafft
        # Genus with slash causes error while tree construction
        genus = genus.strip().replace(" ", "_").replace("/", "_")
        genus = manage_unicode(genus, column="Genus")

        # Check ambiguity
        if self.genus != "" and self.genus != genus:
            error = f"In {self.source}, Colliding genus info found for {self.original_id}, {self.genus} and {genus}"

        # Update original if should
        if self.ori_genus == "":
            self.ori_genus = genus

        # Update genus
        self.genus = genus

        return error

    def update_ori_species(self, species):
        error = None
        # Try to solve illegal unicode characters
        if pd.isnull(species):
            species = ""

        # Genus with space causes error while mafft
        # Genus with slash causes error while tree construction
        species = species.strip().replace(" ", "_").replace("/", "_")
        species = manage_unicode(species, column="Species")

        # Check ambiguity
        if self.ori_species != "" and self.ori_species != species:
            error = f"In {self.source}, Colliding species info found for {self.original_id}, {self.ori_species} and {species}"

        # Update original if should
        if self.ori_species == "":
            self.ori_species = species

        return error

    def update_species(self, gene, species):
        self.bygene_species[gene] = species

    def update_group(self, group):
        error = None
        # Try to solve illegal unicode characters
        if pd.isnull(group):
            group = ""

        # Group with space causes error while mafft
        group = group.strip().replace(" ", "_")
        group = manage_unicode(group, column="Group")

        # Check ambiguity
        if self.group != "" and self.group != group:
            error = f"In {self.source}, Colliding group info found for {self.original_id}, {self.group} and {group}"

        # Update group
        self.group = group

        return error

    def update_color(self, color):
        if pd.isnull(color):
            color = None
            self.color = color
        else:
            color = manage_unicode(str(color).strip(), column="Color")
            if isvalidcolor(color) is True:
                self.color = color
            else:
                logging.warning(
                    f"Color {color} does not seems to be valid svg color nor hex code. Using default color"
                )
                self.color = None

    def update_datatype(self, datatype):
        flag = 0
        error = None
        # Available datatypes : db, query
        if not (datatype in ("db", "query", "outgroup")):
            logging.error(f"DEVELOPMENTAL ERROR: {datatype} is not available datatype")
            raise Exception

        # Check ambiguity
        if self.datatype != "" and self.datatype != datatype:
            error = f"In {self.source}, Colliding datatype found for {self.original_id}, {self.datatype} and {datatype}"

        self.datatype = datatype

        return error

    def update_id(self, id_, regexs=None):
        if not regexs == None:
            id_ = get_id(id_, tuple(regexs))

        # if cannot find id by regex
        # Even for original, new line character makes significant errors while working fasta, so remove it
        id_ = str(id_).replace("\n", " ")
        self.original_id = id_
        id_ = newick_legal(id_)
        self.id = id_

    def update_hash(self, n):
        self.hash = f"HS{n}HE"

    def get_issue_str(self):
        flat_issues = []
        alignfail_issues = []
        polyphyly_issues = []
        other_issues = []
        for issue in self.issues:
            if issue.startswith("flat"):
                flat_issues.append(issue)
            elif issue.startswith("alignfail"):
                alignfail_issues.append(issue)
            elif issue.startswith("polyphyly"):
                polyphyly_issues.append(issue)
            else:
                other_issues.append(issue)

        issues_string = ", ".join(other_issues)
        if len(flat_issues) > 0:
            flat_issue_str = "flat:" + "/".join(
                sorted([issue.split(":")[1] for issue in flat_issues])
            )
            if len(issues_string) == 0:
                issues_string = flat_issue_str
            else:
                issues_string += f", {flat_issue_str}"

        if len(alignfail_issues) > 0:
            alignfail_issue_str = "alignfail:" + "/".join(
                sorted([issue.split(":")[1] for issue in alignfail_issues])
            )
            if len(issues_string) == 0:
                issues_string = alignfail_issue_str
            else:
                issues_string += f", {alignfail_issue_str}"

        if len(polyphyly_issues) > 0:
            polyphyly_issue_str = "polyphyly:" + "/".join(
                sorted([issue.split(":")[1] for issue in polyphyly_issues])
            )
            if len(issues_string) == 0:
                issues_string = polyphyly_issue_str
            else:
                issues_string += f", {polyphyly_issue_str}"

        return issues_string

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
    warnings = []
    errors = []

    # initialize path to use function "get_genus_species"
    initialize_path(path)

    errors = []
    warnings = []

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

        seq_list = list(SeqIO.parse(file, "fasta"))
        for seq in seq_list:
            if not opt.regex == None:
                id_ = get_id(seq.description, tuple(opt.regex))
            else:
                id_ = seq.description

            id_ = newick_legal(id_)

            if id_ in funinfo_dict:
                funinfo_dict[id_].update_source(file)
                errors.append(funinfo_dict[id_].update_seqrecord(seq))
                errors.append(funinfo_dict[id_].update_datatype(datatype))
                errors.append(funinfo_dict[id_].update_group(""))
                if get_genus_species(seq.description)[0] != "":
                    errors.append(
                        funinfo_dict[id_].update_genus(
                            get_genus_species(seq.description)[0]
                        )
                    )

                if get_genus_species(seq.description)[1] != "":
                    errors.append(
                        funinfo_dict[id_].update_ori_species(
                            get_genus_species(seq.description)[1]
                        )
                    )

            # For new Funinfo
            else:
                newinfo = Funinfo()
                newinfo.update_source(file)
                errors.append(newinfo.update_seqrecord(seq))
                errors.append(newinfo.update_datatype(datatype))
                errors.append(
                    newinfo.update_group("")
                )  # because group not designated yet
                # id by regex match
                newinfo.update_id(seq.description, regexs=opt.regex)
                if get_genus_species(seq.description)[0] != "":
                    errors.append(
                        newinfo.update_genus(get_genus_species(seq.description)[0])
                    )

                if get_genus_species(seq.description)[1] != "":
                    errors.append(
                        newinfo.update_ori_species(
                            get_genus_species(seq.description)[1]
                        )
                    )

                funinfo_dict[id_] = deepcopy(newinfo)

        """
        try:
            seq_list = list(SeqIO.parse(file, "fasta"))
            for seq in seq_list:
                if not opt.regex == None:
                    id_ = get_id(seq.description, tuple(opt.regex))
                else:
                    id_ = seq.description

                id_ = newick_legal(id_)

                if id_ in funinfo_dict:
                    funinfo_dict[id_].update_source(file)
                    errors.append(funinfo_dict[id_].update_seqrecord(seq))
                    errors.append(funinfo_dict[id_].update_datatype(datatype))
                    errors.append(funinfo_dict[id_].update_group(""))
                    if get_genus_species(seq.description)[0] != "":
                        errors.append(
                            funinfo_dict[id_].update_genus(
                                get_genus_species(seq.description)[0]
                            )
                        )

                    if get_genus_species(seq.description)[1] != "":
                        errors.append(
                            funinfo_dict[id_].update_ori_species(
                                get_genus_species(seq.description)[1]
                            )
                        )

                # For new Funinfo
                else:
                    newinfo = Funinfo()
                    newinfo[id_].update_source(file)
                    errors.append(newinfo.update_seqrecord(seq))
                    errors.append(newinfo.update_datatype(datatype))
                    errors.append(
                        newinfo.update_group("")
                    )  # because group not designated yet
                    # id by regex match
                    newinfo.update_id(seq.description, regexs=opt.regex)
                    if get_genus_species(seq.description)[0] != "":
                        errors.append(
                            newinfo.update_genus(get_genus_species(seq.description)[0])
                        )

                    if get_genus_species(seq.description)[1] != "":
                        errors.append(
                            newinfo.update_ori_species(
                                get_genus_species(seq.description)[1]
                            )
                        )

                    funinfo_dict[id_] = deepcopy(newinfo)

        except:
            errors.append(f"{file} does not seems to be valid fasta file")
        """

        if len(seq_list) == 0:
            errors.append(f"Fasta file {file} seems to be empty please check")
            # raise Exception

        errors = [err for err in errors if not (err is None)]

    return funinfo_dict, warnings, errors


# getting datafile from excel or tabular file
def input_table(funinfo_dict, path, opt, table_list, datatype):
    # Whether to check if GenMine has run
    GenMine_flag = 0
    string_error = 0

    initialize_path(path)  # this one is ugly
    df_list = []

    warnings = []
    errors = []

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
        # This is part is double validation after options
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

                    # To prevent nan error
                    df = df.fillna("")
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
            # 1 letter + 5 digit regex should be last, because they overlap with 2 letter + 6 digit ids / shotgun
            regex_genbank = r"(([A-Z]{1}[0-9]{5})(\.[0-9]{1}){0,1})|(([A-Z]{2}[\_]{0,1}[0-9]{6}){1}([\.][0-9]){0,1})|(([A-Z]{4}[0-9]{8})(\.[0-9]{1}){0,1})|(([A-Z]{6}[0-9]{9,})(\.[0-9]{1}){0,1})"

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
                GenMine_flag = 1
                if opt.email == "":
                    logging.error(
                        "Your database includes GenBank accession but email not provided. --email is required to connect to GenBank"
                    )
                    raise Exception

                logging.info(
                    f"Running GenMine to download {len(download_set)} sequences from GenBank"
                )

                # Write GenMine input file
                with open(f"{path.GenMine}/Accessions.txt", "w") as fg:
                    for acc in download_set:
                        fg.write(f"{acc.strip()}\n")

                # Run GenMine
                accession_path = f"{path.GenMine}/Accessions.txt"
                # To prevent space errors in windows
                if " " in path.GenMine:
                    accession_path = f'"{accession_path}"'
                    GenMine_path = f'"{path.GenMine}"'
                else:
                    GenMine_path = path.GenMine

                cmd = f"GenMine -c {accession_path} -o {GenMine_path} -e {opt.email}"
                logging.info(cmd)

                sleep(5)  # To run GenMine safetly between run and run
                return_code = subprocess.call(cmd, shell=True)

                if return_code != 0:
                    logging.error(f"GenMine failed with return_code: {return_code}")
                    logging.error(f"This is usually NCBI server connection error")
                    logging.error(
                        f"Check your network problems, or replace all accession numbers in your db file with actual sequences to run locally"
                    )
                    raise Exception

                GenMine_outdir_candidate = [
                    directory
                    for directory in os.listdir(path.GenMine)
                    if not ("." in directory)
                ]

                GenMine_df_list = []
                for directory in GenMine_outdir_candidate:
                    for file in os.listdir(f"{path.GenMine}/{directory}/"):
                        if file.endswith("_transformed.xlsx"):
                            GenMine_df_list.append(f"{path.GenMine}/{directory}/{file}")

                if len(GenMine_df_list) == 1:
                    download_df = pd.read_excel(GenMine_df_list[0])

                    # Generate download_dict (I think this can be done with pandas operation, but a bit tricky. Will be done later)
                    for n, acc in enumerate(download_df["acc"]):
                        download_dict[acc.strip()] = download_df["seq"][n]

                    # replace accession to sequence downloaded
                    def update_from_GenMine(string):
                        string = str(string)
                        accession_wo_version = string.strip().split(".")[0]
                        accession_w_version = string.strip()

                        # Do not consider sequence version
                        if accession_wo_version in download_dict:
                            return download_dict[accession_wo_version]

                        # Remove accessions which failed download to prevent conflict with sequence validation
                        elif not (accession_wo_version in download_dict) and (
                            accession_w_version in download_set
                        ):
                            logging.warning(f"Failed updating {string} using GenMine")
                            return ""

                        # For sequence input
                        else:
                            return string

                    for gene in opt.gene:
                        if gene in df.columns:
                            df[gene] = df[gene].apply(update_from_GenMine)

                    # Remove GenMine results to prevent collision with next set
                    for directory in GenMine_outdir_candidate:
                        shutil.rmtree(f"{path.GenMine}/{directory}")

                elif len(GenMine_df_list) == 0:
                    warnings.append(
                        f"In table {table}, None of the GenMine results were succesfully parsed"
                    )

                else:
                    logging.error(
                        f"DEVELOPMENTAL ERROR: Multiple GenMine result colliding!"
                    )
                    raise Exception

        # Manage unicode for ID
        df["id"] = df["id"].apply(
            lambda x: manage_unicode(str(x), column="ID/Accession")
        )

        # To prevent errors on genus / spcies column
        # if this function operates with query mode, there might be no genus or species column
        if "genus" in df.columns:
            df["genus"] = df["genus"].apply(lambda x: x.replace(" ", "_"))
            df["genus"] = df["genus"].apply(lambda x: x.replace("/", "_"))
        if "species" in df.columns:
            df["species"] = df["species"].apply(lambda x: x.replace(" ", "_"))
            df["species"] = df["species"].apply(lambda x: x.replace("/", "_"))

        # Empty id check
        empty_error = []
        for n, acc in enumerate(df["id"]):
            if df["id"][n].strip() == "" or df["id"][n].strip() == "-":
                empty_error.append(n)

        if len(empty_error) > 0:
            errors.append(f"In table {table}, Empty id found, line {empty_error}!")

        # Generate funinfo by each row
        for n, acc in enumerate(df["id"]):
            # Check if each of the ids are unique
            # Remove non-unicode first
            new_acc = True

            # Generate funinfo for each id
            # Duplicate id check
            if df["id"][n] in funinfo_dict:
                newinfo = funinfo_dict[df["id"][n]]
                # Update source
                newinfo.update_source(table)
                new_acc = False
                warnings.append(
                    f"Among table {newinfo.source}, Duplicate id {df['id'][n]} found!"
                )
            else:
                funinfo_dict[df["id"][n]] = Funinfo()
                newinfo = funinfo_dict[df["id"][n]]
                newinfo.update_id(df["id"][n])
                # Update source
                newinfo.update_source(table)

            # if flag_genus is true, try to parse genus
            if not (flag_genus is None or flag_genus is False):
                errors.append(newinfo.update_genus(df["genus"][n]))

            # if flag_species is true, try to parse species
            if not (flag_species is None or flag_species is False):
                errors.append(newinfo.update_ori_species(df["species"][n]))

            # if flag_level is true, try to parse the optimal taxonomic group
            if not (flag_level is None or flag_level is False):
                errors.append(newinfo.update_group(df[flag_level][n]))

            # if flag_color is true, try to parse color for taxon
            if not (flag_color is None or flag_color is False):
                newinfo.update_color(df[flag_color][n])

            # update datatype
            errors.append(newinfo.update_datatype(datatype))

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

                        # Finding if DNA sequence contains error
                        seq_error_cnt = 0
                        seq_error_list = []

                        seq_string = manage_unicode(seq_string)
                        for x in seq_string:  # x is every character of sequence
                            if not x.lower() in "acgtryswkmbdhvn-.":
                                seq_error_cnt += 1
                                seq_error_list.append(x)

                        if seq_error_cnt > 0:
                            warnings.append(
                                f"In table {table}, Illegal DNA character {seq_error_list} found in {gene} of {datatype} {df['id'][n]}"
                            )
                        elif seq_string.lower().strip() in ("nan", "na"):
                            warnings.append(
                                f"In table {table}, Sequence {df['id'][n]} {seq_string} detected as nan, removing it"
                            )
                        elif seq_error_cnt == 0:
                            # remove gaps for preventing BLAST error
                            errors.append(
                                newinfo.update_seq(
                                    gene, seq_string.replace("-", "").replace(".", "")
                                )
                            )

        # After successfully parsed this table, save it
        save.save_df(
            df,
            f"{path.out_db}/Saved_{'.'.join(table.split('/')[-1].split('.')[:-1])}.{opt.tableformat}",
            fmt=opt.tableformat,
        )

    # Remove non-errors from errors
    errors = [err for err in errors if not (err is None)]

    return funinfo_dict, GenMine_flag, warnings, errors


def db_input(funinfo_dict, opt, path) -> list:
    # Get DB input
    logging.info(f"Input DB list: {opt.db}")

    (
        funinfo_dict,
        GenMine_flag,
        warnings,
        errors,
    ) = input_table(
        funinfo_dict=funinfo_dict,
        path=path,
        opt=opt,
        table_list=opt.db,
        datatype="db",
    )

    for warning in sorted(list(warnings)):
        logging.warning(warning)

    if len(errors) > 0:
        for error in sorted(list(set(errors))):
            logging.error(error)

        logging.error(
            f"FunVIP terminated because error found during input validation. Please check [ERROR] in log.txt"
        )
        raise Exception
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

    return funinfo_dict, GenMine_flag


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

    funinfo_dict, GenMine_flag, table_warnings, table_errors = input_table(
        funinfo_dict=funinfo_dict,
        path=path,
        opt=opt,
        table_list=query_table,
        datatype="query",
    )
    funinfo_dict, fasta_warnings, fasta_errors = input_fasta(
        path=path,
        opt=opt,
        fasta_list=query_fasta,
        funinfo_dict=funinfo_dict,
        datatype="query",
    )

    warnings = table_warnings + fasta_warnings
    errors = table_errors + fasta_errors

    for warning in sorted(list(warnings)):
        logging.warning(warning)

    if len(errors) > 0:
        for error in sorted(list(set(errors))):
            logging.error(error)

        logging.error(
            f"FunVIP terminated because error found during input validation. Please check [ERROR] in log.txt"
        )
        raise Exception

    # Save query initially in raw format, because gene has not been assigned yet
    for file in query_table:
        shutil.copy(f"{file}", f"{path.out_query}")

    for file in query_fasta:
        shutil.copy(f"{file}", f"{path.out_query}")

    logging.info(
        f"Total {len([funinfo_dict[key].datatype =='query' for key in funinfo_dict.keys()])} sequences parsed from query"
    )

    return funinfo_dict, GenMine_flag


# combined db and query input
def data_input(V, R, opt, path):
    funinfo_dict = {}

    # get database input
    funinfo_dict, GenMine_flag_db = db_input(
        funinfo_dict=funinfo_dict, opt=opt, path=path
    )

    # get query input
    funinfo_dict, GenMine_flag_query = query_input(
        funinfo_dict=funinfo_dict, opt=opt, path=path
    )

    if GenMine_flag_db != 0 or GenMine_flag_query != 0:
        GenMine_flag = 1
    else:
        GenMine_flag = 0

    # combine all data
    V.list_FI = [funinfo_dict[key] for key in funinfo_dict]

    # hashing data for safety in tree analysis
    V.list_FI = hash_funinfo_list(V.list_FI)

    # make hash dict
    for FI in V.list_FI:
        V.dict_hash_FI[FI.hash] = FI

    # If all genes are empty, indicate them
    for FI in V.list_FI:
        flag_available_gene = 0
        for gene in FI.seq:
            if FI.seq[gene] != "":
                flag_available_gene = 1

        if flag_available_gene == 0:
            FI.issues.add("noseq")

    return V, R, opt, GenMine_flag


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
        df_out, f"{path.out_query}/Saved_query.{opt.tableformat}", fmt=opt.tableformat
    )

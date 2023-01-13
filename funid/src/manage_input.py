import os
import sys
import shutil
import re
import subprocess
import json
import pandas as pd
import numpy as np
import logging
from unidecode import unidecode
from Bio import SeqIO
from funid.src.tool import (
    initialize_path,
    get_genus_species,
    get_accession,
    manage_unicode,
)
from funid.src.logics import isnewicklegal, isuniquecolumn
from funid.src.hasher import decode, newick_legal, hash_funinfo_list

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

        self.original_accession = ""  # original accession, can be newick illegal
        self.accession = ""  #  newick illegal characters removed
        self.hash = ""  # hash : HSXXHE
        self.description = ""  # full description from fasta
        self.ori_genus = ""  # original genus
        self.genus = ""  # final genus
        self.ori_species = ""  # original species
        self.bygene_species = {}  # species name designated by gene
        self.final_species = ""  # final species designated by concatenated analysis
        self.species_identifier = (
            0  # species identifier if multiple branches with same species exists
        )
        self.source = ""
        self.datatype = ""  # DB or Query
        self.section = ""  # original section
        self.adjusted_section = ""  # adjusted section by sectional clustering
        self.seq = {}
        self.unclassified_seq = []

    def update_seqrecord(self, seq, gene=None):
        self.description = seq.description
        self.genus, self.ori_species = get_genus_species(seq.description)

        if gene in self.seq:
            logging.error(f"More than 1 sequence for {gene} found for {self.accession}")
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
                    f"More than 1 sequence for {gene} found for {self.accession}"
                )
                raise Exception
            else:
                pass
        else:
            self.seq[gene] = seq

        self.bygene_species[gene] = self.ori_species

    def update_description(self, description):
        self.description = description

    def update_genus(self, genus):

        # Try to solve illegal unicode characters
        if pd.isnull(genus):
            genus = ""
        genus = genus.strip()
        genus = manage_unicode(genus)

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
        species = manage_unicode(species)

        # Check ambiguity
        if self.ori_species != "" and self.ori_species != species:
            logging.error(
                f"Colliding species info found for {self.original_accession}, {self.ori_species} and {species}"
            )
            raise Exception

        # Update original if should
        if self.ori_species == "":
            self.ori_species = species

    def update_species(self, gene, species):

        self.bygene_species[gene] = species

    def update_section(self, section):
        # Try to solve illegal unicode characters
        if pd.isnull(section):
            section = ""
        section = section.strip()
        section = manage_unicode(section)

        # Check ambiguity
        if self.section != "" and self.section != section:
            logging.error(
                f"Colliding section info found for {funinfo}, {self.section} and {section}"
            )
            raise Exception

        # Update section
        self.section = section

    def update_datatype(self, datatype):

        # Available datatypes : db, query
        if not (datatype in ("db", "query", "outgroup")):
            logging.error(f"{datatype} is not available datatype")
            raise Exception

        # Check ambiguity
        if self.datatype != "" and self.datatype != datatype:
            logging.error(
                f"Colliding section info found for {funinfo}, {self.datatype} and {datatype}"
            )

        self.datatype = datatype

    def update_accession(self, accession, regexs=None):

        if not regexs == None:
            accession = get_accession(accession, tuple(regexs))

        # if cannot find accession by regex
        if accession == "":
            accession = newick_legal(accession)
        accession = str(accession)
        self.original_accession = accession
        if not (isnewicklegal(accession)):
            for c in NEWICK_ILLEGAL:
                accession = accession.replace(c, "")
            accession = accession.replace(" ", "_")
        self.accession = accession

    def update_hash(self, n):
        self.hash = f"HS{n}HE"

    def __repr__(self):
        return f"FI: {self.accession}"

    def __hash__(self):
        return hash((self.original_accession, self.hash, self.description))

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        return (
            self.original_accession == other.original_accession
            and self.hash == other.hash
            and self.description == other.description
        )


# getting data input from fasta file
def input_fasta(path, opt, fasta_list, datatype):

    # initialize path to use function "get_genus_species"
    initialize_path(path)
    funinfo_list = []

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

        ## Testing
        fasta_file = SeqIO.parse(file, "fasta")
        for seq in fasta_file:
            newinfo = Funinfo()
            newinfo.update_seqrecord(seq)
            newinfo.update_datatype(datatype)
            newinfo.update_section("")  # because not sectioin not designated yet
            # Accession by regex match
            newinfo.update_accession(seq.description, regexs=opt.regex)
            if get_genus_species(seq.description)[0] != "":
                newinfo.update_genus(get_genus_species(seq.description)[0])

            if get_genus_species(seq.description)[1] != "":
                newinfo.update_ori_species(get_genus_species(seq.description)[1])

            tmp_list.append(newinfo)

        funinfo_list += tmp_list

        try:
            fasta_file = SeqIO.parse(file, "fasta")
            for seq in fasta_file:
                newinfo = Funinfo()
                newinfo.update_seqrecord(seq)
                newinfo.update_datatype(datatype)
                newinfo.update_section("")  # because not sectioin not designated yet
                # Accession by regex match
                newinfo.update_accession(seq.description, regexs=opt.regex)
                if get_genus_species(seq.description)[0] != "":
                    newinfo.update_genus(get_genus_species(seq.description)[0])

                if get_genus_species(seq.description)[1] != "":
                    newinfo.update_ori_species(get_genus_species(seq.description)[1])

                tmp_list.append(newinfo)

            funinfo_list += tmp_list
        except:
            logging.warning(f"{file} is not a valid fasta file skipping")

    return funinfo_list


# getting datafile from excel or tabular file
def input_table(path, option, table_list, datatype):

    string_error = 0

    initialize_path(path)  # this one is ugly
    funinfo_dict = {}
    df_list = []

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

        # Lower case
        df.columns = df.columns.str.lower()
        df_list.append(df)

        # Check if accession column exists and unique
        flag_accession = isuniquecolumn(
            list_column=list(df.columns), column="accession", table_name=table
        )

        # Check if genus column exists and unique
        check_none = True if datatype == "db" else False
        flag_genus = isuniquecolumn(
            list_column=list(df.columns),
            column="genus",
            table_name=table,
            check_none=check_none,
        )

        # Check if species column exists and unique
        check_none = True if datatype == "db" else False
        flag_species = isuniquecolumn(
            list_column=list(df.columns),
            column="species",
            table_name=table,
            check_none=check_none,
        )

        # Check if section column exists and unique
        check_none = True if datatype == "db" else False
        flag_section = isuniquecolumn(
            list_column=list(df.columns),
            column="section",
            table_name=table,
            check_none=check_none,
        )

        # Sequence column operations, download sequences with GenMine
        if 1:  # If download on/off option added, change this part
            download_dict = {}  # for downloaded sequences
            download_set = set()
            # 1 letter + 5 digit regex should be last, because they overlap with 2 letter + 6 digit ids
            regex_genbank = r"(([A-Z]{1}[0-9]{5})(\.[0-9]{1}){0,1})|(([A-Z]{2}[\_]{0,1}[0-9]{6}){1}([\.][0-9]){0,1})"

            # if gene name were not designated by user, use seq
            option.gene = list(set([gene.lower().strip() for gene in option.gene]))

            # should be imported after initialize to prevent error
            from .tool import mkdir

            # find all NCBI accessions in seq
            for gene in option.gene:
                if isuniquecolumn(
                    list_column=df.columns,
                    column=gene,
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
                cmd = f"GenMine -c {path.GenMine}/Accessions.txt -o {path.GenMine} -e {option.email}"
                subprocess.call(cmd)

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
                    for n, _ in enumerate(df["accession"]):
                        for gene in option.gene:
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

                else:
                    logging.warning(
                        f"None of the GenMine results were succesfully parsed"
                    )

        # Generate funinfo by each row
        for n, acc in enumerate(df["accession"]):
            # Check if accession is unique
            # Remove non-unicode first
            new_acc = True
            df["accession"][n] = manage_unicode(
                str(df["accession"][n]), column="accssion", row=n
            )
            # Generate funinfo for each accession
            if df["accession"][n] in funinfo_dict:
                newinfo = funinfo_dict[df["accession"][n]]
                new_acc = False
                logging.warning(f"Duplicate accession {df['accession'][n]} found!")
            else:
                funinfo_dict[df["accession"][n]] = Funinfo()
                newinfo = funinfo_dict[df["accession"][n]]
                newinfo.update_accession(df["accession"][n])

            # if flag_genus is true, try to parse genus
            if flag_genus is True:
                newinfo.update_genus(df["genus"][n])

            # if flag_species is true, try to parse species
            if flag_species is True:
                newinfo.update_ori_species(df["species"][n])

            # if flag_section is true, try to parse section
            if flag_section is True:
                newinfo.update_section(df["section"][n])

            # update datatype
            newinfo.update_datatype(datatype)

            # parse each of the genes
            # For each of the gene
            for gene in option.gene:
                seq_error = 0
                if gene in df.columns:
                    if (
                        not (pd.isna(df[gene][n])) or str(df[gene][n]).strip() == ""
                    ):  # skip blank sequences
                        if df[gene][n].startswith(
                            ">"
                        ):  # automatically rearrange for total fasta style
                            seq_string = "".join(df[gene][n].split("\n")[1:])
                        else:
                            seq_string = df[gene][n]

                        # Finding if sequence contains error
                        seq_error_cnt = 0
                        seq_error_list = []
                        for x in seq_string:  # x is every character of sequence
                            if not x.lower() in "acgtryswkmbdhvn-.":
                                seq_error_cnt += 1
                                seq_error_list.append(x)

                        if seq_error_cnt > 0:
                            logging.warning(
                                f'Illegal DNA character {seq_error_list} found in {gene} of DB {df["accession"][n]}'
                            )
                        elif seq_error_cnt == 0:
                            # remove gaps for preventing BLAST error
                            newinfo.update_seq(
                                gene, seq_string.replace("-", "").replace(".", "")
                            )

    # make it to list at last
    list_funinfo = [funinfo_dict[x] for x in funinfo_dict]

    return list_funinfo, df_list


def db_input(option, path) -> list:

    # Get DB input
    logging.info(f"Input DB list: {option.db}")
    db_namelist = [db.split("/")[-1] for db in option.db]

    funinfo_list, df_list = input_table(path, option, option.db, "db")
    print(funinfo_list)

    # do it after save_db option enabled
    for n, db in enumerate(db_namelist):
        df_list[n].to_excel(
            f"{path.out_db}/Saved_{'.'.join(db.split('.')[:-1])}.xlsx", index=False
        )

    # validate dataset
    # if only one section exists, outgroup cannot work
    section_set = set(funinfo.section for funinfo in funinfo_list)
    section_set.discard("")
    if len(section_set) <= 1:
        logging.error(
            f"Only {len(section_set)} detected : {section_set}. Please add outgroup sequences"
        )
        raise Exception

    # check if minimum outgroup number count exceeds minimum section
    section_cnt_dict = {x: 0 for x in section_set}
    for funinfo in funinfo_list:
        if type(funinfo.section) is str:
            if funinfo.section in section_set:
                section_cnt_dict[funinfo.section] += 1

    if any(section_cnt_dict[x] < option.maxoutgroup for x in section_cnt_dict):
        logging.warning(
            f"Sequences in database of some section has lower number than MINIMUM_OUTGROUP_COUNT. It may cause error when outgroup selection, or may select not most appropriate outgroup to section. Please lower number of MINIMUM_OUTGROUP_COUNT in option or add more sequences to these sections"
        )

    return funinfo_list


def query_input(option, path):

    query_fasta = [
        file
        for file in option.query
        if any(file.endswith(x) for x in (".fa", ".fna", ".fas", ".fasta", ".txt"))
    ]
    query_excel = [file for file in option.query if (file.endswith(".xlsx"))]

    query_list = []
    query_list += input_fasta(path, option, query_fasta, "query")
    query_list += input_table(path, option, query_excel, "query")[0]

    for file in query_fasta:
        shutil.copy(f"{file}", f"{path.out_query}/{file}")

    logging.info(f"Total {len(query_list)} sequences parsed")

    return query_list, option


# combined db and query input
def data_input(V, R, opt, path):
    # get database input
    db_funinfo_list = db_input(opt, path)
    # get query input
    query_funinfo_list, opt = query_input(opt, path)
    # combine all data
    V.list_FI = db_funinfo_list + query_funinfo_list
    # hashing data for safety in tree analysis
    V.list_FI = hash_funinfo_list(V.list_FI)

    # make hash dict
    for FI in V.list_FI:
        V.dict_hash_FI[FI.hash] = FI

    # update report
    # R.statistics.update_input_statistics(V, opt)

    return V, R, opt


def save(list_funinfo, path, option):
    def save_originalfasta(list_info, path, filename):
        with open(f"{path}/{filename}", "w") as fp:
            for info in list_info:
                fp.write(f">{info.description}\n{info.seq}\n")

    def save_excel(list_info, path, filename):
        dict_excel = {
            "hash": [],
            "accession": [],
            "genus": [],
            "species": [],
            "source": [],
            "datatype": [],
            "section": [],
            "adjusted_section": [],
        }

        seq_set = set()
        for info in list_info:
            for gene in info.seq:
                seq_set.add(gene)

        for gene in seq_set:
            dict_excel[gene] = []

        for info in list_info:
            dict_excel["hash"].append(info.hash)
            dict_excel["accession"].append(info.original_accession)
            dict_excel["genus"].append(info.genus)
            dict_excel["species"].append(info.ori_species)
            dict_excel["source"].append(info.source)
            dict_excel["datatype"].append(info.datatype)
            dict_excel["section"].append(info.section)
            dict_excel["adjusted_section"].append(info.adjusted_section)
            for gene in seq_set:
                if gene in info.seq:
                    dict_excel[gene].append(info.seq[gene])
                else:
                    dict_excel[gene].append("")

        df = pd.DataFrame(dict_excel)
        df.to_excel(f"{path}/{filename}", index=False)

    save_excel(list_funinfo, path.data, f"{option.runname}_Section Assignment.xlsx")

    # Save by source
    origin_set = set()

    for funinfo in list_funinfo:
        origin_set.add((funinfo.source, funinfo.datatype))

    for origin in origin_set:

        # set path
        if origin[1] in ["DB", "Query", "Outgroup"]:
            outpath = path.data
        else:
            logging.info(origin)
            logging.error("Wrong datatype")
            raise Exception

        tmp_list = []
        for funinfo in list_funinfo:
            if funinfo.datatype == origin[1]:
                if funinfo.source == origin[0]:
                    tmp_list.append(funinfo)


def save_fasta(list_funinfo, gene, filename, by="accession"):

    list_funinfo = list(set(list_funinfo))  # remove ambiguous seqs

    with open(f"{filename}", "w") as fp:
        if gene == "unclassified":  # for unclassified query
            flag = 0
            for info in list_funinfo:
                for n, seq in enumerate(info.unclassified_seq):
                    if by == "hash":
                        fp.write(f">{info.hash}_{n}\n{seq}\n")
                    else:
                        fp.write(f">{info.accession}_{n}\n{seq}\n")
                    flag = 1

        else:
            flag = 0
            for info in list_funinfo:
                if gene in info.seq:
                    if gene in info.seq:
                        if by == "hash":
                            fp.write(f">{info.hash}\n{info.seq[gene]}\n")
                        else:
                            fp.write(f">{info.accession}\n{info.seq[gene]}\n")
                        flag = 1

    # returns 1 if meaningful sequence exists
    return flag


def save_originalfasta(list_info, path, filename):
    with open(f"{path}/{filename}", "w") as fp:
        for info in list_info:
            fp.write(f">{info.description}\n{info.seq}\n")


def save_fastabysection(list_funinfo, path, option, add="Reference", outgroup=False):

    outpath = path.data
    set_section = set()

    for section in set_section:
        tmp_list = []
        for funinfo in list_funinfo:
            if funinfo.adjusted_section == section:
                tmp_list.append(funinfo)

        save_fasta(tmp_list, outpath, f"{add}_{section}.fasta")


# Save tree file to designated path, and decode it
def save_tree(out, hash_dict, hash_file_path, decoded_file_path):

    # print(out)
    file = out.split("/")[-1]
    shutil.move(out, hash_file_path)
    decode(
        hash_dict=hash_dict,
        file=hash_file_path,
        out=decoded_file_path,
    )


def save_mergedfasta(fasta_list, out_path):

    out_fasta_list = []
    for fasta in fasta_list:
        out_fasta_list += list(SeqIO.parse(fasta, "fasta"))

    SeqIO.write(out_fasta_list, out_path, "fasta")


# save dataframe
def save_df(df, out, fmt="csv"):

    if fmt == "csv" or fmt == "tsv":
        df.to_csv(out, index=False)
    elif fmt == "xlsx" or fmt == "excel":
        df.to_excel(out, index=False)
    elif fmt == "parquet":
        df.to_parquet(out, index=False)
    elif fmt == "feather" or fmt == "ftr":
        df.to_feather(out, index=False)
    else:
        logging.warning(
            f"Not appropriate format entered for matrix format, using csv as default"
        )
        df.to_csv(out, index=False)

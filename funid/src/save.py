# For save and load functions
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
from funid.src.logics import isnewicklegal
from funid.src.hasher import decode, newick_legal, hash_funinfo_list
import shelve

# Saving functions should be in main branch because of globals() function
# try move it to other place by sending globals() as variable
# In future, try selectively save to reduce datasize
# Session saving function
def save_session(opt, path, global_var: dict, var: dict) -> None:

    # if opt.save_run is True:
    if 1:  # ADD saving options in further developmental stage
        save = shelve.open(path.save, "n")
        for key in global_var:
            try:
                save[key] = global_var[key]
                if opt.verbose >= 3:
                    logging.debug(f"Saved {key}")
            except:
                if opt.verbose >= 3:
                    logging.debug(f"Failed shelving: {key}")
        logging.info(f"Saved current session")
        save.close()


# Session loading function
def load_session(opt, global_var: dict, savefile: str) -> None:
    save = shelve.open(savefile)
    for key in save:
        try:
            global_var[key] = save[key]
            if opt.verbose >= 3:
                logging.debug(f"Loading {key}")
        except:
            if opt.verbose >= 3:
                logging.debug(f"Failed Loading {key}")

    logging.info(f"Loaded previous session")
    save.close()


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
        if origin[1] in ["db", "query", "outgroup"]:
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

    """
    for funinfo in list_funinfo:
        set_section.add(funinfo.adjusted_section)
        if funinfo.adjusted_section != str:
            print(funinfo)
    """

    # print(set_section)

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
            f"[Warning] Not appropriate format entered for matrix format, using csv as default"
        )
        df.to_csv(out, index=False)

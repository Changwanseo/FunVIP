import time, re
from Bio import SeqIO
from Bio.Seq import Seq

# from .logger import Mes
import logging
from functools import lru_cache
import sys, os
import shutil
import platform

# Maybe we should move this to /funid/src/toolbox/ and split these to multiple files
def initialize_path(path):
    global genus_file
    genus_file = path.genusdb


def mkdir(path):
    if os.path.exists(path) == False:
        os.makedirs(path)


def union_funinfo_list(funinfo_list1, funinfo_list2):

    hash_list = [x.hash for x in funinfo_list1]
    for funinfo in funinfo_list2:
        if not funinfo.hash in hash_list:
            funinfo_list1.append(funinfo)

    return funinfo_list1


# if string is not ascii, automatically solve it with unicode library
def manage_unicode(string, column="", row=""):

    try:
        string.encode("ascii")
        return string
    except:
        logging.warning(
            f"Illegal unicode character found in ID/Accession column, {row}th row : {string}. Trying flexible solve"
        )
        try:
            return unidecode(string)
        except:
            logging.error(
                f"Flexible solve failed to {string}. Please change the cell with available ascii strings"
            )
            raise Exception


@lru_cache(maxsize=10000)
def get_genus_species(
    string,
    endwords=[
        "small",
        "18S",
        "ribosomal",
        "internal",
        "gene",
        "ITS",
        "5.8S",
        "voucher",
        "strain",
        "beta",
        "tubulin",
    ],
    genus_list=None,
):

    return_genus = ""
    return_species = ""

    # en for enumeratable object (splited string)
    if genus_list is None:
        with open(genus_file, "r") as f:
            genus_list = f.read().splitlines()

    en = string.replace(" ", "_").split("_")

    for genus in genus_list:
        for n, i in enumerate(en):
            if i == genus:  # or i == genus[0] or i == genus[0]+".":
                return_genus = i
                try:
                    if en[n + 2] in ["var", "var.", "f", "f.", "nom", "nom."]:
                        try:
                            en[n + 3]
                            if not (en[n + 3] in endwords):
                                return_species = " ".join(
                                    [en[n + 1], en[n + 2], en[n + 3]]
                                )
                            else:
                                return_species = " ".join([en[n + 1], en[n + 2]])
                        except:
                            return_species = " ".join([en[n + 1], en[n + 2]])

                    else:
                        if en[n + 1] in ["aff", "aff.", "cf", "cf."]:
                            try:
                                if not (en[n + 2] in endwords):
                                    return_species = " ".join([en[n + 1], en[n + 2]])
                                else:
                                    return_species = en[n + 1]
                            except:
                                return_species = en[n + 1]

                        elif en[n + 1] in ["sp", "sp."]:
                            try:
                                if not (en[n + 2] in endwords):
                                    try:
                                        int(en[n + 2])
                                        return_species = " ".join(
                                            [en[n + 1], en[n + 2]]
                                        )
                                    except:
                                        return_species = en[n + 1]
                                else:
                                    return_species = en[n + 1]
                            except:
                                return_species = en[n + 1]

                        else:
                            if not (en[n + 1] in endwords):
                                return_species = en[n + 1]
                            else:
                                return_species = "NaN"

                except:
                    try:
                        if en[n + 1] in ["sp", "sp.", "aff", "aff."]:
                            try:
                                en[n + 2]
                                if not (en[n + 2] in endwords):
                                    return_species = " ".join([en[n + 1], en[n + 2]])
                                else:
                                    return_species = en[n + 1]
                            except:
                                return_species = en[n + 1]
                        else:
                            if not (en[n + 1] in endwords):
                                return_species = en[n + 1]
                            else:
                                return_species = "NaN"
                    except:
                        return_species = "NaN"

    return return_genus, return_species


# get id from string by regex match
# return message as error, not print it
@lru_cache(maxsize=10000)
def get_id(string, id_list):
    id_set = set()
    id_ = ""
    for regex in id_list:
        if re.search(regex, string):
            id_set.add(re.search(regex, string).group(0))
            # longest match as a best match
            if len(re.search(regex, string).group(0)) > len(id_):
                id_ = re.search(regex, string).group(0)
        else:
            pass

    if len(id_set) >= 2:
        logging.warning(
            f"[Warning] Ambiguous regex match found. {id_} selected among {id_set}"
        )

    if id_ == "":
        logging.warning(
            f"[Warning] Cannot find regex match from {string}, using default fasta name"
        )
        id_ = string

    id_ = str(id_)
    return id_


# select FI with specific datatype from FI_list
def select(funinfo_list, datatype):

    tmp_list = []
    for funinfo in funinfo_list:
        if funinfo.datatype == datatype:
            tmp_list.append(funinfo)

    return tmp_list


def excel_sheetname_legal(string):
    newick_illegal = ["'", "[", "]", ":", "*", "?", "/", "\\"]
    for i in newick_illegal:
        string = string.replace(i, "")

    string = string.replace("  ", " ")
    string = string.replace(" ", "_")

    return string


# moves newick files generated in tree building process to appropriate location
def cleanup_tree(path):

    # move result files to each of the folders
    mkdir(f"{path.out_tree}/hash")
    hash_files = [
        f
        for f in os.listdir(f"{path.out_tree}")
        if f.startswith("hash_") and f.endswith(".nwk")
    ]
    for file in hash_files:
        try:
            shutil.move(f"{path.out_tree}/{file}", f"{path.out_tree}/hash/{file}")
        except:
            os.remove(f"{path.out_tree}/hash/{file}")
            shutil.move(f"{path.out_tree}/{file}", f"{path.out_tree}/hash/{file}")


# moves svg files generated in CAT-V to appropriate location
def cleanup_tree_image(path):

    # move result files to each of the folders
    mkdir(f"{path.out_tree}/hash")
    hash_files = [
        f
        for f in os.listdir(f"{path.out_tree}")
        if f.startswith("hash_") and f.endswith(".svg")
    ]
    for file in hash_files:
        try:
            shutil.move(f"{path.out_tree}/{file}", f"{path.out_tree}/hash/{file}")
        except:
            os.remove(f"{path.out_tree}/hash/{file}")
            shutil.move(f"{path.out_tree}/{file}", f"{path.out_tree}/hash/{file}")

    mkdir(f"{path.out_tree}/original")
    original_files = [
        f for f in os.listdir(f"{path.out_tree}") if f.endswith("_original.svg")
    ]
    for file in original_files:
        try:
            shutil.move(f"{path.out_tree}/{file}", f"{path.out_tree}/original/{file}")
        except:
            os.remove(f"{path.out_tree}/original/{file}")
            shutil.move(f"{path.out_tree}/{file}", f"{path.out_tree}/original/{file}")


# Change step string into numbered step
def index_step(step):
    # Declaring available steps
    step_list = [
        "setup",
        "search",
        "cluster",
        "align",
        "trim",
        "concatenate",
        "modeltest",
        "tree",
        "visualize",
        "report",
    ]

    try:
        return step_list.index(step.lower().strip())
    except:
        logging.error(f"DEVELOPMENTAL ERROR, INVALID STEP {step} USED")
        raise Exception


def check_avx():
    return platform.machine() == "x86_64"

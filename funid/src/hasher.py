import os, shutil
import re
from Bio import SeqIO
import copy
import pandas as pd


def newick_legal(string: str) -> str:
    """
    Remove all newick illegal strings
    """
    newick_illegal = ["(", ")", "{", "}", "[", "]", ":", ";", "'", '"', ",", "."]
    for i in newick_illegal:
        string = string.replace(i, "")

    string = string.replace("  ", " ")
    string = string.replace(" ", "_")

    return str(string)


def encode(funinfo_list: list, newick: bool = False) -> dict:
    """
    Encode funinfo_list and return hash dict
    """
    hash_dict = {}

    if newick is False:
        for funinfo in funinfo_list:
            hash_dict[funinfo.hash] = f"{funinfo.accession}"
    else:
        for funinfo in funinfo_list:
            try:
                hash_dict[
                    funinfo.hash
                ] = f"{funinfo.accession}_{funinfo.genus}_{funinfo.species}"
            except:
                hash_dict[
                    funinfo.hash
                ] = f"{funinfo.accession}_{funinfo.genus}_{funinfo.ori_species}"

    return hash_dict


def decode(hash_dict: dict, file: str, out: str, newick: bool = True) -> None:
    """
    Decode given file with given hash_dict
    """

    with open(file, "rt") as fp:
        line = fp.read()
        if newick == True:
            hash_dict = dict(
                (re.escape(k), newick_legal(v)) for k, v in hash_dict.items()
            )
        else:
            hash_dict = dict((re.escape(k), v) for k, v in hash_dict.items())

        pattern = re.compile("|".join(hash_dict.keys()))
        line = pattern.sub(lambda m: hash_dict[re.escape(m.group(0))], line)

        with open(out, "w") as fw:
            fw.write(line)


def decode_df(hash_dict: dict, df: pd.DataFrame) -> pd.DataFrame:
    """
    Decode given dataframe with given hash_dict
    """

    hash_dict = dict((re.escape(k), v) for k, v in hash_dict.items())

    df_return = copy.deepcopy(df)

    for column in df.columns:
        df_return[column] = df_return[column].map(lambda x: hash_dict.get(x, x))

    return df_return


def hasher(funinfo_list: list, path, option, outgroup: bool = False):
    """
    # main pipeline
    Hash all
    """
    section_result = {}
    for section in section_list:
        tmp_funinfo_list = []
        for funinfo in funinfo_list:
            if funinfo.adjusted_section == section:
                tmp_funinfo_list.append(funinfo)

        section_result[section] = (
            f"Hashed_{section}.fasta",
            encode(tmp_funinfo_list, f"{path.out_hash}/Hashed_{section}.fasta"),
        )

    return section_result


def hash_funinfo_list(list_funinfo: list) -> list:
    """
    Generate hash numbers
    """
    for n, funinfo in enumerate(list_funinfo):
        funinfo.update_hash(n)

    return list_funinfo

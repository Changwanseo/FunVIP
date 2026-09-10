import os, shutil
import re
from Bio import SeqIO
import copy
import pandas as pd
import sys


# Remove all newick illegal strings
def newick_legal(string: str) -> str:
    newick_illegal = ["(", ")", "{", "}", "[", "]", ":", ";", "'", '"', ",", "."]
    for i in newick_illegal:
        string = string.replace(i, "")

    string = string.replace("  ", " ")
    string = string.replace(" ", "_")

    return str(string)


# Fix "&" sign appropriate for svg
def svg_legal(string: str) -> str:
    return string.replace("&", "&amp;")


# Encode funinfo_list and return hash dict
def encode(funinfo_list: list, newick: bool = False) -> dict:
    hash_dict = {}

    if newick is False:
        return {funinfo.hash: f"{funinfo.id}" for funinfo in funinfo_list}

        """
        for funinfo in funinfo_list:
            hash_dict[funinfo.hash] = f"{funinfo.id}"
        """
    else:
        return {
            funinfo.hash: f"{funinfo.id}_{funinfo.genus}_{getattr(funinfo, 'species', funinfo.ori_species)}"
            for funinfo in funinfo_list
        }

        """
        for funinfo in funinfo_list:
            try:
                hash_dict[
                    funinfo.hash
                ] = f"{funinfo.id}_{funinfo.genus}_{funinfo.species}"
            except:
                hash_dict[
                    funinfo.hash
                ] = f"{funinfo.id}_{funinfo.genus}_{funinfo.ori_species}"
        """

    # return hash_dict


# Decode given file with given hash_dict


def decode(
    hash_dict: dict, file: str, out: str, newick: bool = True, svg: bool = False
) -> None:
    with open(file, "rt", encoding="utf-8") as fp:
        content = fp.read()

    if newick and svg:
        lookup = {k: svg_legal(newick_legal(v)) for k, v in hash_dict.items()}
    elif newick:
        lookup = {k: newick_legal(v) for k, v in hash_dict.items()}
    elif svg:
        lookup = {k: svg_legal(v) for k, v in hash_dict.items()}
    else:
        lookup = dict(hash_dict)

    # Every hash has the fixed, self-delimiting form HS<number>HE, so match that
    # single token and look each occurrence up, instead of compiling an
    # N-alternative regex ("|".join of every hash) and re-escaping each match --
    # the former is O(text * N_hashes) and a hotspot at metabarcoding scale.
    # A HS<n>HE token not present in the dict is left unchanged.
    pattern = re.compile(r"HS\d+HE")
    decoded_content = pattern.sub(lambda m: lookup.get(m.group(0), m.group(0)), content)

    with open(out, "w", encoding="utf-8") as fw:
        fw.write(decoded_content)


"""
def decode(hash_dict: dict, file: str, out: str, newick: bool = True) -> None:
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
"""


# Decode given dataframe with given hash_dict
def decode_df(hash_dict: dict, df: pd.DataFrame) -> pd.DataFrame:
    # hashes are HS<number>HE (no regex metacharacters), so no escaping is needed;
    # each cell is either exactly a hash to replace or is left as-is.
    df_return = df.copy()

    for column in df.columns:
        df_return[column] = df_return[column].map(lambda x: hash_dict.get(x, x))

    return df_return


def hasher(funinfo_list: list, path, option, outgroup: bool = False):
    """
    # main pipeline
    Hash all
    """
    group_result = {}
    for group in group_list:
        tmp_funinfo_list = []
        for funinfo in funinfo_list:
            if funinfo.adjusted_group == group:
                tmp_funinfo_list.append(funinfo)

        group_result[group] = (
            f"Hashed_{group}.fasta",
            encode(tmp_funinfo_list, f"{path.out_hash}/Hashed_{group}.fasta"),
        )

    return group_result


def hash_funinfo_list(list_FI: list) -> list:
    """
    Generate hash numbers
    """

    for n, FI in enumerate(list_FI):
        FI.update_hash(n)

        # To reduce memory

        # sys.intern(FI.hash)
        # sys.intern(FI.id)

    return list_FI

from Bio import Entrez
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
import time
from time import sleep
import pickle
import json
import re
import logging
from datetime import datetime
from itertools import repeat
import xmltodict
import os, sys, subprocess
import random
import pickle
import pandas as pd
import shutil

# from .logger import Mes


def Downloadbyacclist(Email, list_ID, out):

    if len(list_ID) == 0:
        return []

    path_tmp = "tmp_download"

    # print(list_ID)

    def xml2dict(record):
        dict_type = xmltodict.parse(record)
        json_type = json.dumps(dict_type, indent=4)
        dict2_type = json.loads(json_type)
        return dict2_type

    Entrez.email = Email

    logging.info(f"Number of IDs: {len(list_ID)}")

    # Parse
    cnt = 0
    cut = 50  # parse by 50 sequences
    cnt_all = len(list_ID)

    start_time = time.time()

    record_list = []

    for i in range(int((len(list_ID) - 1) / cut) + 1):

        # last chunk
        if i * cut + cut > len(list_ID):
            # Mes(i)
            ID_string = ",".join(list_ID[i * cut :])
            sleep(0.3)
            cnt += len(list_ID[i * cut :])
            logging.info(f"{cnt}/{cnt_all} {100*cnt/cnt_all}% {time.time()-start_time}s")

            try:
                print(ID_string)
                handle = Entrez.efetch(
                    db="nucleotide", id=ID_string, rettype="gb", retmode="xml"
                )
            except:  # Retry
                logging.info("Requesting again...")
                try:
                    sleep(10)
                    handle = Entrez.efetch(
                        db="nucleotide", id=ID_string, rettype="gb", retmode="xml"
                    )
                except:
                    sleep(100)
                    handle = Entrez.efetch(
                        db="nucleotide", id=ID_string, rettype="gb", retmode="xml"
                    )

            pre_record = handle.read()
            json_record = xml2dict(pre_record)
            tmp_record_list = []

            # If only one result available
            if type(json_record["GBSet"]["GBSeq"]) is dict:
                json_record["GBSet"]["GBSeq"] = [json_record["GBSet"]["GBSeq"]]

            if len(list_ID[i * cut :]) != 1:
                for record in json_record["GBSet"]["GBSeq"]:
                    print(record)
                    record_list.append(record)
                    tmp_record_list.append(record)
            else:
                record = json_record["GBSet"]["GBSeq"]
                record_list.append(record)
                tmp_record_list.append(record)

            try:
                with open(f"./{path_tmp}/{i}", "wb") as f:
                    pickle.dump(tmp_record_list, f)
            except:
                logging.error("Saving Error")
                raise Exception

        # non last chunk
        else:
            # Mes(i)
            ID_string = ",".join(list_ID[i * cut : i * cut + cut])
            sleep(0.3)
            cnt += cut
            logging.info(
                f"{cnt}/{cnt_all} {100*cnt/cnt_all}% {time.time()-start_time}s"
            )

            try:
                handle = Entrez.efetch(
                    db="nucleotide", id=ID_string, rettype="gb", retmode="xml"
                )
            except:  # Retry
                logging.info("Requesting again...")
                try:
                    sleep(10)
                    handle = Entrez.efetch(
                        db="nucleotide", id=ID_string, rettype="gb", retmode="xml"
                    )
                except:
                    sleep(100)
                    handle = Entrez.efetch(
                        db="nucleotide", id=ID_string, rettype="gb", retmode="xml"
                    )

            pre_record = handle.read()
            json_record = xml2dict(pre_record)
            tmp_record_list = []

            if len(list_ID[i * cut :]) != 1:
                for record in json_record["GBSet"]["GBSeq"]:
                    record_list.append(record)
                    tmp_record_list.append(record)
            else:
                record = json_record["GBSet"]["GBSeq"]
                record_list.append(record)
                tmp_record_list.append(record)

            try:
                with open(f"./{path_tmp}/{i}", "wb") as f:
                    pickle.dump(tmp_record_list, f)
            except:
                logging.error("Saving Error")
                raise Exception

    with open(out, "w") as fp:
        json_term = json.dump(record_list, fp, indent=4)

    # If success, remove tmp files
    tmp_file_list = [file for file in os.listdir(f"./{path_tmp}/")]
    for file in tmp_file_list:
        os.remove(f"./{path_tmp}/{file}")

    return record_list

from funid.src import save
from Bio import SeqIO
import os
import sys
import numpy as np
import logging
import re
import json


# Datasets are group of sequences with same genes for analysis, including query, database and outgroup
# Each of the dataset varialbe
class Dataset:
    def __init__(self, gene, group, list_qr, list_db, list_og):
        # Which gene is this dataset, gene name or concatenated
        self.gene = gene
        # Which group is this dataset
        self.group = group
        # query funinfo_list for this dataset
        self.list_qr_FI = list_qr
        # db funinfo_list for this dataset
        self.list_db_FI = list_db
        # outgroup funinfo_list for this dataset
        self.list_og_FI = list_og


# FunID all dataset variables bundle
class FunID_var:
    def __init__(self):
        # dataset class to manage datasets used for run
        # FI : funinfo
        # SR : search result
        # prefix c : concatenated
        # group : group
        # dict_A_B : dictionary, {A : B}
        # opt : list of option
        # qr : query
        # db : database
        # og : outgroup
        # rslt : result

        # Funinfo related variables
        # all funinfo_list, old funinfo_list
        self.list_FI = []
        # {funinfo.hash : funinfo}, old hash_dict
        self.dict_hash_FI = {}
        # {funinfo.hash : "{funinfo.id}_{funinfo.genus}_{funinfo.species}"}
        self.dict_hash_name = {}
        # {funinfo.hash : "{funinfo.id}"}
        self.dict_hash_id = {}

        # genus, gene - dataset
        # tuple of genera used for current run, old genus_list
        self.tup_genus = ()
        # list of gene found from db in current run
        self.list_db_gene = []
        # list of gene used for query in current run
        self.list_qr_gene = []
        # old group list
        self.list_group = []

        # Search result datasets
        # {gene : search df}, old df_dict
        self.dict_gene_SR = {}
        # concatenated search df dict, old concatenated_df
        self.cSR = None
        # {groupon : concatenated search df dict}, old concatenated_df in dictionary form by group
        # self.dict_group_cSR = {}

        # dataset dict - dict_dataset[group][gene] = class Dataset
        # using key concat for concatenated
        self.dict_dataset = {}

        # Running options
        # group clustering option
        self.opt_cluster = []
        self.opt_align = []
        self.opt_tree = []

        # Running result
        self.rslt_cluster = []
        self.rslt_align = []
        self.rslt_tree = []

        # multigene_list
        self.multigene_list = None

    def __repr__(self):
        out_dict = {}
        for key in self.dict_dataset:
            out_dict[key] = len(self.dict_dataset[key])
        return f"<< FunID_var object >>\nNumber of FI: {len(self.list_FI)}\nDB genes:{self.list_db_gene}\nQuery genes:{self.list_qr_gene}\ngroups:{self.list_group}\nDataset_dict:{json.dumps(out_dict, indent=2)}"

    def add_dataset(self, group, gene, list_qr, list_db, list_og):
        data = Dataset(group, gene, list_qr, list_db, list_og)
        if not (group) in self.dict_dataset:
            self.dict_dataset[group] = {}

        if not (gene) in self.dict_dataset[group]:
            self.dict_dataset[group][gene] = data
        else:
            logging.warning(f"Overriding dataset on {group} {gene}")
            self.dict_dataset[group][gene] = data

        # If gene is concatenated, update ori_species
        if gene == "concatenated":
            for FI in list_qr + list_db + list_og:
                FI.bygene_species[gene] = FI.ori_species

    def remove_dataset(self, group, gene):
        if not (group) in self.dict_dataset:
            logging.warning(
                f"Passing removing dataset of {group} {gene} because no {group} priorly exists"
            )
        elif not (gene) in self.dict_dataset[group]:
            logging.warning(
                f"Passing removing dataset of {group} {gene} because no {gene} priorly exists"
            )
        else:
            del self.dict_dataset[group][gene]
            # if all gene removed
            if len(self.dict_dataset[group]) == 0:
                del self.dict_dataset[group]

    def exist_dataset(self, group, gene):
        try:
            self.dict_dataset[group][gene]
            return True
        except:
            return False

    # generate dataset by group and gene
    def generate_dataset(self, opt):
        dict_funinfo = {}

        for group in self.list_group:
            logging.info(f"Generating dataset for {group}")
            dict_funinfo[group] = {}

            # For queryonly case
            if opt.queryonly is True:
                group_flag = False  # whether to run this group
                for gene in self.list_db_gene:
                    logging.debug(f"Searching {group} {gene} is good dataset")
                    list_qr = [
                        FI
                        for FI in self.list_FI
                        if (
                            gene in FI.seq
                            and FI.datatype == "query"
                            and FI.adjusted_group == group
                        )
                    ]

                    # do not manage db when query only mode and query does not exists
                    if len(list_qr) > 0:
                        group_flag = True

                if group_flag:  # if decided to run this group
                    logging.info(f"Decided to construct dataset on {group}")
                    for gene in self.list_db_gene:
                        list_qr = [
                            FI
                            for FI in self.list_FI
                            if (
                                gene in FI.seq
                                and FI.datatype == "query"
                                and FI.adjusted_group == group
                            )
                        ]

                        list_db = [
                            FI
                            for FI in self.list_FI
                            if (
                                gene in FI.seq
                                and FI.datatype == "db"
                                and FI.adjusted_group == group
                            )
                        ]
                        self.add_dataset(group, gene, list_qr, list_db, [])

                else:
                    logging.warning(
                        f"group {group} did not passed dataset construction"
                    )

                # for concatenated
                list_qr = [
                    FI
                    for FI in self.list_FI
                    if (FI.datatype == "query" and FI.adjusted_group == group)
                ]

                # do not manage db when query only mode and query does not exists
                if len(list_qr) > 0:
                    list_db = [
                        FI
                        for FI in self.list_FI
                        if (FI.datatype == "db" and FI.adjusted_group == group)
                    ]
                    self.add_dataset(group, "concatenated", list_qr, list_db, [])

            # For opt.queryonly is False -> run all dataset in database
            else:
                for gene in self.list_db_gene:
                    list_qr = [
                        FI
                        for FI in self.list_FI
                        if (
                            gene in FI.seq
                            and FI.datatype == "query"
                            and FI.adjusted_group == group
                        )
                    ]
                    list_db = [
                        FI
                        for FI in self.list_FI
                        if (
                            gene in FI.seq
                            and FI.datatype == "db"
                            and FI.adjusted_group == group
                        )
                    ]
                    self.add_dataset(group, gene, list_qr, list_db, [])

                # for concatenated
                list_qr = [
                    FI
                    for FI in self.list_FI
                    if (FI.datatype == "query" and FI.adjusted_group == group)
                ]
                list_db = [
                    FI
                    for FI in self.list_FI
                    if (FI.datatype == "db" and FI.adjusted_group == group)
                ]
                self.add_dataset(group, "concatenated", list_qr, list_db, [])

    # homogenize list_dataset and dicts from multiple results
    def homogenize_dataset(self):
        for FI in self.list_FI:
            if FI.hash in self.dict_hash_FI:
                h = FI.hash

                # final species
                if FI.final_species != self.dict_hash_FI[h].final_species:
                    if FI.final_species == "":
                        FI.final_species = self.dict_hash_FI[h].final_species
                    elif self.dict_hash_FI[h].final_species == "":
                        self.dict_hash_FI[h].final_species = FI.final_species
                    else:
                        logging.error(
                            f"DEVELOPMNETAL ERROR Both list_FI and dict_hash_FI have conflicting final species, {FI.final_species} and {self.dict_hash_FI[h].final_species}"
                        )
                        raise Exception

                # adjusted group
                if FI.adjusted_group != self.dict_hash_FI[h].adjusted_group:
                    if FI.adjusted_group == "" or FI.adjusted_group == "":
                        FI.adjusted_group = self.dict_hash_FI[h].adjusted_group
                    elif self.dict_hash_FI[h].adjusted_group == "":
                        self.dict_hash_FI[h].adjusted_group = FI.adjusted_group
                    else:
                        logging.error(
                            f"DEVELOPMENTAL ERROR Both list_FI and dict_hash_FI have conflicting final group, {FI.adjusted_group} and {self.dict_hash_FI[h].adjusted_group}, {FI}"
                        )
                        raise Exception

                elif (
                    FI.adjusted_group == ""
                    and self.dict_hash_FI[h].adjusted_group == ""
                ):
                    if FI.group != "":
                        FI.adjusted_group = FI.group
                        self.dict_hash_FI[h].adjusted_group = FI.group

                    elif self.dict_hash_FI[h].group != "":
                        FI.adjusted_group = self.dict_hash_FI[h].group
                        self.dict_hash_FI[h].adjusted_group = self.dict_hash_FI[h].group
                    else:
                        logging.error(f"{FI.id} does not have assigned group!")

                # bygene_species
                if FI.bygene_species != self.dict_hash_FI[h].bygene_species:
                    # if both are empty
                    if not (FI.bygene_species) and not (
                        self.dict_hash_FI[h].bygene_species
                    ):
                        pass
                    elif not (FI.bygene_species):
                        FI.bygene_species = self.dict_hash_FI[h].bygene_species
                    elif self.dict_hash_FI[h].bygene_species:
                        self.dict_hash_FI[h].bygene_species = FI.bygene_species
                    else:
                        logging.error(
                            f"DEVELOPMENTAL ERROR Both list_FI and dict_hash_FI have conflicting gene identification results, {FI.bygene_species} and {self.dict_hash_FI[h].bygene}"
                        )
                        raise Exception

    # Remove invalid dataset to be analyzed
    def remove_invalid_dataset(self):
        # collect remove list : removing after iterating
        list_remove = []
        for group in self.dict_dataset:
            for gene in self.dict_dataset[group]:
                if len(self.dict_dataset[group][gene].list_db_FI) == 0:
                    logging.warning(
                        f"Removing {gene} from analysis in group {group} because there are no corresponding sequences"
                    )
                    list_remove.append((group, gene))
                elif (
                    len(self.dict_dataset[group][gene].list_db_FI)
                    + len(self.dict_dataset[group][gene].list_qr_FI)
                    + len(self.dict_dataset[group][gene].list_og_FI)
                    < 4
                ):
                    logging.warning(
                        f"Removing {group} {gene} from downstream phylogenetic analysis because there are not enough sequences"
                    )
                    list_remove.append((group, gene))

        for x in list_remove:
            self.remove_dataset(*x)

    # save fasta for outgroup adjusted fasta
    def save_dataset(self, path, opt):
        for group in self.dict_dataset:
            for gene in self.dict_dataset[group]:
                if not (gene == "concatenated"):
                    if "concatenated" in self.dict_dataset[group]:
                        fasta_list = list(
                            set(
                                self.dict_dataset[group][gene].list_db_FI
                                + self.dict_dataset[group][gene].list_qr_FI
                                + self.dict_dataset[group][gene].list_og_FI
                                + self.dict_dataset[group]["concatenated"].list_og_FI
                            )
                        )
                    else:
                        fasta_list = (
                            self.dict_dataset[group][gene].list_db_FI
                            + self.dict_dataset[group][gene].list_qr_FI
                            + self.dict_dataset[group][gene].list_og_FI
                        )

                    save.save_fasta(
                        fasta_list,
                        gene,
                        f"{path.out_adjusted}/{opt.runname}_Adjusted_{group}_{gene}.fasta",
                        by="hash",
                    )

    # Validate if any multiple sequence alignment has no overlapping region
    def validate_alignments(self, path, opt):
        for group in self.dict_dataset:
            for gene in self.dict_dataset[group]:
                # Check if dataset exists
                pass

                # Check if alignment corresponding to dataset exists
                if not (
                    os.path.isfile(
                        f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta"
                    )
                ):
                    logger.warning(
                        f"Alignment file {path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta does not exists"
                    )

                else:
                    # If alignment exists, check if alignment does have overlapping regions

                    ## Parse alignment
                    seq_list = list(
                        SeqIO.parse(
                            f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta",
                            "fasta",
                        )
                    )

                    ## Transform alignment to vector
                    ## Change gap to 0, and other characters to 1

                    vectors = [
                        np.fromiter(
                            re.sub(r"[^0]", "1", re.sub(r"[\-]", "0", str(seq.seq))),
                            dtype=np.int,
                        )
                        for seq in seq_list
                    ]

                    ## Multiply vectors
                    vector_products = np.prod(np.vstack(vectors), axis=0)

                    ## If all value of vectors are zero, it means that all regions have at least one gap
                    ## Raise warning for this
                    if np.all((vector_products == 0)):
                        logging.warning(
                            f"Alignment for {group} {gene} does not have any overlapping regions! Please check alignment to find out if some of the sequences were from different regions"
                        )
                    else:
                        logging.info(f"Alignment for {group} {gene} passed validation")

                # If alignment corresponding to dataset does not exists, raise warning or error
                pass

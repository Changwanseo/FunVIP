from funid.src import manage_input
from Bio import SeqIO
import os
import sys
import numpy as np
import logging
import re
import json

# Each of the dataset varialbe
class Dataset:
    def __init__(self, gene, sect, list_qr, list_db, list_og):

        # Which gene is this dataset, gene name or concatenated
        self.gene = gene
        # Which section is this dataset
        self.sect = sect
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
        # sect : section
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
        # {funinfo.hash : "{funinfo.accession}_{funinfo.genus}_{funinfo.species}"}
        self.dict_hash_name = {}
        # {funinfo.hash : "{funinfo.accession}"}
        self.dict_hash_acc = {}

        # genus, gene - dataset
        # tuple of genera used for current run, old genus_list
        self.tup_genus = ()
        # list of gene found from db in current run
        self.list_db_gene = []
        # list of gene used for query in current run
        self.list_qr_gene = []
        # old section list
        self.list_sect = []

        # Search result datasets
        # {gene : search df}, old df_dict
        self.dict_gene_SR = {}
        # concatenated search df dict, old concatenated_df
        self.cSR = None
        # {secton : concatenated search df dict}, old concatenated_df in dictionary form by section
        # self.dict_sect_cSR = {}

        # dataset dict - dict_dataset[section][dict] = class Dataset
        # using key concat for concatenated
        self.dict_dataset = {}

        # Running options
        # sectional clustering option
        self.opt_cluster = []
        self.opt_append_og = []
        self.opt_align = []
        self.opt_tree = []

        # Running result
        self.rslt_cluster = []
        self.rslt_append_og = []
        self.rslt_align = []
        self.rslt_tree = []

        # multigene_list
        self.multigene_list = None

    def __repr__(self):

        out_dict = {}
        for key in self.dict_dataset:
            out_dict[key] = len(self.dict_dataset[key])
        return f"<< FunID_var object >>\nNumber of FI: {len(self.list_FI)}\nDB genes:{self.list_db_gene}\nQuery genes:{self.list_qr_gene}\nSections:{self.list_sect}\nDataset_dict:{json.dumps(out_dict, indent=2)}"

    def add_dataset(self, sect, gene, list_qr, list_db, list_og):
        data = Dataset(sect, gene, list_qr, list_db, list_og)
        if not (sect) in self.dict_dataset:
            self.dict_dataset[sect] = {}

        if not (gene) in self.dict_dataset[sect]:
            self.dict_dataset[sect][gene] = data
        else:
            logging.warning(f"Overriding dataset on {sect} {gene}")
            self.dict_dataset[sect][gene] = data

        # If gene is concatenated, update ori_species
        if gene == "concatenated":
            for FI in list_qr + list_db:
                FI.bygene_species[gene] = FI.ori_species

    def remove_dataset(self, sect, gene):

        if not (sect) in self.dict_dataset:
            logging.warning(
                f"Passing removing dataset of {sect} {gene} because no {sect} priorly exists"
            )
        elif not (gene) in self.dict_dataset[sect]:
            logging.warning(
                f"Passing removing dataset of {sect} {gene} because no {gene} priorly exists"
            )
        else:
            del self.dict_dataset[sect][gene]
            # if all gene removed
            if len(self.dict_dataset[sect]) == 0:
                del self.dict_dataset[sect]

    def exist_dataset(self, sect, gene):
        if not (sect) in self.dict_dataset:
            return False
        elif not (gene) in self.dict_dataset[sect]:
            return False
        else:
            return True

    # generate dataset by section and gene
    def generate_dataset(self, opt):

        dict_funinfo = {}

        for sect in self.list_sect:
            logging.info(f"Generating dataset for {sect}")
            dict_funinfo[sect] = {}

            # For query only option
            if opt.queryonly is True:

                section_flag = 0  # whether to run this section

                for gene in self.list_db_gene:

                    logging.info(f"Searching {gene}")

                    list_qr = [
                        FI
                        for FI in self.list_FI
                        if (
                            gene in FI.seq
                            and FI.datatype == "query"
                            and FI.adjusted_section == sect
                        )
                    ]

                    # do not manage db when query only mode and query does not exists
                    if len(list_qr) > 0:
                        section_flag = 1

                if section_flag == 1:  # if decided to run this section
                    logging.info(f"Section {sect} passed dataset construction")
                    for gene in self.list_db_gene:
                        list_qr = [
                            FI
                            for FI in self.list_FI
                            if (
                                gene in FI.seq
                                and FI.datatype == "query"
                                and FI.adjusted_section == sect
                            )
                        ]

                        list_db = [
                            FI
                            for FI in self.list_FI
                            if (
                                gene in FI.seq
                                and FI.datatype == "db"
                                and FI.adjusted_section == sect
                            )
                        ]
                        self.add_dataset(sect, gene, list_qr, list_db, [])

                else:
                    logging.warning(
                        f"Section {sect} did not passed dataset construction"
                    )

                # for concatenated
                if opt.concatenate is True:
                    list_qr = [
                        FI
                        for FI in self.list_FI
                        if (FI.datatype == "query" and FI.adjusted_section == sect)
                    ]

                    # do not manage db when query only mode and query does not exists
                    if len(list_qr) > 0:
                        list_db = [
                            FI
                            for FI in self.list_FI
                            if (FI.datatype == "db" and FI.adjusted_section == sect)
                        ]
                        self.add_dataset(sect, "concatenated", list_qr, list_db, [])

            # For DB included opt
            else:
                for gene in self.list_db_gene:
                    list_qr = [
                        FI
                        for FI in self.list_FI
                        if (
                            gene in FI.seq
                            and FI.datatype == "query"
                            and FI.adjusted_section == sect
                        )
                    ]
                    list_db = [
                        FI
                        for FI in self.list_FI
                        if (
                            gene in FI.seq
                            and FI.datatype == "db"
                            and FI.adjusted_section == sect
                        )
                    ]
                    self.add_dataset(sect, gene, list_qr, list_db, [])

                # for concatenated
                if opt.concatenate is True:

                    list_qr = [
                        FI
                        for FI in self.list_FI
                        if (FI.datatype == "query" and FI.adjusted_section == sect)
                    ]
                    list_db = [
                        FI
                        for FI in self.list_FI
                        if (FI.datatype == "db" and FI.adjusted_section == sect)
                    ]
                    self.add_dataset(sect, "concatenated", list_qr, list_db, [])

    # homogenize list_dataset and dicts
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
                            f"Both list_FI and dict_FI have conflicting final species, {FI.final_species} and {self.dict_hash_FI[h].final_species}"
                        )
                        raise Exception

                # adjusted section
                if FI.adjusted_section != self.dict_hash_FI[h].adjusted_section:
                    if FI.adjusted_section == "" or FI.adjusted_section == "":
                        FI.adjusted_section = self.dict_hash_FI[h].adjusted_section
                    elif self.dict_hash_FI[h].adjusted_section == "":
                        self.dict_hash_FI[h].adjusted_section = FI.adjusted_section
                    else:
                        logging.error(
                            f"Both list_FI and dict_FI have conflicting final species, {FI.adjusted_section} and {self.dict_hash_FI[h].adjusted_section}, {FI}"
                        )
                        raise Exception

                elif (
                    FI.adjusted_section == ""
                    and self.dict_hash_FI[h].adjusted_section == ""
                ):
                    if FI.section != "":
                        FI.adjusted_section = FI.section
                        self.dict_hash_FI[h].adjusted_section = FI.section

                    elif self.dict_hash_FI[h].section != "":
                        FI.adjusted_section = self.dict_hash_FI[h].section
                        self.dict_hash_FI[h].adjusted_section = self.dict_hash_FI[
                            h
                        ].section
                    else:
                        logging.error(f"{FI.accession} does not have assigned section!")

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
                            f"Both list_FI and dict_FI have conflicting gene identification results, {FI.bygene_species} and {self.dict_hash_FI[h].bygene}"
                        )
                        raise Exception

    # Remove invalid dataset to be analyzed
    def remove_invalid_dataset(self):

        # collect remove list : removing after iterating
        list_remove = []
        for sect in self.dict_dataset:
            for gene in self.dict_dataset[sect]:
                if len(self.dict_dataset[sect][gene].list_db_FI) == 0:
                    logging.warning(
                        f"Removing {gene} from section {sect} because there are no corresponding sequences"
                    )
                    list_remove.append((sect, gene))
                elif (
                    len(self.dict_dataset[sect][gene].list_db_FI)
                    + len(self.dict_dataset[sect][gene].list_qr_FI)
                    + len(self.dict_dataset[sect][gene].list_og_FI)
                    < 4
                ):
                    logging.warning(
                        f"Removing {sect} from downstream phylogenetic analysis because there are not enough sequences"
                    )
                    list_remove.append((sect, gene))

        for x in list_remove:
            self.remove_dataset(*x)

    # save fasta for outgroup adjusted fasta
    def save_dataset(self, path, opt):
        print(self.dict_dataset)
        for sect in self.dict_dataset:
            for gene in self.dict_dataset[sect]:
                if not (gene == "concatenated"):
                    if (
                        opt.concatenate is True
                        and "concatenated" in self.dict_dataset[sect]
                    ):
                        fasta_list = list(
                            set(
                                self.dict_dataset[sect][gene].list_db_FI
                                + self.dict_dataset[sect][gene].list_qr_FI
                                + self.dict_dataset[sect][gene].list_og_FI
                                + self.dict_dataset[sect]["concatenated"].list_og_FI
                            )
                        )
                    else:
                        fasta_list = (
                            self.dict_dataset[sect][gene].list_db_FI
                            + self.dict_dataset[sect][gene].list_qr_FI
                            + self.dict_dataset[sect][gene].list_og_FI
                        )

                    manage_input.save_fasta(
                        fasta_list,
                        gene,
                        f"{path.out_adjusted}/{opt.runname}_Adjusted_{sect}_{gene}.fasta",
                        by="hash",
                    )

    # Validate if any multiple sequence alignment has no overlapping region
    def validate_alignments(self, path, opt):
        for sect in self.dict_dataset:
            for gene in self.dict_dataset[sect]:

                # Check if dataset exists
                pass

                # Check if alignment corresponding to dataset exists
                if not (
                    os.path.isfile(
                        f"{path.out_alignment}/{opt.runname}_trimmed_{sect}_{gene}.fasta"
                    )
                ):
                    logger.warning(
                        f"Alignment file {path.out_alignment}/{opt.runname}_trimmed_{sect}_{gene}.fasta does not exists"
                    )

                else:
                    # If alignment exists, check if alignment does have overlapping regions

                    ## Parse alignment
                    seq_list = list(
                        SeqIO.parse(
                            f"{path.out_alignment}/{opt.runname}_trimmed_{sect}_{gene}.fasta",
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
                            f"Alignment for {sect} {gene} does not have any overlapping regions! Please check alignment to find out if some of the sequences were from different regions"
                        )
                    else:
                        logging.info(f"Alignment for {sect} {gene} passed validation")

                # If alignment corresponding to dataset does not exists, raise warning or error
                pass

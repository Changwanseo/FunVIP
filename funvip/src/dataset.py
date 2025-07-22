from funvip.src import save
from funvip.src import hasher
from funvip.src import ext
from funvip.src.opt_generator import opt_generator
from Bio import SeqIO
import os
import sys
import shutil
import numpy as np
import multiprocessing as mp
import logging
import re
import json


### Datasets are group of sequences with same genes for analysis, including query, database and outgroup
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

    def __repr__(self):
        return f"<< FunVIP_Dataset object >> Gene: {self.gene} | Group: {self.group} | Query: {len(self.list_qr_FI)} | DB: {len(self.list_db_FI)} | Outgroup: {len(self.list_og_FI)}\n"


### FunVIP all dataset variables bundle
# Called as "V" in main pipeline
class FunVIP_var:
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
        # group dict ->  {group:set(genus1, genus2, ... ) format}
        # used in synchronizing to check if the genus is main interest of the group
        self.dict_group = {}

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
        self.multigene_list = []

        # Partition information by dataset
        # group : {"len" : dict, "order" : list} format
        self.partition = {}

    def __repr__(self):
        out_dict = {}
        for key in self.dict_dataset:
            out_dict[key] = len(self.dict_dataset[key])
        return f"<< FunVIP_var object >>\nNumber of FI: {len(self.list_FI)}\nDB genes:{self.list_db_gene}\nQuery genes:{self.list_qr_gene}\ngroups:{self.list_group}\nDataset_dict:{json.dumps(out_dict, indent=2)}"

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

        # Update dict_group
        if not group in self.dict_group:
            self.dict_group[group] = set()

        for FI in list_db:
            self.dict_group[group].add(FI.genus)

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
            # 1 for concatenated
            if len(self.dict_dataset[group]) <= 1:
                del self.dict_dataset[group]
                logging.info(f"Removed {group} from dataset")
            else:
                logging.info(f"Removed {group} {gene} from dataset")

    # Check if given dict_dataset[group][gene] exists
    def exist_dataset(self, group, gene):
        try:
            self.dict_dataset[group][gene]
            return True
        except:
            return False

    # Check if dict_group has been properly generated
    def check_dict_group(self, opt):
        # If level is lower than genus, check if each group includes only one genus
        if opt.level in ["subseries", "series", "subsection", "section", "genus"]:
            for group in self.dict_group:
                if len(self.dict_group[group]) > 1:
                    logging.warning(
                        f"{opt.level} {group} includes more than one genus, {self.dict_group[group]}.\n",
                        f"This may cause unusual behaviour in synchronizing",
                    )

        # If level is higher than genus, check if each genus belongs to more than one group
        elif opt.level in [
            "subtribe",
            "tribe",
            "subfamily",
            "family",
            "suborder",
            "order",
            "subclass",
            "class",
            "subdivision",
            "division",
            "subphylum",
            "phylum",
            "subkingdom",
            "kingdom",
        ]:
            reverse_dict = {}
            for group in self.dict_group:
                for genus in self.dict_group[group]:
                    reverse_dict[genus] = set()
                reverse_dict[genus].add(group)

            for genus in reverse_dict:
                if len(reverse_dict[genus]) > 1:
                    logging.warning(
                        f"genus {group} belongs to more than one {opt.level}, {reverse_dict[genus]}.\n",
                        f"This may cause unusual behaviour in synchronizing",
                    )

        else:
            logging.error(f"DEVELOPMENTAL ERROR, UNEXPECTED LEVEL {opt.level} selected")
            raise Exception

    # generate dataset by group and gene
    def generate_dataset(self, opt):
        # Format : dict_funinfo = {group: {gene : [FI]}}
        dict_funinfo = {}

        for group in self.list_group:
            logging.info(f"Generating dataset for {group}")

            # print(f"opt.queryonly: {opt.queryonly}")

            dict_funinfo[group] = {}

            # For queryonly case
            if opt.queryonly is True:
                # whether to run this group
                group_flag = False
                for gene in self.list_db_gene:
                    logging.debug(
                        f"Searching dataset {group} {gene} includes query sequences"
                    )
                    list_qr = [
                        FI
                        for FI in self.list_FI
                        if (
                            gene in FI.seq
                            and FI.datatype == "query"
                            and FI.adjusted_group == group
                        )
                    ]

                    # do not manage db when --queryonly True (--all False) and query does not exists
                    if len(list_qr) > 0:
                        group_flag = True

                # if decided to run this group
                if group_flag is True:
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

            # For opt.queryonly is False -> run all dataset in database if possible
            else:
                for gene in self.list_db_gene:
                    list_qr = []
                    for FI in self.list_FI:
                        if (
                            FI.datatype == "query"
                            and FI.adjusted_group == group
                            and gene in FI.seq
                        ):
                            if FI.seq[gene] != "":
                                list_db.append(FI)

                    list_db = []
                    for FI in self.list_FI:
                        if (
                            FI.datatype == "db"
                            and FI.adjusted_group == group
                            and gene in FI.seq
                        ):
                            if FI.seq[gene] != "":
                                list_db.append(FI)

                    # If none of the database is possible for this group and gene pair, it should be excluded
                    if len(list_db) > 0:
                        self.add_dataset(group, gene, list_qr, list_db, [])
                    else:
                        logging.warning(
                            f"dataset {group} {gene} did not passed dataset construction due to lack of sequences"
                        )

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

        self.check_dict_group(opt)

    # homogenize list_dataset and dict_hash_FI from multiple results
    def homogenize_dataset(self):
        for FI in self.list_FI:
            if FI.hash in self.dict_hash_FI:
                h = FI.hash
                # final species
                if FI.final_species != self.dict_hash_FI[h].final_species:
                    # If final species was empty
                    if FI.final_species == "":
                        FI.final_species = self.dict_hash_FI[h].final_species
                    # If final species in hash dict was empty
                    elif self.dict_hash_FI[h].final_species == "":
                        self.dict_hash_FI[h].final_species = FI.final_species
                    # If they collides, it is error
                    else:
                        logging.error(
                            f"DEVELOPMNETAL ERROR Both list_FI and dict_hash_FI have conflicting final species, {FI.final_species} and {self.dict_hash_FI[h].final_species} for hash {h}"
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
                        logging.warning(f"{FI.id} does not have assigned group!")

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
        # Before save_dataset, remove pre-existing results
        for file in os.listdir(path.out_adjusted):
            if file.endswith(".fasta"):
                os.remove(f"{path.out_adjusted}/{file}")

        for file in os.listdir(f"{path.out_adjusted}/hash/"):
            if file.endswith(".fasta"):
                os.remove(f"{path.out_adjusted}/hash/{file}")

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

                    # Remove no seqs
                    for fasta in fasta_list:
                        save.save_fasta(
                            fasta_list,
                            gene,
                            f"{path.out_adjusted}/{opt.runname}_Adjusted_{group}_{gene}.fasta",
                            by="hash",
                        )

    # Validate if any multiple sequence alignment has no overlapping region
    def validate_alignments(self, path, opt):
        # 1. Manual validation for illegal alignments
        fail_list = []
        remove_dict = {}
        tree_hash_dict = hasher.encode(self.list_FI, newick=True)
        for group in self.dict_dataset:
            remove_dict[group] = {}
            for gene in self.dict_dataset[group]:
                if gene != "concatenated":
                    remove_dict[group][gene] = []
                    # Check if alignment corresponding to dataset exists
                    if not (
                        os.path.isfile(
                            f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta"
                        )
                    ):
                        logger.warning(
                            f"Alignment file {path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta does not exists"
                        )

                        fail_list.append(group_gene)

                    else:
                        # If alignment exists, check if alignment does have overlapping regions
                        ## Parse alignment
                        seq_list = list(
                            SeqIO.parse(
                                f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta",
                                "fasta",
                            )
                        )

                        # Remove sequences that has not been existed during alignment stage
                        seq_id_list = [seq.id for seq in seq_list]
                        for FI in self.dict_dataset[group][gene].list_db_FI:
                            if not (FI.hash in seq_id_list):
                                self.dict_dataset[group][gene].list_db_FI.remove(FI)

                        for FI in self.dict_dataset[group][gene].list_qr_FI:
                            if not (FI.hash in seq_id_list):
                                self.dict_dataset[group][gene].list_qr_FI.remove(FI)

                        for FI in self.dict_dataset[group][gene].list_og_FI:
                            if not (FI.hash in seq_id_list):
                                self.dict_dataset[group][gene].list_og_FI.remove(FI)

                        # Remove empty sequences
                        remove_hash = []
                        for seq in seq_list:
                            if len(str(seq.seq).replace("-", "")) == 0:
                                logging.debug(
                                    f"{group} {gene} {seq.id} : {len(str(seq.seq).replace('-', ''))}"
                                )
                                remove_hash.append(seq.id)

                        remove_dict[group][gene] = remove_hash

                        # Remove unusable sequence and re-read it
                        ## db_list, query_list, outgroup_list might has to be changed
                        seq_list = [
                            seq for seq in seq_list if not seq.id in remove_hash
                        ]
                        SeqIO.write(
                            seq_list,
                            f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta",
                            "fasta",
                        )

                        ## Transform alignment to vector
                        ## Change gap to 0, and other characters to 1
                        vectors = [
                            np.fromiter(
                                re.sub(
                                    r"[^0]", "1", re.sub(r"[\-]", "0", str(seq.seq))
                                ),
                                dtype=np.int32,
                            )
                            for seq in seq_list
                        ]
                        ## Multiply vectors
                        vector_products = np.prod(np.vstack(vectors), axis=0)
                        logging.debug(
                            f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta resulted multiplied vector {vector_products}"
                        )

                        ## If all value of vectors are zero, it means that all regions have at least one gap
                        ## Raise warning for this
                        if np.all((vector_products == 0)):
                            logging.critical(
                                f"Alignment for {group} {gene} does not have any overlapping regions! Removing from analysis"
                            )

                            # Report one non_overlapping pair
                            for i in range(len(seq_list) - 1):
                                seq1 = seq_list[i]
                                seq1_hash = seq1.id
                                seq1_str = str(seq.seq)

                                for j in range(i + 1, len(seq_list)):
                                    seq2 = seq_list[j]
                                    seq2_hash = seq2.id
                                    seq2_str = str(seq.seq)

                                    has_overlap = False

                                    for k in range(len(seq1_str)):
                                        if seq1_str[k] != "-" and seq2_str[k] != "-":
                                            has_overlap = True
                                            break

                            seq1_id = self.dict_hash_FI[seq1_hash]
                            seq2_id = self.dict_hash_FI[seq2_hash]

                            logging.critical(
                                f"At least one pair of sequence does not overlap, such as {seq1_id} and {seq2_id}"
                            )

                            fail_list.append((group, gene))
                            # for tree, use hash dict with genus and species information
                            # Decoding process in done in tree building processes, so these alignments cannot be decoded. So decode them here
                            shutil.move(
                                f"{path.out_alignment}/{opt.runname}_MAFFT_{group}_{gene}.fasta",
                                f"{path.out_alignment}/hash/{opt.runname}_hash_MAFFT_{group}_{gene}.fasta",
                            )

                            hasher.decode(
                                tree_hash_dict,
                                f"{path.out_alignment}/hash/{opt.runname}_hash_MAFFT_{group}_{gene}.fasta",
                                f"{path.out_alignment}/{opt.runname}_MAFFT_{group}_{gene}.fasta",
                            )

                            # Move failed alignment to failed directory
                            shutil.move(
                                f"{path.out_alignment}/{opt.runname}_MAFFT_{group}_{gene}.fasta",
                                f"{path.out_alignment}/failed/{opt.runname}_MAFFT_{group}_{gene}.fasta",
                            )
                        else:
                            logging.debug(
                                f"Alignment for {group} {gene} passed validation"
                            )

                    # If alignment corresponding to dataset does not exists, raise warning or error
                    pass

        # Remove bad datasets
        for fail in fail_list:
            group = fail[0]
            gene = fail[1]
            self.dict_dataset[group].pop(gene)
            logging.warning(f"Alignment for {group} {gene} has removed from analysis")

        # Add issue
        for fail in fail_list:
            for FI in self.list_FI:
                if (
                    FI.adjusted_group == fail[0]
                    and fail[1] in FI.seq
                    and FI.seq[fail[1]] != ""
                ):
                    FI.issues.add(f"alignfail:{fail[1]}")

        logging.debug("Remove dict")
        logging.debug(remove_dict)

        # Remove removed sequences from dataset
        for group in remove_dict:
            for gene in remove_dict[group]:
                if group in self.dict_dataset:
                    if gene in self.dict_dataset[group]:
                        for _hash in remove_dict[group][gene]:
                            if _hash in self.dict_dataset[group][gene].list_qr_FI:
                                self.dict_dataset[group][gene].list_qr_FI.remove(
                                    self.dict_hash_FI[_hash]
                                )
                                logging.warning(
                                    f"{self.dict_hash_ID[_hash]} removed from dataset {group} {gene}. Please check the alignment and see the region is correct"
                                )
                            if _hash in self.dict_dataset[group][gene].list_db_FI:
                                self.dict_dataset[group][gene].list_db_FI.remove(
                                    self.dict_hash_FI[_hash]
                                )
                                logging.warning(
                                    f"{self.dict_hash_ID[_hash]} removed from dataset {group} {gene}. Please check the alignment and see the region is correct"
                                )
                            if _hash in self.dict_dataset[group][gene].list_og_FI:
                                self.dict_dataset[group][gene].list_og_FI.remove(
                                    self.dict_hash_FI[_hash]
                                )
                                logging.warning(
                                    f"{self.dict_hash_ID[_hash]} removed from dataset {group} {gene}. Please check the alignment and see the region is correct"
                                )

        # Finally, check again if the datasets meet criteria
        final_fail_list = []
        for group in self.dict_dataset:
            for gene in self.dict_dataset[group]:
                logging.debug(f"Validating alignment for {group} {gene}")
                logging.debug(
                    f"list_qr_FI : {len(self.dict_dataset[group][gene].list_qr_FI)}"
                )
                logging.debug(
                    f"list_db_FI : {len(self.dict_dataset[group][gene].list_db_FI)}"
                )
                logging.debug(
                    f"list_og_FI : {len(self.dict_dataset[group][gene].list_og_FI)}"
                )

                if (
                    len(self.dict_dataset[group][gene].list_qr_FI)
                    + len(self.dict_dataset[group][gene].list_db_FI)
                    + len(self.dict_dataset[group][gene].list_og_FI)
                    < 4
                ):
                    logging.warning(
                        f"After removing invalid datasets from dataset {group} {gene}, the number of remaining sequences are under 4, removing from anlaysis."
                    )
                    final_fail_list.append((group, gene))

                elif len(self.dict_dataset[group][gene].list_og_FI) < 1:
                    logging.warning(
                        f"After removing invalid datasets from dataset {group} {gene}, outgroup of the dataset has completely removed, removing from anlaysis."
                    )
                    final_fail_list.append((group, gene))

                elif (
                    len(self.dict_dataset[group][gene].list_qr_FI)
                    + len(self.dict_dataset[group][gene].list_db_FI)
                    < 1
                ):
                    logging.warning(
                        f"After removing invalid datasets from dataset {group} {gene}, dataset has completely removed, removing from anlaysis."
                    )
                    final_fail_list.append((group, gene))

        # Remove bad datasets
        for fail in final_fail_list:
            group = fail[0]
            gene = fail[1]
            self.dict_dataset[group].pop(gene)

        # If concatenated is the only left dataset, remove entire group
        group_pop_list = []
        for group in self.dict_dataset:
            if len(self.dict_dataset[group].keys()) == 0:
                group_pop_list.append(group)
            elif (
                len(self.dict_dataset[group].keys()) == 1
                and "concatenated" in self.dict_dataset[group].keys()
            ):
                group_pop_list.append(group)

        for group in group_pop_list:
            self.dict_dataset.pop(group)
            logging.warning(
                f"No alignment left for group {group}, removed from analysis"
            )

        # Add issue: the number of sequences are insufficient
        for fail in final_fail_list:
            for FI in self.list_FI:
                if FI.adjusted_group == fail[0]:
                    if FI.seq[fail[1]] != "":
                        FI.issues.add(f"lackseq")

        # Validate multiple sequence alignment with TCS score from T-COFFEE
        # As T-COFFEE build is only available in Mac and Linux, should check if it is available
        if opt.method.tcs is True:
            bad_cnt = 0

            if opt.verbose < 3:
                tcs_opt = opt_generator(
                    V=self, opt=opt, path=path, step="tcs", thread=1
                )
                p = mp.Pool(opt.thread)
                p.starmap(ext.TCS, tcs_opt)
                p.close()
                p.join()

            else:
                tcs_opt = opt_generator(V=self, opt=opt, path=path, step="tcs")
                for option in tcs_opt:
                    ext.TCS(*option)

                # non-multithreading mode for debugging
                for group in self.dict_dataset:
                    for gene in self.dict_dataset[group]:
                        # Running TCS for concatenated alignment is duplicate
                        if gene != "concatenated":
                            tcs_out = f"{path.out_alignment}/alignment/{opt.runname}_{group}_{gene}.tcs"
                            # Parse tcs result
                            with open(tcs_out, "r") as f_tcs:
                                tcs_result_raw = f_tcs.read()
                                tcs_result = tcs_result_raw.split("*")[2].split("cons")[
                                    0
                                ]
                                for line in tcs_result.split("\n")[1:-1]:
                                    _hash = line.split(":")[0].strip()
                                    tcs_score = int(line.split(":")[1].strip())
                                    if (
                                        tcs_score < 50
                                    ):  # cutoff 50 comes from TCS documentation
                                        FI_id = self.dict_hash_FI[_hash].id
                                        logging.warning(
                                            f"{FI_id} has poor alignment score in {group} {gene}"
                                        )
                                        bad_cnt += 1

            if bad_cnt == 0:
                logging.info(f"All sequences in alignment passed TCS validation")
            else:
                logging.warning(
                    f"{bad_cnt} sequences in alignment failed TCS validation. Please check sequneces"
                )

    # check inconsistency exists along identification result of each genes
    def check_inconsistent(self):
        set_gene = set(self.list_db_gene + self.list_qr_gene)

        for _hash in self.dict_hash_FI:
            FI = self.dict_hash_FI[_hash]

            # Collect result from only used sequences
            if FI.adjusted_group in self.dict_dataset:
                # If any of appropriate gene used
                if any(
                    key in self.dict_dataset[FI.adjusted_group] for key in FI.seq.keys()
                ):
                    inconsistent_flag = 0

                    for gene in set_gene:
                        # Check if data analysis had performed for specific FI, group, gene combination
                        if (
                            gene in FI.bygene_species
                            and len(FI.seq[gene]) > 0
                            and gene in self.dict_dataset[FI.adjusted_group]
                        ):
                            # Check inconsistent identification across genes
                            if not (
                                any(
                                    _sp in FI.final_species
                                    for _sp in FI.bygene_species[gene].split("/")
                                )
                            ):
                                inconsistent_flag = 1

                    if inconsistent_flag == 1:
                        FI.issues.add("inconsistent")

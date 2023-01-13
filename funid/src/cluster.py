import multiprocessing as mp
import pandas as pd
import numpy as np
import shutil
import os, sys, subprocess
import copy
from functools import lru_cache
from time import sleep

from Bio.Blast import NCBIXML
from Bio import SeqIO

import logging
from funid.src.visualize import plot_pca, plot_heatmap
from funid.src.ext import mmseqs


def check_gene_availability(V, opt):

    available_gene = []
    unavailable_gene = []

    if opt.queryonly is True:
        available_gene = V.list_qr_gene
    else:
        for gene in opt.gene:
            if gene in V.dict_gene_SR:
                available_gene.append(gene)
            else:
                unavailable_gene.append(gene)

    return available_gene, unavailable_gene


def get_naive_section(V):
    return list(
        set([funinfo.section for funinfo in V.list_FI if type(funinfo.section) == str])
    )


def append_concatenated_query_section(V):

    section_dict = {}
    for funinfo in V.list_FI:
        section_dict[funinfo.hash] = funinfo.adjusted_section

    V.cSR["query_section"] = V.cSR["qseqid"].apply(lambda x: section_dict.get(x))

    return V


def append_query_section(V):

    section_dict = {}
    for FI in V.list_FI:
        section_dict[FI.hash] = FI.adjusted_section

    for gene in V.dict_gene_SR:
        df = V.dict_gene_SR[gene]
        df["query_section"] = df["qseqid"].apply(lambda x: section_dict.get(x))

    return V


# assign gene to unclassified gene by search result
def assign_gene(result_dict, V, cutoff=0.99):

    # if no query exists to assign, return dataset without operations
    if len(result_dict.keys()) == 0:
        return V

    for result in result_dict:
        # add gene column to result_dict
        result_dict[result]["gene"] = result

    print(result_dict)

    for result in result_dict:
        print(result_dict[result])

    # combine all result in single dataframe
    gene_result_all = pd.concat([result_dict[result] for result in result_dict], axis=0)

    # split by query
    gene_result_grouped = gene_result_all.groupby(gene_result_all["qseqid"])
    # print(sorted(list(set(gene_result_all["qseqid"]))))

    for funinfo in V.list_FI:

        # for seq of funinfo
        for n, seq in enumerate(funinfo.unclassified_seq):

            try:
                # get each of the dataframe for each funinfo
                current_df = gene_result_grouped.get_group(f"{funinfo.hash}_{n}")
                # sort by bitscore
                # sorting has peformed after split for better performance
                current_df.sort_values(by=["bitscore"], inplace=True, ascending=False)
                # reset index to easily get maximum
                current_df.reset_index(inplace=True, drop=True)
                # get result stasifies over cutoff
                cutoff_df = current_df[
                    current_df["bitscore"] > current_df["bitscore"][0] * cutoff
                ]
                gene_count = len(set(cutoff_df["gene"]))
                gene_list = list(set(cutoff_df["gene"]))

            except:
                gene_count = 0
                gene_list = []

            # if only 1 gene available, take it
            if gene_count == 1:
                funinfo.update_seq(gene_list[0], seq)
            # if no gene matched, warn it
            elif gene_count == 0:
                logging.warning(
                    f"[Warning] Query seq in {funinfo.accession} cannot be assigned to gene. Check sequence. Skipping {funinfo.accession}"
                )
            elif gene_count >= 2:
                logging.warning(
                    f"[Warning] Query seq in {funinfo.accession} has multiple matches to gene."
                )
                funinfo.update_seq(gene_list[0], seq)

            else:
                logging.error("[ERROR] DEVELOPMENTAL ERROR IN GENE ASSIGN")
                raise Exception
    return V


def cluster(FI, df_search, V, path, opt):

    list_sect = copy.deepcopy(V.list_sect)
    # delete V to reduce memory assumption
    del V

    # If no evidence available, return it
    if df_search is None:
        FI.adjusted_section = FI.section
        return FI, None

    # Clustering
    # for sequence with section, retain it
    elif not (FI.section == ""):  # or type(FI.section) != str):
        FI.adjusted_section = FI.section
        return FI, FI.adjusted_section

    # update section if sequence does not have section
    else:
        # sorting has peformed after split for better performance
        df_search.sort_values(by=["bitscore"], inplace=True, ascending=False)
        # reset index to easily get maximum
        df_search.reset_index(inplace=True, drop=True)
        # get result stasifies over cutoff
        cutoff_df = df_search[
            df_search["bitscore"] > df_search["bitscore"][0] * opt.cluster.cutoff
        ]

        section_count = len(set(cutoff_df["subject_section"]))
        list_sect = list(set(cutoff_df["subject_section"]))

        # if first time of section update
        if FI.adjusted_section == "" or FI.adjusted_section == "":
            # if only 1 section available, take it
            if section_count == 1:
                FI.adjusted_section = list_sect[0]
            # if no section matched, warn it
            elif section_count == 0:
                logging.warning(
                    f"[Warning] Query seq in {FI.accession} cannot be assigned to section. Check sequence"
                )
            elif section_count >= 2:
                logging.warning(
                    f"[Warning] Query seq in {FI.accession} has multiple matches to section, {list_sect}."
                )
                FI.adjusted_section = list_sect[0]
            else:
                logging.error("[ERROR] DEVELOPMENTAL ERROR IN SECTION ASSIGN")
                raise Exception

            logging.info(
                f"{FI.accession} {FI.description} has clustered to {FI.adjusted_section}"
            )

        # if section already updated
        else:
            if not (FI.adjusted_section == ""):
                if not (FI.adjusted_section in list_sect):
                    logging.warning(
                        f"[Warning] Clustering result colliding in {FI.accession}"
                    )

        # return funinfo object with adjusted section, and selected section
        return FI, list_sect[0]


# appending outgroup
def append_outgroup(V, df_search, gene, section, path, opt):

    logging.info(f"Appending outgroup on Section:{section}, Gene:{gene}")

    # funinfo for designated sections
    list_FI_return = [
        funinfo for funinfo in V.list_FI if funinfo.adjusted_section == section
    ]

    list_FI = copy.deepcopy(V.list_FI)

    del V

    # ready for by sseqid hash, which section to append
    # this time, append adjusted section
    section_dict = {}
    funinfo_dict = {}

    for funinfo in list_FI:
        section_dict[funinfo.hash] = funinfo.adjusted_section
        funinfo_dict[funinfo.hash] = funinfo

    if gene != "concatenated":
        # generate minimal bitscore cutoff that does not overlaps to query inside value

        cutoff_set_df = df_search[df_search["subject_section"] == section]

        try:
            bitscore_cutoff = min(cutoff_set_df["bitscore"])
        except:
            bitscore_cutoff = 999999  # Infinite

        # get result stasifies over cutoff
        cutoff_df = df_search[df_search["bitscore"] < bitscore_cutoff]

        cutoff_df = cutoff_df[cutoff_df["bitscore"] > 0]
        # split that same section to include all to alignment, and left other sections for outgroup selection
        cutoff_df = cutoff_df[cutoff_df["subject_section"] != section]

        if cutoff_df.groupby(["subject_section"]).count().empty:
            logging.warning(
                f"[Warning] Not enough outgroup sequences matching for {section} {gene}. There might be outlier sequence that does not matches to section. Trying flexible cutoff"
            )
            cutoff_df = df_search[df_search["bitscore"] > 0]
            cutoff_df = cutoff_df[cutoff_df["subject_section"] != section]

        elif (
            cutoff_df.groupby(["subject_section"]).count()["sseqid"].max()
            < opt.maxoutgroup
        ):
            logging.warning(
                f"[Warning] Not enough outgroup sequences matching for {section} {gene}. There might be outlier sequence that does not matches to section. Trying flexible cutoff"
            )
            cutoff_df = df_search[df_search["bitscore"] > 0]
            cutoff_df = cutoff_df[cutoff_df["subject_section"] != section]
        # else:
        #    Mes(f"Passed outgroup counting")

    # in concatenated dataset, upper method makes problem when gene dataset were biased. using old method instead
    # Work from here

    else:
        cutoff_df = df_search[df_search["bitscore"] > 0]
        cutoff_df = cutoff_df[cutoff_df["subject_section"] != section]

    # sort by bitscore
    cutoff_df.sort_values(by=["bitscore"], inplace=True, ascending=False)
    # reset index to easily get maximum
    cutoff_df.reset_index(inplace=True, drop=True)

    # iterate until designated number of sequences from the most closest section selected
    # we sould get outgroup from columns
    outgroup_dict = {}
    max_cnt = 0
    max_section = ""

    for n, subject_section in enumerate(cutoff_df["subject_section"]):
        # if first sequence of the section found, make new key to dict
        if not (subject_section) in outgroup_dict:
            # print(funinfo_dict)
            outgroup_dict[subject_section] = [funinfo_dict[cutoff_df["sseqid"][n]]]
        # if section record already exists, append it
        else:
            if (
                not (funinfo_dict[cutoff_df["sseqid"][n]])
                in outgroup_dict[subject_section]
            ):
                outgroup_dict[subject_section].append(
                    funinfo_dict[cutoff_df["sseqid"][n]]
                )

        # if enough outgroup sequences found
        if len(outgroup_dict[subject_section]) >= opt.maxoutgroup:
            logging.info(f"Outgroup {subject_section} were selected to {section}")
            logging.info(
                f"Outgroup selection for section {section} : {outgroup_dict[subject_section]}"
            )

            # include better blast match to intrasectional sequences to analysis
            ## should be added as option in further development!
            """
            if opt.include_ambiguous_seqs is True:
                for sect in outgroup_dict:
                    if not sect == subject_section:
                        list_FI_return += outgroup_dict[sect]
            """

            return (
                gene,
                section,
                outgroup_dict[subject_section],
                outgroup_dict[subject_section] + list_FI_return,
            )
        else:
            if len(outgroup_dict[subject_section]) > max_cnt:
                max_cnt = len(outgroup_dict[subject_section])
                max_section = subject_section

    logging.warning(
        f"[Warning] Not enough sequences are available for outgroup number {opt.maxoutgroup} in {section}, using {max_section} despite of lower number"
    )

    if not (max_section) == "":
        logging.info(
            f"Final outgroup selection for section {section} : {outgroup_dict[max_section]}"
        )

        # include better blast match to intrasectional sequences to analysis
        ## should be added as option in further development!
        """
        if opt.include_ambiguous_seqs is True:
            for sect in outgroup_dict:
                if not sect == subject_section:
                    list_FI_return += outgroup_dict[sect]
        """

        return (
            gene,
            section,
            outgroup_dict[max_section],
            outgroup_dict[max_section] + list_FI_return,
        )
    else:
        logging.warning(f"[Warning] No outgroup available for {section}")
        return (gene, section, [], list_FI_return)


def section_cluster_opt_generator(V, opt, path):

    if len(V.list_qr_gene) == 0:
        logging.error(
            "[ERROR] In section cluster option generator, no possible gene were selected"
        )
        raise Exception

    # For non concatenated analysis or if only one gene exists
    elif len(V.list_qr_gene) <= 1 or V.cSR is None:

        # For caching
        df_group_dict = {}
        qseqid_dict = {}
        for gene in V.list_qr_gene:
            qseqid_dict[gene] = list(set(V.dict_gene_SR[gene]["qseqid"]))
            df_group_dict[gene] = V.dict_gene_SR[gene].groupby(
                V.dict_gene_SR[gene]["qseqid"]
            )

        if opt.queryonly is True:
            list_FI = [FI for FI in V.list_FI if FI.datatype == "query"]
        else:
            list_FI = V.list_FI

        for FI in list_FI:
            if len(list(FI.seq.keys())) == 1:
                gene = list(FI.seq.keys())[0]
                appropriate_df = df_group_dict[gene].get_group(FI.hash)

            elif len(list(FI.seq.keys())) == 0:
                appropriate_df = None

            else:
                logging.warning(
                    f"[Warning] {FI} has matches to multiple genes, but concatenation option not selected"
                )

                appropriate_df = None

                for gene in list(FI.seq.keys()):
                    if appropriate_df is None:
                        appropriate_df = df_group_dict[gene].get_group(FI.hash)
                    else:
                        if max(appropriate_df["bitscore"]) < max(
                            df_group_dict[gene].get_group(FI.hash)["bitscore"]
                        ):
                            appropriate_df = df_group_dict[gene].get_group(FI.hash)

            V.opt_cluster.append((FI, appropriate_df, V, path, opt))

    # For concatenated analysis
    else:
        # cluster section by concatenated search result
        df_group = V.cSR.groupby(V.cSR["qseqid"])
        list_accession = list(set(V.cSR["qseqid"]))
        for FI in V.list_FI:
            if FI.hash in list_accession:
                df_search = df_group.get_group(FI.hash)
            else:
                df_search = None
            V.opt_cluster.append((FI, df_search, V, path, opt))

    return V


# opts ready for multithreading in outgroup append
def outgroup_append_opt_generator(V, path, opt):

    # if concatenated analysis is false
    # Assign different outgroup for each dataset
    for gene in opt.gene:
        for sect in V.list_sect:
            if V.exist_dataset(sect, gene) is True:
                try:
                    df = V.dict_gene_SR[gene]
                    df_group = df.groupby(df["query_section"])
                    df_sect = df_group.get_group(sect)
                    # Generating outgroup opt for multiprocessing
                    V.opt_append_og.append((V, df_sect, gene, sect, path, opt))
                except:
                    logging.warning(
                        f"[Warning] {sect} / {gene} dataset exists, but cannot append outgroup due to no corresponding search result"
                    )
                    pass

    # if concatenated analysis is true
    # concatenated
    if opt.concatenate is True:
        for sect in V.list_sect:
            if V.exist_dataset(sect, "concatenated") is True:
                try:
                    df = V.cSR
                    df_group = df.groupby(df["query_section"])
                    df_sect = df_group.get_group(sect)
                    # Generating outgroup opt for multiprocessing
                    V.opt_append_og.append(
                        (V, df_sect, "concatenated", sect, path, opt)
                    )
                except:
                    logging.warning(
                        f"[Warning] {sect} / concatenated dataset exists, but cannot append outgroup due to no corresponding search result"
                    )
                    pass

    return V


# multiprocessing run result collector for outgrouping
def outgroup_result_collector(V):
    for result in V.rslt_append_og:
        gene = result[0]
        section = result[1]
        result_outgroup_list = result[2]
        result_list_sect = result[3]

        # append outgroup lists
        V.dict_dataset[section][gene].list_og_FI = result_outgroup_list

    return V


## Main cluster pipe
def pipe_cluster(V, opt, path):
    # If clustering enabled
    if opt.method.search in ("blast", "mmseqs"):

        logging.info("[INFO] Sectional clustering")

        # cluster opt generation for multiprocessing
        V = section_cluster_opt_generator(V, opt, path)

        # run multiprocessing start
        if opt.verbose < 3:
            p = mp.Pool(opt.thread)
            V.rslt_cluster = p.starmap(cluster, V.opt_cluster)
            p.close()
            p.join()

        else:
            # non-multithreading mode for debugging
            V.rslt_cluster = [cluster(*o) for o in V.opt_cluster]

        # gather cluster result
        for r in V.rslt_cluster:
            FI = r[0]

        # replace section assigning result
        replace_FI = [r[0] for r in V.rslt_cluster]
        replace_hash_FI = [FI.hash for FI in replace_FI]
        V.list_FI = [
            FI for FI in V.list_FI if not (FI.hash in replace_hash_FI)
        ] + replace_FI

        V.list_sect = list(set([r[1] for r in V.rslt_cluster if r[1] is not None]))

        if opt.queryonly is True:
            for FI in V.list_FI:
                if FI.datatype == "db":
                    FI.adjusted_section = FI.section

    # If not, try to use original sections in tabled format
    else:
        logging.info(
            "[INFO] No searching method designated. Trying to use designated section"
        )
        section_list = get_naive_section(V)
        ## [WIP] Need to make validation process and warn if the input does not have
        for FI in V.list_FI:
            FI.adjusted_section = FI.section

    return V, opt, path


## Main outgroup appending pipeline
def pipe_append_outgroup(V, path, opt):

    V = outgroup_append_opt_generator(V, path, opt)

    # run multiprocessing start
    if opt.verbose < 3:
        p = mp.Pool(opt.thread)
        V.rslt_append_og = p.starmap(append_outgroup, V.opt_append_og)
        p.close()
        p.join()

    else:
        # non-multithreading mode for debugging
        V.rslt_append_og = [append_outgroup(*o) for o in V.opt_append_og]

    return V, path, opt

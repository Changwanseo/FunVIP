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


def get_naive_group(V):
    return list(
        set([funinfo.group for funinfo in V.list_FI if type(funinfo.group) == str])
    )


# Append query_group information column to dict_gene_SR
def append_query_group(V):

    # For normal gene matrix
    group_dict = {}
    for FI in V.list_FI:
        group_dict[FI.hash] = FI.adjusted_group

    for gene in V.dict_gene_SR:
        V.dict_gene_SR[gene]["query_group"] = V.dict_gene_SR[gene]["qseqid"].apply(
            lambda x: group_dict.get(x)
        )

    # For concatenated gene matrix
    group_dict = {}
    for funinfo in V.list_FI:
        group_dict[funinfo.hash] = funinfo.adjusted_group

    V.cSR["query_group"] = V.cSR["qseqid"].apply(lambda x: group_dict.get(x))

    return V


# assign gene to unclassified gene by search result
def assign_gene(result_dict, V, cutoff=0.99):

    # if no query exists to assign, return dataset without operations
    if len(result_dict.keys()) == 0:
        return V

    for result in result_dict:
        # add gene column to result_dict
        result_dict[result]["gene"] = result

    # combine all result in single dataframe
    gene_result_all = pd.concat([result_dict[result] for result in result_dict], axis=0)

    # split by query
    gene_result_grouped = gene_result_all.groupby(gene_result_all["qseqid"])

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
                logging.info(
                    f" Query seq in {funinfo.id} has assigned to {gene_list[0]}."
                )
            # if no gene matched, warn it
            elif gene_count == 0:
                logging.warning(
                    f" Query seq in {funinfo.id} cannot be assigned to any gene. Check sequence. Skipping {funinfo.id}"
                )
            elif gene_count >= 2:
                logging.warning(
                    f" Query seq in {funinfo.id} has multiple matches to gene."
                )
                funinfo.update_seq(gene_list[0], seq)

            else:
                logging.error("DEVELOPMENTAL ERROR IN GENE ASSIGN")
                raise Exception
    return V


# cluster each of Funinfo object and assign group
def cluster(FI, df_search, V, path, opt):

    list_group = copy.deepcopy(V.list_group)

    # delete V to reduce memory assumption
    del V

    # If no evidence available, return it
    if df_search is None:
        FI.adjusted_group = FI.group
        return FI, None
    # for db sequence with group, retain it
    elif not (FI.group == "") and FI.datatype == "db":  # or type(FI.group) != str):
        FI.adjusted_group = FI.group
        return FI, FI.adjusted_group

    # update group if sequence does not have group
    else:
        # sorting has peformed after split for better performance
        df_search.sort_values(by=["bitscore"], inplace=True, ascending=False)
        # reset index to easily get maximum
        df_search.reset_index(inplace=True, drop=True)
        # get result stasifies over cutoff
        cutoff_df = df_search[
            df_search["bitscore"] > df_search["bitscore"][0] * opt.cluster.cutoff
        ]

        group_count = len(set(cutoff_df["subject_group"]))
        list_group = list(set(cutoff_df["subject_group"]))

        # if first time of group update
        if FI.adjusted_group == "" or FI.adjusted_group == "":
            # if only 1 group available, take it
            if group_count == 1:
                FI.adjusted_group = list_group[0]
            # if no group matched, warn it
            elif group_count == 0:
                logging.warning(
                    f" Query seq in {FI.id} cannot be assigned to group. Check sequence"
                )
            elif group_count >= 2:
                logging.warning(
                    f" Query seq in {FI.id} has multiple matches to group, {list_group}."
                )
                FI.adjusted_group = list_group[0]
            else:
                logging.error("DEVELOPMENTAL ERROR IN GROUP ASSIGN")
                raise Exception

            logging.info(
                f"{FI.id} {FI.description} has clustered to {FI.adjusted_group}"
            )

        # if group already updated
        else:
            if not (FI.adjusted_group == ""):
                if not (FI.adjusted_group in list_group):
                    logging.warning(f" Clustering result colliding in {FI.id}")

        # return funinfo object with adjusted group, and selected group

        return FI, list_group[0]


# Append outgroup to given group-gene dataset by search matrix
def append_outgroup(V, df_search, gene, group, path, opt):

    logging.info(f"Appending outgroup on group: {group}, Gene: {gene}")

    list_FI = copy.deepcopy(V.list_FI)

    # In multiprocessing, delete V for memory
    del V

    # ready for by sseqid hash, which group to append
    # this time, append adjusted group
    group_dict = {}
    funinfo_dict = {}

    for funinfo in list_FI:
        group_dict[funinfo.hash] = funinfo.adjusted_group
        funinfo_dict[funinfo.hash] = funinfo

    # For non-concatenated analysis
    if gene != "concatenated":
        # generate minimal bitscore cutoff that does not overlaps to query-query bitscore value range
        cutoff_set_df = df_search[df_search["subject_group"] == group]
        try:
            bitscore_cutoff = min(cutoff_set_df["bitscore"])
        except:
            bitscore_cutoff = 999999  # use infinite if failed

        # get result stasifies over cutoff
        cutoff_df = df_search[df_search["bitscore"] < bitscore_cutoff]
        cutoff_df = cutoff_df[cutoff_df["bitscore"] > 0]
        # split that same group to include all to alignment, and left other groups for outgroup selection
        cutoff_df = cutoff_df[cutoff_df["subject_group"] != group]

        # If no or fewer than designated number of outgroup matches to condition, use flexible criteria
        if cutoff_df.groupby(["subject_group"]).count().empty:
            logging.warning(
                f"Not enough outgroup sequences matched for group {group} | gene {gene}. There might be outlier sequence that does not matches to group. Trying flexible cutoff"
            )
            cutoff_df = df_search[df_search["bitscore"] > 0]
            cutoff_df = cutoff_df[cutoff_df["subject_group"] != group]

        elif (
            cutoff_df.groupby(["subject_group"]).count()["sseqid"].max()
            < opt.maxoutgroup
        ):
            logging.warning(
                f"Not enough outgroup sequences matched for group {group} | gene {gene}. There might be outlier sequence that does not matches to group. Trying flexible cutoff"
            )
            cutoff_df = df_search[df_search["bitscore"] > 0]
            cutoff_df = cutoff_df[cutoff_df["subject_group"] != group]

    # in concatenated dataset, upper method makes problem when gene dataset were biased. using old method instead
    else:
        cutoff_df = df_search[df_search["bitscore"] > 0]
        cutoff_df = cutoff_df[cutoff_df["subject_group"] != group]

    # sort by bitscore
    cutoff_df.sort_values(by=["bitscore"], inplace=True, ascending=False)
    # reset index to easily get maximum
    cutoff_df.reset_index(inplace=True, drop=True)

    # iterate until designated number of sequences from the most closest group selected
    # we should get outgroup from columns
    outgroup_dict = {}
    max_cnt = 0
    max_group = ""

    for n, subject_group in enumerate(cutoff_df["subject_group"]):
        # if first sequence of the group found, make new key to dict
        if not (subject_group) in outgroup_dict:
            outgroup_dict[subject_group] = [funinfo_dict[cutoff_df["sseqid"][n]]]
        # if group record already exists, append it
        else:
            if (
                not (funinfo_dict[cutoff_df["sseqid"][n]])
                in outgroup_dict[subject_group]
            ):
                outgroup_dict[subject_group].append(
                    funinfo_dict[cutoff_df["sseqid"][n]]
                )

        # if enough outgroup sequences found while running
        if len(outgroup_dict[subject_group]) >= opt.maxoutgroup:
            text_outgroup_list = "\n ".join(
                [FI.id for FI in outgroup_dict[subject_group]]
            )
            logging.info(
                f"Outgroup [{subject_group}] selected to [{group}]\n {text_outgroup_list}"
            )

            # Get database sequences in ambiguous position
            ambiguous_group = []
            for _group in outgroup_dict:
                if _group != subject_group:
                    ambiguous_group += outgroup_dict[_group]

            return (
                group,
                gene,
                outgroup_dict[subject_group],
                ambiguous_group,
            )
        else:
            if len(outgroup_dict[subject_group]) > max_cnt:
                max_cnt = len(outgroup_dict[subject_group])
                max_group = subject_group

    logging.warning(
        f"Not enough sequences are available for outgroup number {opt.maxoutgroup} in {group}, using '{max_group}' despite of lower number"
    )

    if not (max_group) == "":
        logging.info(
            f"Final outgroup selection for group {group} : {outgroup_dict[max_group]}"
        )

        # Get database sequences in ambiguous position
        ambiguous_group = []
        for _group in outgroup_dict:
            if _group != subject_group:
                ambiguous_group += outgroup_dict[_group]

        return (group, gene, outgroup_dict[max_group], ambiguous_group)
    else:
        logging.warning(f"No outgroup sequence available for {group}")
        return (group, gene, [], [])


def group_cluster_opt_generator(V, opt, path):

    # cluster(FO, df_search, V, path, opt)
    if len(V.list_qr_gene) == 0:
        logging.error(
            "In group_cluster_option_generator, no possible query genes were selected"
        )
        raise Exception

    # For non concatenated analysis or if only one gene exists
    elif len(V.list_qr_gene) <= 1 or V.cSR is None:

        # For caching
        df_group_dict = {}
        qseqid_dict = {}

        # Remove database sequences if queryonly is True
        if opt.queryonly is True:
            list_FI = [FI for FI in V.list_FI if FI.datatype == "query"]
        else:
            list_FI = V.list_FI

        for gene in V.list_qr_gene:
            qseqid_dict[gene] = list(set(V.dict_gene_SR[gene]["qseqid"]))
            df_group_dict[gene] = V.dict_gene_SR[gene].groupby(
                V.dict_gene_SR[gene]["qseqid"]
            )

        list_multigene_FI = []

        for FI in list_FI:
            # If only one gene
            if len(list(FI.seq.keys())) == 1:
                gene = list(FI.seq.keys())[0]
                appropriate_df = df_group_dict[gene].get_group(FI.hash)
                V.opt_cluster.append((FI, appropriate_df, V, path, opt))

            # If no sequence available
            elif len(list(FI.seq.keys())) == 0:
                appropriate_df = None
                V.opt_cluster.append((FI, appropriate_df, V, path, opt))

            else:
                # Get FI with multigene for reporting
                list_multigene_FI.append(FI)
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

        if len(list_multigene_FI) > 0:
            logging.warning(
                f"{' '.join([FI.id for FI in list_multigene_FI])} has multiple genes, but concatenation option not selected"
            )

    # For concatenated analysis
    else:
        # cluster group by concatenated search result
        df_group = V.cSR.groupby(V.cSR["qseqid"])
        list_id = list(set(V.cSR["qseqid"]))
        for FI in V.list_FI:
            if FI.hash in list_id:
                df_search = df_group.get_group(FI.hash)
            else:
                df_search = None
            V.opt_cluster.append((FI, df_search, V, path, opt))

    return V


# opts ready for multithreading in outgroup append
def outgroup_append_opt_generator(V, path, opt):

    opt_append_outgroup = []

    # if concatenated analysis is false
    # Assign different outgroup for each dataset
    for group in V.dict_dataset:
        for gene in V.dict_dataset[group]:
            if gene != "concatenated":
                try:
                    df = V.dict_gene_SR[gene]
                    df_group = df.groupby(df["query_group"])
                    df_group_ = df_group.get_group(group)
                    # Generating outgroup opt for multiprocessing
                    opt_append_outgroup.append((V, df_group_, gene, group, path, opt))
                except:
                    logging.warning(
                        f"{group} / {gene} dataset exists, but cannot append outgroup due to no corresponding search result"
                    )

    # if concatenated analysis is true
    # concatenated
    for group in V.dict_dataset:
        if "concatenated" in V.dict_dataset[group]:
            try:
                df = V.cSR
                df_group = df.groupby(df["query_group"])
                df_group_ = df_group.get_group(group)
                # Generating outgroup opt for multiprocessing
                opt_append_outgroup.append(
                    (V, df_group_, "concatenated", group, path, opt)
                )
            except:
                logging.warning(
                    f"{group} / concatenated dataset exists, but cannot append outgroup due to no corresponding search result"
                )

    return opt_append_outgroup


## Main cluster pipe
def pipe_cluster(V, opt, path):
    # If clustering enabled
    if opt.method.search in ("blast", "mmseqs"):

        logging.info("group clustering")

        # cluster opt generation for multiprocessing
        V = group_cluster_opt_generator(V, opt, path)

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
        for cluster_result in V.rslt_cluster:
            FI = cluster_result[0]
            logging.debug((FI.id, FI.datatype, FI.group, FI.adjusted_group))

        # replace group assigning result
        # collect FI from cluster result
        replace_FI = [r[0] for r in V.rslt_cluster]
        # collect hash
        replace_hash_FI = [FI.hash for FI in replace_FI]
        # maintain not clustered result and append clustered result
        V.list_FI = [
            FI for FI in V.list_FI if not (FI.hash in replace_hash_FI)
        ] + replace_FI
        # For syncyhronizing FI in dict_hash_FI to prevent error
        for FI in replace_FI:
            V.dict_hash_FI[FI.hash] = FI

        V.list_group = list(set([r[1] for r in V.rslt_cluster if (not (r[1] is None))]))

        if opt.queryonly is True:
            for FI in V.list_FI:
                if FI.datatype == "db":
                    FI.adjusted_group = FI.group

            # Update dict_hash_FI
            for FI in V.list_FI:
                V.dict_hash_FI[FI.hash] = FI

        for FI in V.list_FI:
            logging.debug((FI.id, FI.datatype, FI.group, FI.adjusted_group))

    # If not, try to use original groups in tabled format
    else:
        logging.info(
            "[INFO] No searching method designated. Trying to use designated group"
        )
        group_list = get_naive_group(V)
        ## [WIP] Need to make validation process and warn if the input does not have
        for FI in V.list_FI:
            FI.adjusted_group = FI.group

        # Update dict_hash_FI
        for FI in V.list_FI:
            V.dict_hash_FI[FI.hash] = FI

    return V, opt, path


## Main outgroup appending pipeline
def pipe_append_outgroup(V, path, opt):

    opt_append_outgroup = outgroup_append_opt_generator(V, path, opt)

    # run multiprocessing start
    if opt.verbose < 3:
        p = mp.Pool(opt.thread)
        result_append_outgroup = p.starmap(append_outgroup, opt_append_outgroup)
        p.close()
        p.join()

    else:
        # non-multithreading mode for debugging
        result_append_outgroup = [append_outgroup(*o) for o in opt_append_outgroup]

    # append outgroup by running result
    # (group, gene, outgroup, ambiguous_group)
    for result in result_append_outgroup:
        group = result[0]
        gene = result[1]
        outgroup = result[2]
        ambiguous_group = result[3]

        V.dict_dataset[group][gene].list_og_FI = outgroup
        V.dict_dataset[group][gene].list_db_FI += ambiguous_group

    return V, path, opt

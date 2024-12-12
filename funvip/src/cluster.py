import multiprocessing as mp
import pandas as pd
import numpy as np
import shutil
import os, sys, subprocess
from copy import deepcopy
from functools import lru_cache
from time import sleep
from time import time

# from Bio.Blast import NCBIXML
from Bio import SeqIO

import logging
import gc
from funvip.src.ext import mmseqs

import sys


def cluster_unpack(args):
    return cluster(*args)


def batched_generator(generator, batch_size):
    batch = []
    for item in generator:
        batch.append(item)
        if len(batch) == batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


# return list of original group of given FI
def get_naive_group(V):
    return list(set([FI.group for FI in V.list_FI if type(FI.group) == str]))


# Append query_group information column to concatenated search result
def append_query_group(V):
    # For concatenated gene matrix
    group_dict = {}
    for FI in V.list_FI:
        group_dict[FI.hash] = FI.adjusted_group

    V.cSR["query_group"] = V.cSR["qseqid"].apply(lambda x: group_dict.get(x))

    logging.debug(V.cSR["query_group"])

    # Indicate queries without any corresponding group
    for FI in V.list_FI:
        if FI.adjusted_group == "" and not ("noseq" in FI.issues):
            FI.issues.add("nodb")

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

    for FI in V.list_FI:
        # for seq of FI
        for n, seq in enumerate(FI.unclassified_seq):
            try:
                # get each of the dataframe for each FI
                current_df = gene_result_grouped.get_group(f"{FI.hash}_{n}")
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
                FI.update_seq(gene_list[0], seq)
                logging.info(f" Query seq in {FI.id} has assigned to {gene_list[0]}.")
            # if no gene matched, warn it
            elif gene_count == 0:
                logging.warning(
                    f" Query seq in {FI.id} cannot be assigned to any gene. Check sequence. Skipping {FI.id}"
                )
            elif gene_count >= 2:
                logging.warning(f" Query seq in {FI.id} has multiple matches to gene.")
                FI.update_seq(gene_list[0], seq)

            else:
                logging.error("DEVELOPMENTAL ERROR IN GENE ASSIGN")
                raise Exception
    return V


def cluster(FI, V_list_group, V_cSR, path, opt):
    # Reduce memory by focusing on relevant rows
    df_search = V_cSR[V_cSR["qseqid"] == FI.hash]

    if df_search.empty:
        # If confident is False and FI datatype is "db"
        if opt.confident is False and FI.datatype is "db":
            logging.warning(f"No adjusted_group assigned to {FI}")
        FI.adjusted_group = FI.group
        return FI

    # For db sequence with group, retain it
    elif not (FI.group == "") and FI.datatype == "db":
        FI.adjusted_group = FI.group
        return FI

    # Update group if sequence doesn't have one
    else:
        df_search = df_search.sort_values(by=["bitscore"], ascending=False)

        # Apply cutoff filter to reduce DataFrame size
        cutoff = df_search["bitscore"].iloc[0] * opt.cluster.cutoff
        cutoff_df = df_search[df_search["bitscore"] > cutoff]

        # Clear unused data
        del df_search
        gc.collect()

        # Extract group information
        unique_groups = set(cutoff_df["subject_group"])
        group_count = len(unique_groups)

    if FI.adjusted_group == "":
        if group_count == 1:
            FI.adjusted_group = unique_groups.pop()
        elif group_count == 0:
            logging.warning(
                f"Query seq in {FI.id} cannot be assigned to group. Check sequence."
            )
        elif group_count >= 2:
            logging.warning(
                f"Query seq in {FI.id} has multiple matches to groups within {opt.cluster.cutoff} cutoff: {list(unique_groups)}"
            )

            # Among the multiple matches, get the best group
            FI.adjusted_group = cutoff_df["subject_group"].iloc[0]

        else:
            logging.error("DEVELOPMENTAL ERROR IN GROUP ASSIGN")
            raise Exception

        logging.info(f"{FI.id} has clustered to {FI.adjusted_group}")

    # if group already updated
    else:
        if not (FI.adjusted_group in unique_groups):
            logging.warning(f"Clustering result collides for {FI.id}")

    if unique_groups:
        return FI
    else:
        return FI


### Append outgroup to given group-gene dataset by search matrix
def append_outgroup(V_list_FI, df_search, gene, group, path, opt):
    logging.info(f"Appending outgroup on group: {group}, Gene: {gene}")
    list_FI = V_list_FI

    # In multiprocessing, delete V to reduce memory consumption
    # del V
    # gc.collect()

    # ready for by sseqid hash, which group to append
    # this time, append adjusted group
    group_dict = {}
    FI_dict = {}

    for FI in list_FI:
        group_dict[FI.hash] = FI.adjusted_group
        FI_dict[FI.hash] = FI

    # For non-concatenated analysis
    # if gene != "concatenated":
    # generate minimal bitscore cutoff that does not overlaps to query-query bitscore value range

    ## For getting inner group
    cutoff_set_df = df_search[df_search["subject_group"] == group]
    try:
        # offset will prevent selecting outgroup too close to ingroup, which may confuse the monophyly of outgroup
        bitscore_cutoff = max(
            1, min(cutoff_set_df["bitscore"]) - opt.cluster.outgroupoffset
        )
    except:
        bitscore_cutoff = 999999  # use infinite if failed

    # print(f"Ingroup cutoff {bitscore_cutoff} selected for group {group} gene {gene}")

    ## get result stasifies over cutoff
    # outgroup should be outside of ingroup
    cutoff_df = df_search[df_search["bitscore"] < bitscore_cutoff]

    # Remove malformat result, which bitscore is under 0
    cutoff_df = cutoff_df[cutoff_df["bitscore"] > 0]

    # split that same group to include all to alignment, and leave other groups for outgroup selection
    cutoff_df = cutoff_df[cutoff_df["subject_group"] != group]

    ## For ambiugous database, mostly because of contaminated database
    # For each of the input, should use different cutoff
    ambiguous_db = set()
    for qseqid, _df in cutoff_set_df.groupby(["qseqid"]):
        # Select dataframe corresponding to current qseqid
        df_qseqid = df_search[df_search["qseqid"] == qseqid]
        """
        print(
            f"Ambiguous ingroup cutoff selected for query {qseqid} group {group} gene {gene} cutoff {min(list(_df['bitscore']))}"
        )
        """
        # Get the list of subjects, which is closer than furtest ingroup
        ambiguous_df = df_qseqid[df_qseqid["bitscore"] >= min(list(_df["bitscore"]))]
        # Within the furthest match, get possible ingroups with ambiguous group
        ambiguous_df = ambiguous_df[ambiguous_df["subject_group"] != group]
        # Add inner ambiugities to ambiguous db
        ambiguous_db.update([FI_dict[i] for i in list(ambiguous_df["sseqid"])])

    ambiguous_db = list(ambiguous_db)

    # If no or fewer than designated number of outgroup matches to condition, use flexible criteria
    if cutoff_df.groupby(["subject_group"]).count().empty:
        logging.warning(
            f"Not enough outgroup sequences matched for group {group} | gene {gene}. There might be outlier sequence that does not matches to group. Trying flexible cutoff"
        )
        cutoff_df = df_search[df_search["bitscore"] > 0]
        cutoff_df = cutoff_df[cutoff_df["subject_group"] != group]

    elif cutoff_df.groupby(["subject_group"]).count()["sseqid"].max() < opt.maxoutgroup:
        logging.warning(
            f"Not enough outgroup sequences matched for group {group} | gene {gene}. There might be outlier sequence that does not matches to group. Trying flexible cutoff"
        )
        cutoff_df = df_search[df_search["bitscore"] > 0]
        cutoff_df = cutoff_df[cutoff_df["subject_group"] != group]

    # Garbage collection to reduce memory consumption in this process
    del df_search
    gc.collect()

    # sort by bitscore
    cutoff_df.sort_values(by=["bitscore"], inplace=True, ascending=False)
    # reset index to easily get maximum
    cutoff_df.reset_index(inplace=True, drop=True)

    # iterate until designated number of sequences from the most closest group selected
    # we should get outgroup from columns
    # outgroup_dict = {group1 : [FI1, FI2, ...]}
    outgroup_dict = {}
    max_cnt = 0
    max_group = ""

    for n, subject_group in enumerate(cutoff_df["subject_group"]):
        ## Check if outgroup sequence that we're going to use actually exists
        # If gene is not "concatenated", the outgroup sequence should actually exists
        # keep in mind if case of certain gene is blank

        cond1 = gene in FI_dict[cutoff_df["sseqid"][n]].seq
        if cond1 is True:
            cond1 = FI_dict[cutoff_df["sseqid"][n]].seq[gene] != ""

        # If gene is "concatenated", any of the genes should exists
        cond2 = (
            gene == "concatenated"
            and len(FI_dict[cutoff_df["sseqid"][n]].seq.keys()) > 0
        )
        if cond1 or cond2:
            # if first sequence belonging to the group found, make new key to dict
            if not (subject_group) in outgroup_dict:
                outgroup_dict[subject_group] = [FI_dict[cutoff_df["sseqid"][n]]]
            # if group key already exists in dict, append it
            else:
                if (
                    not (FI_dict[cutoff_df["sseqid"][n]])
                    in outgroup_dict[subject_group]
                ):
                    outgroup_dict[subject_group].append(FI_dict[cutoff_df["sseqid"][n]])

            # if enough outgroup sequences found while running
            if len(outgroup_dict[subject_group]) >= opt.maxoutgroup:
                text_outgroup_list = "\n ".join(
                    [FI.id for FI in outgroup_dict[subject_group]]
                )
                logging.info(
                    f"Outgroup [{subject_group}] selected to [{group}]\n {text_outgroup_list}"
                )

                return (
                    group,
                    gene,
                    outgroup_dict[subject_group],
                    ambiguous_db,
                )
            else:
                if len(outgroup_dict[subject_group]) > max_cnt:
                    max_cnt = len(outgroup_dict[subject_group])
                    max_group = subject_group

    # If not enough outgroup sequences found while running
    logging.warning(
        f"Not enough sequences are available for outgroup number {opt.maxoutgroup} in {group}, using '{max_group}' despite of lower number"
    )

    # If outgroup are selected
    if not (max_group) == "":
        logging.info(
            f"Final outgroup selection for group {group} : {outgroup_dict[max_group]}"
        )

        return (group, gene, outgroup_dict[max_group], ambiguous_db)

    # If outgroup cannot be selected
    else:
        logging.warning(
            f"No outgroup sequence available for {group}. Try to use\n A. Higher --cluster-evalue\n B. Lower --outgroupoffset\n C. Add closer sequence of {group} to database"
        )
        return (group, gene, [], [])


def group_cluster_opt_generator(V, opt, path):
    # opt_cluster = []

    # cluster(FO, df_search, V, path, opt)
    if len(V.list_qr_gene) == 0:
        logging.error(
            "In group_cluster_option_generator, no available query genes were selected"
        )
        raise Exception

    # For concatenated analysis
    else:
        # Use a set for faster membership testing
        set_id = set(V.cSR["qseqid"])

        # Pre-select the relevant DataFrame slice once
        cSR_subset = V.cSR[["qseqid", "bitscore", "subject_group"]]

        # Filter `V.list_FI` to relevant entries
        relevant_FIs = (FI for FI in V.list_FI if FI.hash in set_id)

        # Yield results as a generator
        for FI in relevant_FIs:
            yield (
                FI,
                V.list_group,  # Assuming V.list_group is small and static
                cSR_subset,  # Pre-selected DataFrame slice
                path,
                opt,
            )


# opts ready for multithreading in outgroup append
def outgroup_append_opt_generator(V, path, opt):
    opt_append_outgroup = []

    #### Pararellize this part
    # Assign different outgroup for each dataset

    # if concatenated analysis is true
    # concatenated
    for group in V.dict_dataset:
        if "concatenated" in V.dict_dataset[group]:
            try:
                df = V.cSR
                df_group = df.groupby(df["query_group"])
                df_group_ = df_group.get_group(group)
                # Generating outgroup opt for multiprocessing
                for gene in V.dict_dataset[group]:
                    opt_append_outgroup.append(
                        (V.list_FI, df_group_, gene, group, path, opt)
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

        rslt_cluster = []

        # cluster opt generation for multiprocessing
        # (FI, V, path, opt)
        opt_cluster = group_cluster_opt_generator(V, opt, path)

        # print(f"opt_cluster : {sys.getsizeof(opt_cluster)}")

        batch_size = opt.thread * 100
        opt_cluster_batches = batched_generator(opt_cluster, batch_size)

        # run multiprocessing start
        if opt.verbose < 3:
            for batch in opt_cluster_batches:
                batch_list = list(batch)

                with mp.Pool(opt.thread) as p:
                    rslt_cluster.extend(p.starmap(cluster, batch_list))

        else:
            # non-multithreading mode for debugging
            rslt_cluster = [cluster(*o) for o in opt_cluster]
        # gather cluster result
        for cluster_result in rslt_cluster:
            FI = cluster_result
            logging.debug((FI.id, FI.datatype, FI.group, FI.adjusted_group))

        # replace group assigning result
        # collect FI from cluster result
        # somethings been duplicated here
        replace_FI = [r for r in rslt_cluster]

        # collect hash
        replace_hash_FI = [FI.hash for FI in replace_FI]

        # maintain not clustered result and append clustered result
        V.list_FI = [FI for FI in V.list_FI if not (FI.hash in replace_hash_FI)]
        V.list_FI += replace_FI

        # For syncyhronizing FI in dict_hash_FI to prevent error
        for FI in replace_FI:
            V.dict_hash_FI[FI.hash] = FI

        V.list_group = list(
            set(
                [
                    r.adjusted_group
                    for r in rslt_cluster
                    if (not (r.adjusted_group == ""))
                ]
            )
        )

        if opt.queryonly is True:
            for FI in V.list_FI:
                if FI.datatype == "db":
                    FI.adjusted_group = FI.group

            # Update dict_hash_FI
            for FI in V.list_FI:
                V.dict_hash_FI[FI.hash] = FI

        # For debugging
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
        # Parsing result
        group = result[0]
        gene = result[1]
        outgroup = result[2]
        ambiguous_group = result[3]

        # Add outgroup and ambiguous groups to dataset
        # Ambiguous groups are strains locating between outgroup and ingroups, so cannot be decided
        if len(outgroup) == 0 and len(ambiguous_group) == 0:
            logging.warning(
                f"Removing {group} {gene} from analysis because outgroup cannot be selected"
            )
            V.dict_dataset[group].pop(gene, None)

        else:
            V.dict_dataset[group][gene].list_og_FI = outgroup
            # Add ambiguous group to FI
            if opt.ambiguous is True:
                V.dict_dataset[group][gene].list_db_FI += ambiguous_group
            # Add outgroup and db in to dict_hash_FI

    groups = deepcopy(list(V.dict_dataset.keys()))

    for group in groups:
        try:
            if len(V.dict_dataset[group]) == 0:
                logging.warning(
                    f"Removing {group} from analysis because outgroup cannot be selected to all genes"
                )
                V.dict_dataset.pop(group, None)
        except:
            pass

    return V, path, opt

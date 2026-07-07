# Performing multiple tree interpretation
from ete4 import Tree
from funvip.src import tree_interpretation
from funvip.src.tool import initialize_path, get_genus_species
from funvip.src.tool import sizeof_fmt
from funvip.src.hasher import encode, decode
from funvip.src.reporter import Singlereport
from funvip.src.exceptions import TreeError
import traceback
from copy import deepcopy
import pandas as pd
import re
import sys
import os
import shutil
import psutil
import logging
import multiprocessing as mp
from time import time


### For single dataset
# Input : out, group, gene, V, path, opt
# Whole-run FI collections (V.dict_hash_FI / V.list_FI) shared with
# interpretation-pool workers via fork copy-on-write, set before the Pool is
# created, instead of being pickled into every (group, gene) task tuple.
_INTERP_SHARED = {}


def pipe_module_tree_interpretation(
    out,
    group,
    gene,
    V_tup_genus,
    hash_dict,
    query_list,
    outgroup,
    partition,
    path,
    opt,
):
    # time_start = time()

    # Read the whole-run FI collections from the fork-inherited shared store rather
    # than receiving a freshly pickled copy of the entire universe per task.
    funinfo_dict = _INTERP_SHARED["funinfo_dict"]
    funinfo_list = _INTERP_SHARED["funinfo_list"]

    # for unexpectedly included sequence during clustering
    db_list = list(
        set([FI for FI in funinfo_list if FI.datatype == "db"])
        - set(outgroup)
        - set(query_list)
    )
    genus_list = V_tup_genus

    # For get_genus_species
    initialize_path(path)

    # Tree name selection for tree construction software
    tree_name = f"{path.out_tree}/hash/hash_{opt.runname}_{group}_{gene}.nwk"
    logging.debug(f"{tree_name} entered tree interpretation")

    # Check validity of file while importing
    if os.path.isfile(tree_name):
        try:
            # If iqtree, missing supports are not shown
            # ete4 autodetects newick format regardless of software
            Tree(tree_name)
        except Exception as e:
            raise TreeError(f"failed to parse tree file {tree_name}: {e}") from e
    else:
        raise TreeError(f"cannot find tree file {tree_name}")

    # initialize before analysis
    Tree_style = tree_interpretation.Tree_style()

    # Read tree
    tree_info = tree_interpretation.Tree_information(
        tree_name, Tree_style, group, gene, opt
    )

    # Give necessary variables parsed from dataset
    tree_info.db_list = tuple(db_list)
    tree_info.query_list = tuple(query_list)
    tree_info.outgroup = tuple(outgroup)
    tree_info.funinfo_dict = funinfo_dict

    # Main phase
    # calculate zero distance with alignment
    if gene == "concatenated":
        tree_info.calculate_zero(
            alignment_file=f"{path.out_alignment}/hash/{opt.runname}_hash_trimmed_{group}_{gene}.fasta",
            gene=gene,
            partition_dict=partition,
        )
    else:
        tree_info.calculate_zero(
            alignment_file=f"{path.out_alignment}/hash/{opt.runname}_hash_trimmed_{group}_{gene}.fasta",
            gene=gene,
            partition_dict=None,
        )

    # print(f"Calculate zero {time() - time_start}")

    # Reroot outgroup and save original tree into image
    tree_info.reroot_outgroup(
        f"{path.out_tree}/hash_{opt.runname}_{group}_{gene}_original.svg"
    )

    # print(f"Reroot outgroup {time() - time_start}")

    # Decode hash of image
    # Should work more on non-safe characters
    tree_hash_dict = encode(funinfo_list, newick=True)

    decode(
        tree_hash_dict,
        f"{path.out_tree}/hash_{opt.runname}_{group}_{gene}_original.svg",
        f"{path.out_tree}/{opt.runname}_{group}_{gene}_original.svg",
        svg=True,
    )

    # print(f"Decode {time() - time_start}")

    # In validation mode, use original sp. number
    if opt.mode == "validation":
        tree_info.reserve_sp()

    # Reconstruct flat branches if option given
    if opt.solveflat is True:
        # ete4: copy("newick") uses parser=1 which drops support values;
        # use deepcopy to preserve branch support for visualization.
        _clade = tree_info.t.copy("deepcopy")
        for _n in _clade.traverse():
            if _n.dist is None:
                _n.dist = 0.0
            # ete4 compat: ete3 DEFAULT_SUPPORT=1.0 became 100 after scale
            # conversion; ete4 leaves/root get support=None. Normalize only
            # INTERNAL nodes None→100 (when scale conversion occurred) so
            # intermediate nodes created by reconstruct/solve_flat carry
            # support=100 like ete3. Leaves must stay None: concat_clade maps
            # None→1 for them, which is < bscutoff (correct — ete3 also sets
            # leaf DEFAULT_SUPPORT=1.0 → 1 after newick round-trip → not shown).
            if not _n.is_leaf and _n.support is None and tree_info.support_scaled:
                _n.support = 100
        tree_info.t = tree_info.reconstruct(
            clade=_clade, gene=gene, opt=opt
        )
        # reconstruct() replaced tree_info.t with a shallow rebuild, but
        # tree_info.outgroup_clade still points at a node of the PRE-reconstruct
        # tree. reroot_outgroup's resolve_polytomy() combs a conserved gene's giant
        # star (5.8S: ~3670-leaf polytomy) into a ~3681-deep chain, and pickling any
        # ete4 node drags in its whole tree via .up/.children -- so this stale ref
        # re-inflates that deep tree when the worker pickles tree_info back to the
        # parent (multiprocessing), raising RecursionError -> MaybeEncodingError.
        # outgroup_clade is unused after rerooting, so drop it. (Third part of the
        # 5.8S fix, with seperate_clade deepcopy + balanced concat_all.)
        tree_info.outgroup_clade = None

    # print(f"Reconstruct {time() - time_start}")

    # reorder tree for pretty look (ete3-compatible tie-breaking)
    tree_interpretation._ladderize_ete3_compat(tree_info.t)

    # print(f"Ladderize {time() - time_start}")

    # save current status into save version of tree
    # Is not currently used
    # tree_info.t_publish = deepcopy(tree_info.t)

    # Search tree and delimitate species.
    # Enable subtree taxon-count memoization for this (static, post-solve_flat) tree only.
    tree_interpretation._tc_cache_begin()
    try:
        tree_info.tree_search(tree_info.t, gene)
    finally:
        tree_interpretation._tc_cache_end()

    # print(f"Tree search {time() - time_start}")

    # Move original newick and replace with adjusted ones
    shutil.move(
        f"{path.out_tree}/{opt.runname}_{group}_{gene}.nwk",
        f"{path.out_tree}/{opt.runname}_{group}_{gene}_original.nwk",
    )
    tree_info.t.write(
        outfile=f"{path.out_tree}/{opt.runname}_{group}_{gene}.nwk"
    )
    decode(
        tree_hash_dict,
        f"{path.out_tree}/{opt.runname}_{group}_{gene}.nwk",
        f"{path.out_tree}/{opt.runname}_{group}_{gene}.nwk",
        newick=True,
    )

    # print(f"Tree interpretation ended {time() - time_start}")

    return tree_info


### synchronize sp. numbers from multiple dataset
# to use continuous sp numbers over trees
# Seperated from multithreading, because this step should traverse over multiple trees, therefore cannot be done simultaneously
def synchronize(V, path, tree_info_list):
    # Gets hash dict, and returns taxon name of hash_dict
    # Generate final taxon name for synchronizing
    def get_new_taxon(hash_list, hash_taxon_dict):
        # Get candidate taxons from hash_taxon_dict
        taxon_candidates = set()
        for _hash in hash_list:
            if _hash in hash_taxon_dict:
                taxon_candidates.add(hash_taxon_dict[_hash])
            else:
                logging.debug(
                    f"{_hash} does not seems to be analyzed from concatenated dataset"
                )

        list_taxon_candidates = sorted(list(taxon_candidates))

        # Merge genus
        genus = "/".join(sorted(list(set(t[0][0] for t in list_taxon_candidates))))

        # Merge species
        species_list = []

        clade_cnt_set = set()
        for t in list_taxon_candidates:
            # If the taxon is unique
            if not (t[0], 2) in hash_taxon_dict.values():
                species_list.append(t[0][1])
            else:
                species_list.append(f"{t[0][1]} {t[1]}")
            clade_cnt_set.add(t[1])

        # Work with species with numbers
        dict_species = {}
        for s in species_list:
            splited_species = s.split(" ")
            try:
                # Collect with numbers
                int(splited_species[-1])
                if not (" ".join(splited_species[:-1]) in dict_species):
                    dict_species[" ".join(splited_species[:-1])] = [
                        int(splited_species[-1])
                    ]
                else:
                    dict_species[" ".join(splited_species[:-1])].append(
                        int(splited_species[-1])
                    )
            except (ValueError, IndexError):
                # species label does not end in an integer -> treat as unnumbered
                dict_species[s] = [0]

        species = ""

        for key in sorted(list(dict_species.keys())):
            if len(set(dict_species[key]) - set([0])) == 0:
                species += key
                species += "/"
            else:
                species_numbers = [
                    str(x) for x in sorted(list(set(dict_species[key]) - set([0])))
                ]
                species += key
                species += " "
                species += "/".join(species_numbers)

        # Remove last slash
        if species.endswith("/"):
            species = species[:-1]

        if len(clade_cnt_set) == 1:
            clade_cnt = list(clade_cnt_set)[0]
        else:
            clade_cnt = 0

        return (genus, species), clade_cnt
        ### End of get_new_taxon

    ## Initialize
    # get available groups per genus
    tree_info_dict = {}
    # hash : corresponding taxon
    hash_taxon_dict = {}
    # genus : cnt, counting sp. numbers
    sp_cnt_dict = {}

    # To synchronize sp. number by genus, generate by-group dataset
    for tree_info in tree_info_list:
        if not (tree_info.group in tree_info_dict):
            tree_info_dict[tree_info.group] = {tree_info.gene: tree_info}
        elif not (tree_info.gene in tree_info_dict[tree_info.group]):
            tree_info_dict[tree_info.group][tree_info.gene] = tree_info
        else:
            raise TreeError(
                f"duplicated tree_info for group {tree_info.group} gene {tree_info.gene}"
            )

    # Memoize iterative calling
    # For each group list
    valid_hash_dict = {}
    # DB
    for group in tree_info_dict:
        valid_hash_dict[group] = [
            _hash
            for _hash in V.dict_hash_FI
            if V.dict_hash_FI[_hash].datatype == "db"
            and V.dict_hash_FI[_hash].adjusted_group == group
        ]
    # Query
    query_hash_list = [
        _hash for _hash in V.dict_hash_FI if V.dict_hash_FI[_hash].datatype == "query"
    ]

    ## Starting with concatenated
    # In priority, count corresponding group taxa first
    for group in tree_info_dict:
        if "concatenated" in tree_info_dict[group]:
            # Catch concatenated tree
            tree_info = tree_info_dict[group]["concatenated"]
            # Get list of hash in interest
            valid_hash_list = valid_hash_dict[group]
            for taxon in tree_info.collapse_dict:
                # Get all monophyletic clades from tree and make it list
                clade_list = tree_info.collapse_dict[taxon]
                for n, clade in enumerate(clade_list):
                    # list of hash in clade
                    hash_list = [leaf[0] for leaf in clade.leaf_list]
                    # If any of the leaf consisting clade is included to valid_hash_list
                    if any(_h in valid_hash_list for _h in hash_list):
                        for _hash in hash_list:
                            if not (_hash) in hash_taxon_dict:
                                if len(clade_list) == 1:
                                    hash_taxon_dict[_hash] = (taxon, 0)
                                else:
                                    hash_taxon_dict[_hash] = (taxon, n + 1)
                            elif _hash in hash_taxon_dict and hash_taxon_dict[
                                _hash
                            ] != (taxon, n):
                                logging.debug(
                                    f"{_hash} collided while putting in hash_taxon_dict"
                                )

                    # If leaf consisting with only queries, that won't collide with other groups
                    # This is about new species clade
                    if all(_h in query_hash_list for _h in hash_list):
                        # If this clade is first sp. species for the genus, start counting sp. number
                        if not (taxon[0] in sp_cnt_dict):
                            sp_cnt_dict[taxon[0]] = 1
                        for _hash in hash_list:
                            if not (_hash) in hash_taxon_dict:
                                hash_taxon_dict[_hash] = (
                                    (
                                        taxon[0],
                                        f"sp. {sp_cnt_dict[taxon[0]]}",
                                    ),
                                    0,
                                )

                            elif (
                                _hash in hash_taxon_dict
                                and hash_taxon_dict[_hash] != taxon
                            ):
                                logging.debug(
                                    f"{_hash} collided while putting in hash_taxon_dict"
                                )

                        sp_cnt_dict[taxon[0]] += 1

    # Next, taxa that doesn't belongs to any of the group
    # concatenated first
    for group in tree_info_dict:
        if "concatenated" in tree_info_dict[group]:
            tree_info = tree_info_dict[group]["concatenated"]

            all_hash = [
                _hash
                for _hash in V.dict_hash_FI
                if V.dict_hash_FI[_hash].datatype == "db"
                and V.dict_hash_FI[_hash].adjusted_group != group
            ]

            # Get list of hash not in interest
            invalid_hash_list = list(set(all_hash) - set(valid_hash_dict[group]))
            for taxon in tree_info.collapse_dict:
                clade_list = tree_info.collapse_dict[taxon]
                for n, clade in enumerate(clade_list):
                    hash_list = [leaf[0] for leaf in clade.leaf_list]
                    # If the hash has not been counted in any of the tree,
                    if not (any(_h in invalid_hash_list for _h in hash_list)):
                        for _hash in hash_list:
                            if not (_hash in hash_taxon_dict):
                                logging.debug(
                                    f"New hash: {_hash} {_hash.adjusted_group} from {group}"
                                )
                                hash_taxon_dict[_hash] = (taxon, n)
                            elif _hash in hash_taxon_dict and hash_taxon_dict[
                                _hash
                            ] != (taxon, n):
                                logging.debug(
                                    f"{_hash} collided while putting in hash_taxon_dict. Tried to put {(taxon, n+1)}, but existing {hash_taxon_dict[_hash]}"
                                )

    # Then, non-concatenated
    for group in tree_info_dict:
        # A group whose concatenated tree failed interpretation (dropped by the
        # per-item guard) can still have gene entries here; skip it rather than
        # KeyError on the missing "concatenated" key and abort the whole run.
        if "concatenated" not in tree_info_dict[group]:
            continue
        for gene in tree_info_dict[group]:
            if gene != "concatenated":
                tree_info = tree_info_dict[group]["concatenated"]
                # Get list of hash in interest
                valid_hash_list = valid_hash_dict[group]
                for taxon in tree_info.collapse_dict:
                    clade_list = tree_info.collapse_dict[taxon]
                    for n, clade in enumerate(clade_list):
                        hash_list = [leaf[0] for leaf in clade.leaf_list]
                        # If the hash has not been counted in any of the tree,
                        if not (any(_h in valid_hash_list for _h in hash_list)):
                            for _hash in hash_list:
                                if not (_hash in hash_taxon_dict):
                                    hash_taxon_dict[_hash] = (taxon, n + 1)
                                elif _hash in hash_taxon_dict and hash_taxon_dict[
                                    _hash
                                ] != (taxon, n + 1):
                                    logging.debug(
                                        f"{_hash} collided while putting in hash_taxon_dict. Tried to put {(taxon, n+1)}, but existing {hash_taxon_dict[_hash]}"
                                    )

                                # print(group, _hash, taxon, n + 1)

    """
    # Now update from concatenated
    # Remove original taxon, and add by clade taxon
    for group in tree_info_dict:
        for gene in tree_info_dict[group]:
            if gene == "concatenated":
                tree_info = tree_info_dict[group][gene]
                # Before taxon, after taxon update list
                remove_list = []  # [taxon1, taxon2, taxon3 ...]
                add_list = {}  # [taxon1 : [clade1], taxon2 : [clade2] ...]
                for taxon in tree_info.collapse_dict:
                    clade_list = tree_info.collapse_dict[taxon]
                    remove_list.append(taxon)
                    for clade in clade_list:
                        hash_list = [leaf[0] for leaf in clade.leaf_list]
                        new_taxon, clade_cnt = get_new_taxon(hash_list, hash_taxon_dict)
                        clade.taxon = new_taxon
                        clade.clade_cnt = clade_cnt
                        if not (new_taxon in add_list):
                            add_list[new_taxon] = [clade]
                        else:
                            add_list[new_taxon].append(clade)

                # Remove previous taxon
                for taxon in remove_list:
                    tree_info.collapse_dict.pop(taxon)

                # Add synchronized taxon
                for taxon in add_list:
                    tree_info.collapse_dict[taxon] = add_list[taxon]
            else:
                tree_info = tree_info_dict[group][gene]
                # Before taxon, after taxon update list
                remove_list = []  # [taxon1, taxon2, taxon3 ...]
                add_list = {}  # [taxon1 : [clade1], taxon2 : [clade2] ...]
                for taxon in tree_info.collapse_dict:
                    clade_list = tree_info.collapse_dict[taxon]
                    remove_list.append(taxon)
                    for clade in clade_list:
                        hash_list = [leaf[0] for leaf in clade.leaf_list]
                        new_taxon, clade_cnt = get_new_taxon(hash_list, hash_taxon_dict)
                        clade.taxon = new_taxon
                        clade.clade_cnt = clade_cnt
                        if not (new_taxon in add_list):
                            add_list[new_taxon] = [clade]
                        else:
                            add_list[new_taxon].append(clade)

                # Remove previous taxon
                for taxon in remove_list:
                    tree_info.collapse_dict.pop(taxon)

                # Add synchronized taxon
                for taxon in add_list:
                    tree_info.collapse_dict[taxon] = add_list[taxon]


    """

    # raise Exception

    # Return sp number fixed tree_info_list
    return tree_info_list


### Visualization after synchronization
def pipe_module_tree_visualization(
    tree_info,
    V_tup_genus,
    V_dict_hash_FI,
    path,
    opt,
):
    # time_start = time()

    ######### Fix collapse_dict.keys()
    # V.tup_genus
    # V.dict_hash_FI

    group = tree_info.group
    gene = tree_info.gene
    genus_list = list(V_tup_genus)
    genus_list.append("AMBIGUOUSGENUS")
    genus_list = tuple(genus_list)

    del V_tup_genus

    # Collapse tree branches for visualization
    taxon_string_dict = tree_info.collapse_tree()

    # print(f"Visualize Collapse tree {time() - time_start}")

    # print(f"taxon_string_list | {group} {gene}:\n {taxon_string_list}\n")

    # Polish tree image
    tree_info.polish_image(
        f"{path.out_tree}/{opt.runname}_{group}_{gene}.svg",
        taxon_string_dict,
    )

    # print(f"Visualize polish image {time() - time_start}")

    # sort taxon order
    list_taxon_1 = [
        taxon
        for taxon in tree_info.collapse_dict.keys()
        if not (taxon[1].startswith("sp."))
    ]

    list_taxon_2 = [
        taxon for taxon in tree_info.collapse_dict.keys() if taxon[1].startswith("sp.")
    ]
    list_taxon_1.sort(key=lambda x: x[1])
    list_taxon_2.sort(key=lambda x: x[1])
    list_taxon = list_taxon_1 + list_taxon_2

    # Declare report collection
    report_list = []
    for taxon in list_taxon:
        # If only one taxon exists, enumerate does not work properly
        if len(tree_info.collapse_dict[taxon]) <= 1:
            collapse_info = tree_info.collapse_dict[taxon][0]
            # Get each of the leaf result to report
            for leaf in collapse_info.leaf_list:
                report = Singlereport()
                report.id = V_dict_hash_FI[leaf[0]][0]
                report.hash = V_dict_hash_FI[leaf[0]][1]
                report.update_group(V_dict_hash_FI[leaf[0]][2])
                report.update_group_analysis(group)
                report.update_gene(gene)
                report.update_species_original(
                    get_genus_species(leaf[2], genus_list=genus_list)
                )
                # joining genus and species
                report.update_species_assigned(" ".join(taxon))
                # report.update_species_assigned(taxon[1])
                report.ambiguous = collapse_info.clade_cnt
                report.flat = collapse_info.flat

                report_list.append(report)

        # If more than one taxon exists,
        else:
            for n, collapse_info in enumerate(tree_info.collapse_dict[taxon]):
                for leaf in collapse_info.leaf_list:
                    report = Singlereport()
                    report.id = V_dict_hash_FI[leaf[0]][0]
                    report.hash = V_dict_hash_FI[leaf[0]][1]
                    report.update_group(V_dict_hash_FI[leaf[0]][2])
                    report.update_group_analysis(group)
                    report.update_gene(gene)
                    report.update_species_original(
                        get_genus_species(leaf[2], genus_list=genus_list)
                    )
                    report.update_species_assigned((f"{taxon[0]} {taxon[1]} {n+1}"))

                    report.ambiguous = collapse_info.clade_cnt
                    report.flat = collapse_info.flat

                    report_list.append(report)

        # print(f"Visualize report {time() - time_start}")

    """
    for name, size in sorted(
        ((name, sys.getsizeof(value)) for name, value in list(locals().items())),
        key=lambda x: -x[1],
    )[:10]:
        print("{:>30}: {:>8}".format(name, sizeof_fmt(size)))
    print("==============================")
    """

    if opt.verbose >= 3:
        logging.debug(f"End of pipe module visualization")
        process = psutil.Process(os.getpid())
        memory_info = process.memory_info()
        logging.debug(f"RAM usage: {memory_info.rss / 1000 / 1000} MB")

    # raise Exception

    return report_list


### For all datasets, multiprocessing part
def _safe_pipe_module_tree_interpretation(*args):
    # Per-item guard: a crash in one (group, gene) must not kill the whole batch
    # (starmap re-raises the first worker exception -> no result.csv for the 900+
    # genera that DID succeed). Log and skip the bad one instead.
    try:
        return pipe_module_tree_interpretation(*args)
    except Exception as e:
        group = args[1] if len(args) > 1 else "?"
        gene = args[2] if len(args) > 2 else "?"
        logging.error(
            f"[TREE INTERPRETATION FAILED] {group} {gene}: {e!r}\n{traceback.format_exc()}"
        )
        return None


def _safe_pipe_module_tree_visualization(*args):
    try:
        return pipe_module_tree_visualization(*args)
    except Exception as e:
        ti = args[0] if args else None
        gg = (
            f"{getattr(ti, 'group', '?')} {getattr(ti, 'gene', '?')}"
            if ti is not None
            else "?"
        )
        logging.error(
            f"[TREE VISUALIZATION FAILED] {gg}: {e!r}\n{traceback.format_exc()}"
        )
        return None


def pipe_tree_interpretation(V, path, opt):
    # Generate tree_interpretation opt to run
    # tree_interpretation_opt = []

    # Reset bygene_species for rerun
    for key in V.dict_hash_FI:  # funinfo_dict
        for gene in V.dict_hash_FI[key].bygene_species:
            V.dict_hash_FI[key].bygene_species[gene] = V.dict_hash_FI[key].ori_species

    # common variable preparation
    funinfo_dict = V.dict_hash_FI
    funinfo_list = V.list_FI
    hash_dict = V.dict_hash_name

    # Share the whole-run FI collections with pool workers via fork copy-on-write
    # (must be set before the Pool below is created) instead of pickling them into
    # every task; workers read them from _INTERP_SHARED.
    _INTERP_SHARED["funinfo_dict"] = funinfo_dict
    _INTERP_SHARED["funinfo_list"] = funinfo_list

    # Generate options using generator
    def generate_interpretation_opt():
        # make option variables
        for group in V.dict_dataset:
            partition = V.partition[group]
            for gene in V.dict_dataset[group]:
                query_list = V.dict_dataset[group][gene].list_qr_FI
                outgroup = V.dict_dataset[group][gene].list_og_FI
                logging.debug(f"pipe_tree_interpretation {group} {gene}")
                # Condition 1 : draw all trees
                cond1 = opt.queryonly is False
                # Condition 2 : When query included
                cond2 = len(V.dict_dataset[group][gene].list_qr_FI) > 0
                # Condition 3 : When any of the branches of the tree is valid in concatenated analysis
                cond3 = (
                    len(V.dict_dataset[group]["concatenated"].list_qr_FI) > 0
                ) and any(
                    FI.hash
                    in [
                        x.hash
                        for x in V.dict_dataset[group][gene].list_qr_FI
                        + V.dict_dataset[group][gene].list_db_FI
                        + V.dict_dataset[group][gene].list_og_FI
                    ]
                    for FI in V.dict_dataset[group]["concatenated"].list_qr_FI
                    + V.dict_dataset[group]["concatenated"].list_db_FI
                    + V.dict_dataset[group]["concatenated"].list_og_FI
                )

                # Interpret tree when valid condition
                if cond1 or cond2 or cond3:
                    if len(V.dict_dataset[group][gene].list_og_FI) > 0:
                        # Generating tree_interpretation opts for multithreading support
                        yield (
                            (
                                f"{opt.runname}_{group}_{gene}",
                                group,
                                gene,
                                V.tup_genus,
                                hash_dict,
                                query_list,
                                outgroup,
                                partition,
                                path,
                                opt,
                            )
                        )
                    # However, if outgroup does not exists, warn it
                    else:
                        logging.warning(
                            f"Failed interpreting tree {group} {gene} because no outgroup available"
                        )

    tree_interpretation_opt = generate_interpretation_opt()

    tree_info_list = []

    ## Tree interpretation - outgroup, reconstruction(solve_flat), collapsing
    if opt.verbose < 3:
        with mp.Pool(opt.thread) as p:
            tree_info_list.extend(
                p.starmap(_safe_pipe_module_tree_interpretation, tree_interpretation_opt)
            )

    else:
        # non-multithreading mode for debugging
        tree_info_list = [
            pipe_module_tree_interpretation(*option)
            for option in tree_interpretation_opt
        ]

    # Drop genera that failed interpretation (guarded above) so one bad tree
    # does not sink the whole run; failures are logged as [TREE INTERPRETATION FAILED].
    _n_failed = sum(1 for ti in tree_info_list if ti is None)
    if _n_failed:
        logging.warning(f"{_n_failed} (group, gene) trees failed interpretation and were skipped")
    tree_info_list = [ti for ti in tree_info_list if ti is not None]

    # Gather flat branch issues
    for tree_info in tree_info_list:
        for flat_hash in tree_info.flat_clades:
            FI = V.dict_hash_FI[flat_hash]
            FI.issues.add(f"flat:{tree_info.gene}")

    synchronized_tree_info_list = synchronize(V, path, tree_info_list)
    tree_info_list = synchronized_tree_info_list

    # V.dict_hash_FI is too large for multiprocessing, extract only necessary part
    V_dict_hash_FI = {}
    for key in V.dict_hash_FI:
        FI = V.dict_hash_FI[key]
        V_dict_hash_FI[key] = (FI.original_id, FI.hash, FI.adjusted_group)

    # Generate options using generator
    def generate_visualization_opt():
        for tree_info in tree_info_list:
            yield (tree_info, V.tup_genus, V_dict_hash_FI, path, opt)

    # Generate visualization option to run
    tree_visualization_opt = generate_visualization_opt()
    ## Tree visualization
    if opt.verbose < 3:
        with mp.Pool(opt.thread) as p:
            tree_visualization_result = p.starmap(
                _safe_pipe_module_tree_visualization, tree_visualization_opt
            )

    else:
        # non-multithreading mode for debugging
        tree_visualization_result = [
            pipe_module_tree_visualization(*option) for option in tree_visualization_opt
        ]

    # Drop genera that failed visualization (guarded) before flattening
    tree_visualization_result = [r for r in tree_visualization_result if r is not None]

    ### Collect identifiation result to V for reporting
    # Merge report list
    all_report_list = [
        x for report_list in tree_visualization_result for x in report_list
    ]

    # hash_dict_analysis to prevent overwrite analysis from other tree
    # hash : group_analysis
    hash_dict_analysis = {}

    # Add final identification result
    # Concatenated first
    for singlereport in all_report_list:
        FI = V.dict_hash_FI[singlereport.hash]
        # Concatenated
        if singlereport.gene == "concatenated":
            cond = 0
            # If the strain has not been reported
            if not singlereport.hash in hash_dict_analysis:
                cond = 1
            # If the strain has reported, but not in major tree
            else:
                if hash_dict_analysis[singlereport.hash] != singlereport.group:
                    cond = 2

            if cond > 0:
                FI.final_species = singlereport.species_assigned
                FI.species_identifier = singlereport.ambiguous

                if singlereport.ambiguous > 0:
                    FI.issues.add("polyphyly:concatenated")

                if singlereport.flat is True:
                    FI.flat.append("concatenated")

                hash_dict_analysis[singlereport.hash] = singlereport.group_analysis

    # Non-concatenated
    for singlereport in all_report_list:
        FI = V.dict_hash_FI[singlereport.hash]
        if singlereport.gene != "concatenated":
            if singlereport.hash in hash_dict_analysis:
                # In each gene tree, follow the group which taxon analyzed from concatenated tree
                if hash_dict_analysis[singlereport.hash] == singlereport.group_analysis:
                    FI.bygene_species[singlereport.gene] = singlereport.species_assigned

                    if singlereport.ambiguous > 0:
                        FI.issues.add(f"polyphyly:{singlereport.gene}")

                    if singlereport.flat is True:
                        FI.flat.append(singlereport.gene)
            else:
                # If not found, FI would be additional database sequences outside the group added for ambiguities.
                pass

    return V, path, opt

# Performing multiple tree interpretation
from ete3 import Tree
from funvip.src import tree_interpretation
from funvip.src.tool import initialize_path, get_genus_species
from funvip.src.tool import sizeof_fmt
from funvip.src.hasher import encode, decode
from funvip.src.reporter import Singlereport
from copy import deepcopy
import pandas as pd
import re
import sys
import os
import shutil
import logging
import multiprocessing as mp


### For single dataset
# Input : out, group, gene, V, path, opt
def pipe_module_tree_interpretation(
    out,
    group,
    gene,
    V_tup_genus,
    funinfo_dict,
    funinfo_list,
    hash_dict,
    query_list,
    outgroup,
    partition,
    path,
    opt,
):
    # To reduce memory usage in multithreaded performance, copy necessary objects and then remove V

    """
    funinfo_dict = V.dict_hash_FI
    funinfo_list = V.list_FI
    hash_dict = V.dict_hash_name
    query_list = V.dict_dataset[group][gene].list_qr_FI
    outgroup = V.dict_dataset[group][gene].list_og_FI
    partition = V.partition[group]
    """

    # for unexpectively included sequence during clustering
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
            if opt.method.tree.lower() == "iqtree":
                Tree(tree_name, format=0)
            else:
                Tree(tree_name, format=2)
        except:
            logging.error(f"[DEVELOPMENTAL ERROR] Failed on importing tree {tree_name}")
            raise Exception
    else:
        logging.error(f"Cannot find {tree_name}")
        raise Exception

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

    # Reroot outgroup and save original tree into image
    tree_info.reroot_outgroup(
        f"{path.out_tree}/hash_{opt.runname}_{group}_{gene}_original.svg"
    )
    # Decode hash of image
    # Should work more on non-safe characters
    tree_hash_dict = encode(funinfo_list, newick=True)

    decode(
        tree_hash_dict,
        f"{path.out_tree}/hash_{opt.runname}_{group}_{gene}_original.svg",
        f"{path.out_tree}/{opt.runname}_{group}_{gene}_original.svg",
        newick=True,
    )

    # In validation mode, use original sp. number
    if opt.mode == "validation":
        tree_info.reserve_sp()

    # Reconstruct flat branches if option given
    if opt.solveflat is True:
        tree_info.t = tree_info.reconstruct(
            clade=tree_info.t.copy("newick"), gene=gene, opt=opt
        )

    # reorder tree for pretty look
    tree_info.t.ladderize(direction=1)

    # save current status into save version of tree
    tree_info.t_publish = deepcopy(tree_info.t)

    # Search tree and delimitate species
    tree_info.tree_search(tree_info.t, gene)

    # Move original newick and replace with adjusted ones
    shutil.move(
        f"{path.out_tree}/{opt.runname}_{group}_{gene}.nwk",
        f"{path.out_tree}/{opt.runname}_{group}_{gene}_original.nwk",
    )
    tree_info.t.write(
        format=0, outfile=f"{path.out_tree}/{opt.runname}_{group}_{gene}.nwk"
    )
    decode(
        tree_hash_dict,
        f"{path.out_tree}/{opt.runname}_{group}_{gene}.nwk",
        f"{path.out_tree}/{opt.runname}_{group}_{gene}.nwk",
        newick=True,
    )

    return tree_info


### synchronize sp. numbers from multiple dataset
# to use continuous sp numbers over trees
# Seperated from multithreading, because this step should traverse over multiple trees, therefore cannot be done simultaniously
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
            except:
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
            logging.error("DEVELOPMENTAL ERROR, DUPLICATED TREE_INFO")
            raise Exception

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

    # print(f"taxon_string_list | {group} {gene}:\n {taxon_string_list}\n")

    # Polish tree image
    tree_info.polish_image(
        f"{path.out_tree}/{opt.runname}_{group}_{gene}.svg",
        taxon_string_dict,
        genus_list,
    )

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

    for name, size in sorted(
        ((name, sys.getsizeof(value)) for name, value in list(locals().items())),
        key=lambda x: -x[1],
    )[:10]:
        print("{:>30}: {:>8}".format(name, sizeof_fmt(size)))

    raise Exception

    return report_list


### For all datasets, multiprocessing part
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
                                funinfo_dict,
                                funinfo_list,
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

    ## Tree interpretation - outgroup, reconstruction(solve_flat), collapsing
    if opt.verbose < 3:
        with mp.Pool(opt.thread) as p:
            tree_info_list = p.starmap(
                pipe_module_tree_interpretation, tree_interpretation_opt
            )

    else:
        # non-multithreading mode for debugging
        tree_info_list = [
            pipe_module_tree_interpretation(*option)
            for option in tree_interpretation_opt
        ]

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
                pipe_module_tree_visualization, tree_visualization_opt
            )

    else:
        # non-multithreading mode for debugging
        tree_visualization_result = [
            pipe_module_tree_visualization(*option) for option in tree_visualization_opt
        ]

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

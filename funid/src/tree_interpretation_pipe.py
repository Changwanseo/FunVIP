# Performing multiple tree interpretation
from ete3 import Tree
from funid.src import tree_interpretation
from funid.src.tool import initialize_path, get_genus_species, mkdir
from funid.src.hasher import encode, decode
from funid.src.reporter import Singlereport

import pandas as pd
import re
import sys, os
import shutil
import logging
import multiprocessing as mp


### For single dataset
# Input : out, group, gene, V, path, opt
def pipe_module_tree_interpretation(
    out,
    group,
    gene,
    V,
    path,
    opt,
):
    # To reduce memory usage in multithreaded performance, copy necessary objects and then remove V
    funinfo_dict = V.dict_hash_FI
    funinfo_list = V.list_FI
    hash_dict = V.dict_hash_name
    query_list = V.dict_dataset[group][gene].list_qr_FI
    outgroup = V.dict_dataset[group][gene].list_og_FI

    # for unexpectively included sequence during clustering
    db_list = list(
        set([FI for FI in V.list_FI if FI.datatype == "db"])
        - set(outgroup)
        - set(query_list)
    )
    genus_list = V.tup_genus

    del V

    # For get_genus_species
    initialize_path(path)

    # Tree name selection for tree construction software
    tree_name = f"{path.out_tree}/hash/hash_{opt.runname}_{group}_{gene}.nwk"
    logging.debug(tree_name)

    if os.path.isfile(tree_name):
        try:
            Tree(tree_name, format=2)
        except:
            for FI in funinfo_list:
                print(FI.original_id)
                print(FI.datatype)
            print(f"[DEVELOPMENTAL ERROR] Failed on tree {tree_name}")
            raise Exception
    else:
        logging.warning(f"Cannot find {tree_name}")
        raise Exception

    # initialize before analysis
    Tree_style = tree_interpretation.Tree_style()

    # Read tree
    tree_info = tree_interpretation.Tree_information(
        tree_name, Tree_style, group, gene, opt
    )

    # Give necessary variables parsed from dataset
    tree_info.db_list = db_list
    tree_info.query_list = query_list
    tree_info.outgroup = outgroup
    tree_info.funinfo_dict = funinfo_dict

    # Main phase
    # calculate zero distance with alignment
    tree_info.calculate_zero(
        f"{path.out_alignment}/{opt.runname}_hash_trimmed_{group}_{gene}.fasta"
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
        tree_info.t = tree_info.reconstruct(tree_info.t.copy("newick"), gene, opt)

    # reorder tree for pretty look
    tree_info.t.ladderize(direction=1)

    # Search tree and delimitate species
    tree_info.tree_search(tree_info.t, gene)

    return tree_info


### synchronize sp. numbers from multiple dataset
# to use continuous sp numbers over trees
# Seperated because this step cannot be done by multiprocessing
def synchronize(V, path, tree_info_list):
    # get available groups per genus
    tree_info_dict = {}

    # To synchronize sp. number by genus, generate by genus dataset
    for tree_info in tree_info_list:
        # Get possible genus from tree_info
        possible_genus = set([taxon[0] for taxon in tree_info.collapse_dict])
        logging.debug(f"Synchronize for genus : {possible_genus}")

        # Generating tree_info_dict - genus - group - tree_info cascade
        for genus in possible_genus:
            if not genus in tree_info_dict:
                tree_info_dict[genus] = {tree_info.group: {tree_info.gene: tree_info}}
            elif not (tree_info.group in tree_info_dict[genus]):
                tree_info_dict[genus][tree_info.group] = {tree_info.gene: tree_info}
            elif not (tree_info.gene in tree_info_dict[genus][tree_info.group]):
                tree_info_dict[genus][tree_info.group][tree_info.gene] = tree_info
            else:
                logging.error("DEVELOPMENTAL ERROR, DUPLICATED TREE_INFO")
                raise Exception

    # By gene sp. number dictionary
    # ex) {Penicillium: {Aspergillus: {(Penicillium, sp. 1):(Penicillium, sp. 1)}, Penicillium:{(Penicillium, sp. 1):(Penicillium, sp. 2)}}}
    sp_convert_dict = {}
    sp_cnt_dict = {}
    # Make one hash - set of available concatenated taxon dict pair
    # {HS0HE : {(Penicillium sp. 1), (Penicillium sp. 2)}}
    hash_dict = {}

    for genus in tree_info_dict.keys():
        # initialize gene sp. dictionary
        if not (genus in sp_convert_dict):
            sp_convert_dict[genus] = {}  # conversion pair
            sp_cnt_dict[genus] = 1  # counter

        # sp. numbers should be counted by genus
        for group in sorted(list(tree_info_dict[genus].keys())):
            sp_convert_dict[genus][group] = {}
            logging.debug(f"Synchronizing sp. number of {genus}, {group}")

            # Now concatenated analysis gets mandatory
            # However, some genera can only exist in certain gene analysis - they are not correlated to synchronizing
            if not ("concatenated" in tree_info_dict[genus][group]):
                logging.info(
                    f"No concatenated dataset for {genus}, {group}. Passing synchronizing"
                )
            else:
                # Start with concatenated to define standard sp numbers
                tree_info = tree_info_dict[genus][group]["concatenated"]

                for taxon in tree_info.collapse_dict.keys():
                    # Get sp. nov taxon
                    if taxon[0] == genus and re.fullmatch(r"sp. [0-9]+", taxon[1]):
                        logging.debug(f"Sp. checking : {taxon[1]}")
                        # Save sp. nov renaming pair to dict
                        sp_convert_dict[genus][group][taxon] = (
                            taxon[0],
                            f"sp. {sp_cnt_dict[genus]}",
                        )
                        # Increase counter
                        sp_cnt_dict[genus] += 1

                        # Add given sp. taxons, generate hash dict
                        # leaf[0] should be hash
                        for leaf in tree_info.collapse_dict[taxon][0].leaf_list:
                            if not (leaf[0]) in hash_dict:
                                hash_dict[leaf[0]] = {
                                    sp_convert_dict[genus][group][taxon]
                                }
                            else:
                                hash_dict[leaf[0]].add(
                                    sp_convert_dict[genus][group][taxon]
                                )

                    # Else, just add to hash_dict to get available taxon list for hash
                    elif taxon[0] == genus:
                        for leaf in tree_info.collapse_dict[taxon][0].leaf_list:
                            if not (leaf[0]) in hash_dict:
                                hash_dict[leaf[0]] = {taxon}
                            else:
                                hash_dict[leaf[0]].add(taxon)

                ## Update tree_info.collapse_dict with cnt_sp_adder
                # Perform in two steps in order to take collapse_dict safe
                # First, pop original name and add with tmp
                tmp_taxon_list = []

                for taxon in sp_convert_dict[genus][group].keys():
                    tmp_taxon = (
                        taxon[0],
                        f"tmp {sp_convert_dict[genus][group][taxon][1]}",
                    )
                    tree_info.collapse_dict[tmp_taxon] = tree_info.collapse_dict.pop(
                        taxon
                    )
                    tmp_taxon_list.append(tmp_taxon)

                # Then, remove tmp
                for taxon in tmp_taxon_list:
                    tree_info.collapse_dict[
                        (
                            taxon[0],
                            taxon[1][4:],
                        )
                    ] = tree_info.collapse_dict.pop(taxon)

                # Update taxon
                for taxon in tree_info.collapse_dict:
                    for n, collapse_info in enumerate(tree_info.collapse_dict[taxon]):
                        collapse_info.taxon = taxon
                        if len(tree_info.collapse_dict[taxon]) == 1:
                            collapse_info.clade_cnt = 0
                        else:
                            collapse_info.clade_cnt = n + 1

                ### CONCATENATED PART DONE
                # Now solve other genes
                for gene in tree_info_dict[genus][group]:
                    if gene != "concatenated":
                        # bygene_taxon_dict : {Penicillium citrinum : {(Penicillium, sp. 1), (Penicillium, sp. 2), (Penicillium, citrinum)}}}
                        # bygene_taxon_string_dict : {Penicillium citrinum : (Penicillium, citrinum/sp.1/sp.2)}
                        bygene_taxon_dict = {}
                        bygene_taxon_string_dict = {}

                        # Grab the gene tree
                        tree_info = tree_info_dict[genus][group][gene]

                        # Synchronize bygene taxon name to concatenated
                        # logging.debug(f"Collapse_dict: {tree_info.collapse_dict}")
                        for taxon in tree_info.collapse_dict:
                            if taxon[0] == genus:
                                # Get all hash of designated taxon leaf
                                clade = tree_info.collapse_dict[taxon][0]

                                try:
                                    pass
                                    # print(clade.leaf_list)
                                except:
                                    logging.debug(
                                        f"Collapse_dict: {tree_info.collapse_dict}"
                                    )
                                    raise Exception

                                hash_list = [leaf[0] for leaf in clade.leaf_list]

                                # Name of the collapsed clade should include all taxon names of hash
                                for _hash in hash_list:
                                    if _hash in hash_dict:
                                        # Only work with related taxon
                                        if not taxon in bygene_taxon_dict:
                                            bygene_taxon_dict[taxon] = set()

                                        bygene_taxon_dict[taxon] = bygene_taxon_dict[
                                            taxon
                                        ].union(hash_dict[_hash])

                        # Make bygene_taxon_dict to string
                        # Only work with current genus
                        for taxon in bygene_taxon_dict:
                            # get sp. and non sp. taxons
                            sp_list = [
                                _taxon
                                for _taxon in bygene_taxon_dict[taxon]
                                if _taxon[0] == genus
                                and re.fullmatch(r"sp. [0-9]+", _taxon[1])
                            ]

                            nonsp_list = [
                                _taxon
                                for _taxon in bygene_taxon_dict[taxon]
                                if _taxon[0] == genus
                                and not (re.fullmatch(r"sp. [0-9]+", _taxon[1]))
                            ]

                            # Make all taxon into one tuple
                            sp_species_epithet_list = sorted([x[1] for x in sp_list])
                            nonsp_species_epithet_list = sorted(
                                x[1] for x in nonsp_list
                            )

                            # Sort by nonsp - sp order
                            sp_string = "/".join(
                                nonsp_species_epithet_list + sp_species_epithet_list
                            )

                            # Remove unexpected slashes when species epithet is empty
                            sp_string.replace("//", "/")
                            if sp_string.startswith("/"):
                                sp_string = sp_string[1:]
                            if sp_string.endswith("/"):
                                sp_string = sp_string[:-1]

                            bygene_taxon_string_dict[taxon] = (genus, sp_string)

                        ## Update tree_info.collapse_dict with designated name
                        # Perform in two steps in order to take collapse_dict safe

                        ## Update tree_info.collapse_dict with cnt_sp_adder
                        # Perform in two steps in order to take collapse_dict safe
                        # First, pop original name and add with tmp
                        change_taxon_list = [
                            taxon
                            for taxon in tree_info.collapse_dict.keys()
                            if (
                                taxon[0] == genus
                                and taxon in tree_info.collapse_dict
                                and taxon in bygene_taxon_string_dict
                            )
                        ]
                        tmp_taxon_list = []

                        for taxon in change_taxon_list:
                            tmp_taxon = (
                                taxon[0],
                                f"tmp {bygene_taxon_string_dict[taxon][1]}",
                            )
                            if not (tmp_taxon in tree_info.collapse_dict):
                                tree_info.collapse_dict[
                                    tmp_taxon
                                ] = tree_info.collapse_dict.pop(taxon)
                                tmp_taxon_list.append(tmp_taxon)
                            else:
                                tree_info.collapse_dict[
                                    tmp_taxon
                                ] += tree_info.collapse_dict[taxon]
                                tree_info.collapse_dict.pop(taxon)

                        for taxon in tmp_taxon_list:
                            tree_info.collapse_dict[
                                (taxon[0], taxon[1][4:])
                            ] = tree_info.collapse_dict.pop(taxon)

                        # Update taxon
                        for taxon in tree_info.collapse_dict:
                            for n, collapse_info in enumerate(
                                tree_info.collapse_dict[taxon]
                            ):
                                collapse_info.taxon = taxon
                                if len(tree_info.collapse_dict[taxon]) == 1:
                                    collapse_info.clade_cnt = 0
                                else:
                                    collapse_info.clade_cnt = n + 1

    # Return sp number fixed tree_info_list
    return tree_info_list


### Visualization after synchronization
def pipe_module_tree_visualization(
    tree_info,
    V,
    path,
    opt,
):
    group = tree_info.group
    gene = tree_info.gene
    genus_list = V.tup_genus

    # Collapse tree branches for visualization
    taxon_string_list = tree_info.collapse_tree()

    # Polish tree image
    tree_info.polish_image(
        f"{path.out_tree}/{opt.runname}_{group}_{gene}.svg",
        taxon_string_list,
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
                report.id = V.dict_hash_FI[leaf[0]].original_id
                report.hash = V.dict_hash_FI[leaf[0]].hash
                report.update_group(group)
                report.update_gene(gene)
                report.update_species_original(
                    get_genus_species(leaf[2], genus_list=genus_list)
                )
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
                    report.id = V.dict_hash_FI[leaf[0]].original_id
                    report.hash = V.dict_hash_FI[leaf[0]].hash
                    report.update_group(group)
                    report.update_gene(gene)
                    report.update_species_original(
                        get_genus_species(leaf[2], genus_list=genus_list)
                    )
                    report.update_species_assigned((f"{taxon[0]} {taxon[1]} {n+1}"))

                    report.ambiguous = collapse_info.clade_cnt
                    report.flat = collapse_info.flat

                    report_list.append(report)

    return report_list


### For all datasets, multiprocessing part
def pipe_tree_interpretation(V, path, opt):
    # Generate tree_interpretation opt to run
    tree_interpretation_opt = []
    for group in V.dict_dataset:
        for gene in V.dict_dataset[group]:
            # draw tree only when query sequence exists
            if (
                len(V.dict_dataset[group][gene].list_qr_FI) > 0
                or opt.queryonly is False
            ) and len(V.dict_dataset[group][gene].list_og_FI) > 0:
                # Generating tree_interpretation opts for multithreading support
                tree_interpretation_opt.append(
                    (
                        f"{opt.runname}_{group}_{gene}",
                        group,
                        gene,
                        V,
                        path,
                        opt,
                    )
                )

    ## Tree interpretation - outgroup, reconstruction(solve_flat), collapsing
    if opt.verbose < 3:
        p = mp.Pool(opt.thread)
        tree_info_list = p.starmap(
            pipe_module_tree_interpretation, tree_interpretation_opt
        )
        p.close()
        p.join()

    else:
        # non-multithreading mode for debugging
        tree_info_list = [
            pipe_module_tree_interpretation(*option)
            for option in tree_interpretation_opt
        ]

    tree_info_list = synchronize(V, path, tree_info_list)
    # Generate visualization option to run
    tree_visualization_opt = []
    for tree_info in tree_info_list:
        tree_visualization_opt.append((tree_info, V, path, opt))

    ## Tree visualization
    if opt.verbose < 3:
        p = mp.Pool(opt.thread)
        tree_visualization_result = p.starmap(
            pipe_module_tree_visualization, tree_visualization_opt
        )
        p.close()
        p.join()

    else:
        # non-multithreading mode for debugging
        tree_visualization_result = [
            pipe_module_tree_visualization(*option) for option in tree_visualization_opt
        ]

    ### After this

    # Collect identifiation result to V
    for report_list in tree_visualization_result:
        for report in report_list:
            FI = V.dict_hash_FI[report.hash]
            # Concatenated
            if report.gene == "concatenated":
                FI.final_species = report.species_assigned
                FI.species_identifier = report.ambiguous
                if report.flat is True:
                    FI.flat.append("concatenated")
            # Non concatenated
            else:
                FI.bygene_species[report.gene] = report.species_assigned
                if report.flat is True:
                    FI.flat.append(report.gene)

    ### Before this

    return V, path, opt

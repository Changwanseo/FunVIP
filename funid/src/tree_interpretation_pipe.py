# Performing multiple tree interpretation
from ete3 import Tree
from funid.src import tree_interpretation
from funid.src.tool import initialize_path, get_genus_species
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

    # Check validity of file while importing
    if os.path.isfile(tree_name):
        try:
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
# Seperated because this step should traverse over multiple trees, therefore cannot be done simultaniously
def synchronize(V, path, tree_info_list):
    ## Initialize

    # get available groups per genus
    tree_info_dict = {}
    # hash : corresponding taxon
    hash_taxon_dict = {}
    # genus : cnt
    sp_cnt_dict = {}

    # To synchronize sp. number by genus, generate by-genus dataset
    for tree_info in tree_info_list:
        if not (tree_info.group in tree_info_dict):
            tree_info_dict[tree_info.group] = {tree_info.gene: tree_info}
        elif not (tree_info.gene in tree_info_dict[tree_info.group]):
            tree_info_dict[tree_info.group][tree_info.gene] = tree_info
        else:
            logging.error("DEVELOPMENTAL ERROR, DUPLICATED TREE_INFO")
            raise Exception

    # Memoize iterative calling
    valid_hash_dict = {}
    for group in tree_info_dict:
        valid_hash_dict[group] = [
            _hash
            for _hash in V.dict_hash_FI
            if V.dict_hash_FI[_hash].datatype == "db"
            and V.dict_hash_FI[_hash].group == group
        ]

    query_hash_list = [
        _hash for _hash in V.dict_hash_FI if V.dict_hash_FI[_hash].datatype == "query"
    ]

    ## Starting with concatenated
    # In priority, count corresponding group taxa first
    for group in tree_info_dict:
        if "concatenated" in tree_info_dict[group]:
            tree_info = tree_info_dict[group]["concatenated"]
            # Get list of hash in interest
            valid_hash_list = valid_hash_dict[group]
            for taxon in tree_info.collapse_dict:
                clade_list = tree_info.collapse_dict[taxon]
                for n, clade in enumerate(clade_list):
                    hash_list = [leaf[0] for leaf in clade.leaf_list]
                    # If any of the leaf consisting clade is included to valid_hash_list
                    if any(_h in valid_hash_list for _h in hash_list):
                        for _hash in hash_list:
                            if not (_hash) in hash_taxon_dict:
                                hash_taxon_dict[_hash] = (taxon, n + 1)
                            else:
                                print(
                                    f"{_hash} collided while putting in hash_taxon_dict"
                                )

                    # If leaf consisting with only query, that won't collide with other groups
                    if all(_h in query_hash_list for _h in hash_list):
                        if not (taxon[0] in sp_cnt_dict):
                            sp_cnt_dict[taxon[0]] = 1
                        for _hash in hash_list:
                            if not (_hash) in hash_taxon_dict:
                                hash_taxon_dict[_hash] = (
                                    (taxon[0], "sp."),
                                    sp_cnt_dict[taxon[0]],
                                )
                            else:
                                print(
                                    f"{_hash} collided while putting in hash_taxon_dict"
                                )

                        sp_cnt_dict[taxon[0]] += 1

    # Next, taxa that doesn't belongs to any of the group
    for group in tree_info_dict:
        if "concatenated" in tree_info_dict[group]:
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
                            else:
                                print(
                                    f"{_hash} collided while putting in hash_taxon_dict. Tried to put {(taxon, n+1)}, but existing {hash_taxon_dict[_hash]}"
                                )

                                # print(group, _hash, taxon, n + 1)

    # DEBUGGING
    """
    for _hash in sorted(list(hash_taxon_dict.keys())):
        print(_hash, hash_taxon_dict[_hash])

    raise Exception
    """

    # Generate final taxon name for synchronizing
    def get_new_taxon(hash_list, hash_taxon_dict):
        taxon_candidates = set()
        for _hash in hash_list:
            taxon_candidates.add(hash_taxon_dict[_hash])

        list_taxon_candidates = sorted(list(taxon_candidates))
        genus = "/".join(t[0][0] for t in list_taxon_candidates)
        species_list = []
        for t in list_taxon_candidates:
            # If the taxon is unique
            if not (t[0], 2) in hash_taxon_dict.values():
                species_list.append(t[0][1])
            else:
                species_list.append(f"{t[0][1]} {t[1]}")

        species = "/".join([s for s in species_list])

        return (genus, species)

    # Now update from concatenated
    # Remove original taxon, and add by clade taxon
    for group in tree_info_dict:
        for gene in tree_info_dict[group]:
            tree_info = tree_info_dict[group][gene]
            # Before taxon, after taxon update list
            remove_list = []  # [taxon1, taxon2, taxon3 ...]
            add_list = {}  # [taxon1 : [clade1], taxon2 : [clade2] ...]
            for taxon in tree_info.collapse_dict:
                clade_list = tree_info.collapse_dict[taxon]
                remove_list.append(taxon)
                for clade in clade_list:
                    hash_list = [leaf[0] for leaf in clade.leaf_list]
                    new_taxon = get_new_taxon(hash_list, hash_taxon_dict)
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
    # By gene sp. number dictionary
    # ex) {Penicillium: {Aspergillus: {(Penicillium, sp. 1):(Penicillium, sp. 1)}, Penicillium:{(Penicillium, sp. 1):(Penicillium, sp. 2)}}}
    sp_convert_dict = {}
    sp_cnt_dict = {}
    # Make one hash - set of available concatenated taxon dict pair
    # {HS0HE : {(Penicillium sp. 1), (Penicillium sp. 2)}}
    hash_dict = {}

    # for genus in tree_info_dict:

    for genus in tree_info_dict.keys():
        # initialize gene sp. dictionary
        if not (genus in sp_convert_dict):
            sp_convert_dict[genus] = {}  # conversion pair
            sp_cnt_dict[genus] = 1  # sp. nov. number counter

        # sp. numbers should be counted by genus
        for group in sorted(list(tree_info_dict[genus].keys())):
            # Number counts for each taxa are generated, according to their main group
            if genus in V.dict_group[group]:
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
                        ## sp. nov taxa
                        # Get sp. nov taxon
                        if taxon[0] == genus and re.fullmatch(r"sp. [0-9]+", taxon[1]):
                            logging.debug(f"sp. checking : {taxon[1]}")
                            # Save sp. nov renaming pair to dict
                            sp_convert_dict[genus][group][taxon] = (
                                taxon[0],
                                f"sp. {sp_cnt_dict[genus]}",
                            )
                            # Increase sp. nov number counter
                            sp_cnt_dict[genus] += 1

                            # Add given sp. taxons, generate {hash: taxon} dict
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
                        ## non sp. nov taxa
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
                        tree_info.collapse_dict[
                            tmp_taxon
                        ] = tree_info.collapse_dict.pop(taxon)
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
                        for n, collapse_info in enumerate(
                            tree_info.collapse_dict[taxon]
                        ):
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

                                            bygene_taxon_dict[
                                                taxon
                                            ] = bygene_taxon_dict[taxon].union(
                                                hash_dict[_hash]
                                            )

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
                                sp_species_epithet_list = sorted(
                                    [x[1] for x in sp_list]
                                )
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

    # Working for unintendly included groups
    for genus in tree_info_dict.keys():
        for group in sorted(list(tree_info_dict[genus].keys())):
            # In this time, concatenated does not matter
            for gene in tree_info_dict[genus][group]:
                changes = []
                tree_info = tree_info_dict[genus][group][gene]
                for taxon in tree_info.collapse_dict.keys():
                    if taxon[0] != genus:
                        print(
                            f"Line 392 debugging! genus {genus} group {group} gene {gene} taxon {taxon}"
                        )
                        for clade in tree_info.collapse_dict[taxon]:
                            hash_list = [leaf[0] for leaf in clade.leaf_list]
                            # Name of the collapsed clade should include all taxon names of hash
                            hash_taxon_set = set()
                            for _hash in hash_list:
                                if _hash in hash_dict:
                                    print(hash_dict[_hash])
                                    hash_taxon_set = hash_taxon_set.union(
                                        hash_dict[_hash]
                                    )

                        # Change the name of taxon as designated
                        hash_taxon_list = sorted(list(hash_taxon_set))

                        tmp_taxon = (
                            "/".join([t[0] for t in hash_taxon_list]),
                            "/".join([t[1] for t in hash_taxon_list]),
                        )

                        changes.append((tmp_taxon, taxon))
                for change in changes:
                    tree_info.collapse_dict[change[0]] = tree_info.collapse_dict.pop(
                        change[1]
                    )

        pass
    """

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

    print(f"taxon_string_list | {group} {gene}:\n {taxon_string_list}\n")

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
    print("-----------------------------------------------")
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
                print(f"Conc {FI.hash} {FI.id} {FI.datatype} {FI.final_species}")
            # Non concatenated
            else:
                FI.bygene_species[report.gene] = report.species_assigned
                if report.flat is True:
                    FI.flat.append(report.gene)
                print(f"Else {FI.hash} {FI.id} {FI.datatype} {FI.final_species}")

    # raise Exception

    ### Before this

    return V, path, opt

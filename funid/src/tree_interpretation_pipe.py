from ete3 import Tree
from funid.src import tree_interpretation
from funid.src.tool import initialize_path, get_genus_species, mkdir
from funid.src.hasher import encode, decode
from funid.src.reporter import Singlereport
import pandas as pd
import sys, os
import shutil
import logging
import multiprocessing as mp

# Single dataset
def pipe_module_tree_interpretation(
    out,
    sect,
    gene,
    V,
    path,
    opt,
):

    # to reduce memory usage in multithreaded performance, copy necessary objects and then remove V
    funinfo_dict = V.dict_hash_FI
    funinfo_list = V.list_FI
    hash_dict = V.dict_hash_name
    db_list = V.dict_dataset[sect][gene].list_db_FI
    query_list = V.dict_dataset[sect][gene].list_qr_FI
    outgroup = V.dict_dataset[sect][gene].list_og_FI
    genus_list = V.tup_genus

    del V

    # for get_genus_species
    initialize_path(path)

    # Tree name selection for tree construction software
    try:
        tree_name = f"{path.out_tree}/hash/hash_{out}.nwk"
        if os.path.isfile(tree_name):
            pass
        else:
            print(f"Cannot find {tree_name}")
        Tree(tree_name)

    except:
        logging.warning(f"Cannot read {tree_name}")
        return None

    try:
        # initialize before analysis
        Tree_style = tree_interpretation.Tree_style()

        # Read tree
        tree_info = tree_interpretation.Tree_information(tree_name, Tree_style, opt)

        # give necessary variables parsed from dataset
        tree_info.db_list = db_list
        tree_info.query_list = query_list
        tree_info.outgroup = outgroup
        tree_info.funinfo_dict = funinfo_dict

        # make zero with alignment
        if gene != "concatenated":
            tree_info.calculate_zero(
                f"{path.out_alignment}/{opt.runname}_hash_trimmed_{sect}_{gene}.fasta"
            )
        else:
            tree_info.calculate_zero(
                f"{path.out_alignment}/{opt.runname}_hash_trimmed_{sect}_concatenated.fasta"
            )

        # Reroot outgroup and save original tree into image
        tree_info.reroot_outgroup(f"{path.out_tree}/hash_{out}_original.svg")
        tree_hash_dict = encode(funinfo_list, newick=True)
        decode(
            tree_hash_dict,
            f"{path.out_tree}/hash_{out}_original.svg",
            f"{path.out_tree}/{out}_original.svg",
            newick=True,
        )

        if opt.mode == "validation":
            tree_info.reserve_sp()

        if tree_info.opt.solveflat is True:
            tree_info.t = tree_info.reconstruct(tree_info.t.copy("newick"), gene, opt)

        tree_info.t.ladderize(direction=1)  # reorder tree
        tree_info.tree_search(tree_info.t, gene, opt=opt)
        # synchronize() # to use continuous sp numbers over trees
        tree_info.collapse_tree()
        tree_info.polish_image(f"{path.out_tree}/{out}.svg", genus_list)

        tmp_dict = {}

        # sort taxon order
        list_taxon_1 = [
            taxon
            for taxon in tree_info.collapse_dict.keys()
            if not (taxon[1].startswith("sp."))
        ]
        list_taxon_2 = [
            taxon
            for taxon in tree_info.collapse_dict.keys()
            if taxon[1].startswith("sp.")
        ]
        list_taxon_1.sort(key=lambda x: x[1])
        list_taxon_2.sort(key=lambda x: x[1])
        list_taxon = list_taxon_1 + list_taxon_2

        report_list = []
        for taxon in list_taxon:
            if len(tree_info.collapse_dict[taxon]) <= 1:
                collapse_info = tree_info.collapse_dict[taxon][0]
                tmp_dict[" ".join(taxon)] = [
                    collapse_info.n_db,
                    collapse_info.n_query,
                    collapse_info.n_others,
                    collapse_info.n_db + collapse_info.n_query + collapse_info.n_others,
                ]

                # leaf structure:
                for leaf in collapse_info.leaf_list:
                    if 1:  # tree_info.option.highlight:
                        report = Singlereport()
                        report.accession = funinfo_dict[leaf[0]].original_accession
                        report.hash = funinfo_dict[leaf[0]].hash
                        report.update_genussection(out, gene)
                        report.update_inputtaxon(
                            get_genus_species(leaf[2], genus_list=genus_list)
                        )
                        report.taxon_cnt = collapse_info.clade_cnt
                        report.update_identifiedtaxon(taxon)
                        report_list.append(report)

            else:
                for n, collapse_info in enumerate(tree_info.collapse_dict[taxon]):
                    tmp_dict[f'{" ".join(taxon)} {n+1}'] = [
                        collapse_info.n_db,
                        collapse_info.n_query,
                        collapse_info.n_others,
                        collapse_info.n_db
                        + collapse_info.n_query
                        + collapse_info.n_others,
                    ]

                    for leaf in collapse_info.leaf_list:
                        if 1:  # tree_info.option.highlight:
                            report = Singlereport()
                            report.accession = funinfo_dict[leaf[0]].original_accession
                            report.hash = funinfo_dict[leaf[0]].hash
                            report.update_genussection(out, gene)
                            report.update_inputtaxon(
                                get_genus_species(leaf[2], genus_list=genus_list)
                            )
                            report.taxon_cnt = collapse_info.clade_cnt
                            report.update_identifiedtaxon((" ".join(taxon), f"{n+1}"))
                            report_list.append(report)

        df = pd.DataFrame(tmp_dict, index=["db", "query", "others", "total"])
        df = df.transpose()

        return (out, df, report_list)

    except:  # for debugging

        if opt.verbose >= 3:
            logging.error(f"Error occured on {tree_name}, running debugging mode")
            # initialize before analysis
            Tree_style = tree_interpretation.Tree_style()
            tree_info = tree_interpretation.Tree_information(tree_name, Tree_style, opt)
            tree_info.db_list = db_list
            tree_info.query_list = query_list
            tree_info.outgroup = outgroup
            tree_info.funinfo_dict = funinfo_dict

            # make zero with alignment
            # make zero with alignment
            if gene != "concatenated":
                tree_info.calculate_zero(
                    f"{path.out_alignment}/{opt.runname}_hash_trimmed_{sect}_{gene}.fasta"
                )
            else:
                tree_info.calculate_zero(
                    f"{path.out_alignment}/{opt.runname}_hash_trimmed_{sect}_concatenated.fasta"
                )

            tree_info.reroot_outgroup(f"{path.out_tree}/hash_{out}_original.svg")
            tree_hash_dict = encode(funinfo_list, newick=True)
            decode(
                tree_hash_dict,
                f"{path.out_tree}/hash_{out}_original.svg",
                f"{path.out_tree}/{out}_original.svg",
                newick=True,
            )

            if opt.mode == "validation":
                tree_info.reserve_sp()

            if tree_info.opt.solveflat is True:
                tree_info.t = tree_info.reconstruct(
                    tree_info.t.copy("newick"), gene, opt
                )

            tree_info.t.ladderize(direction=1)  # reorder tree
            tree_info.tree_search(tree_info.t, gene, opt=opt)
            # synchronize() # to use continuous sp numbers over trees
            tree_info.collapse_tree()
            tree_info.polish_image(f"{path.out_tree}/{out}.svg", genus_list)

            tmp_dict = {}

            # sort taxon order
            list_taxon_1 = [
                taxon
                for taxon in tree_info.collapse_dict.keys()
                if not (taxon[1].startswith("sp."))
            ]
            list_taxon_2 = [
                taxon
                for taxon in tree_info.collapse_dict.keys()
                if taxon[1].startswith("sp.")
            ]
            list_taxon_1.sort(key=lambda x: x[1])
            list_taxon_2.sort(key=lambda x: x[1])
            list_taxon = list_taxon_1 + list_taxon_2

            report_list = []
            for taxon in list_taxon:
                if len(tree_info.collapse_dict[taxon]) <= 1:
                    collapse_info = tree_info.collapse_dict[taxon][0]
                    tmp_dict[" ".join(taxon)] = [
                        collapse_info.n_db,
                        collapse_info.n_query,
                        collapse_info.n_others,
                        collapse_info.n_db
                        + collapse_info.n_query
                        + collapse_info.n_others,
                    ]

                    # leaf structure:
                    for leaf in collapse_info.leaf_list:
                        if 1:  # tree_info.option.highlight:
                            report = Singlereport()
                            report.accession = funinfo_dict[leaf[0]].original_accession
                            report.hash = funinfo_dict[leaf[0]].hash
                            report.update_genussection(out, gene)
                            report.update_inputtaxon(
                                get_genus_species(leaf[2], genus_list=genus_list)
                            )
                            report.taxon_cnt = collapse_info.clade_cnt
                            report.update_identifiedtaxon(taxon)
                            report_list.append(report)

                else:
                    for n, collapse_info in enumerate(tree_info.collapse_dict[taxon]):
                        tmp_dict[f'{" ".join(taxon)} {n+1}'] = [
                            collapse_info.n_db,
                            collapse_info.n_query,
                            collapse_info.n_others,
                            collapse_info.n_db
                            + collapse_info.n_query
                            + collapse_info.n_others,
                        ]

                        for leaf in collapse_info.leaf_list:
                            if 1:  # tree_info.option.highlight:
                                report = Singlereport()
                                report.accession = funinfo_dict[
                                    leaf[0]
                                ].original_accession
                                report.hash = funinfo_dict[leaf[0]].hash
                                report.update_genussection(out, gene)
                                report.update_inputtaxon(
                                    get_genus_species(leaf[2], genus_list=genus_list)
                                )
                                report.taxon_cnt = collapse_info.clade_cnt
                                report.update_identifiedtaxon(
                                    (" ".join(taxon), f"{n+1}")
                                )
                                report_list.append(report)

            df = pd.DataFrame(tmp_dict, index=["db", "query", "others", "total"])
            df = df.transpose()

            return (out, df, report_list)

        logging.warning(f"Failed visualizing tree on {out}, passing")


# for all datasets, multiprocessing part
def pipe_tree_interpretation(V, path, opt):
    tree_interpretation_opt = []

    for sect in V.dict_dataset:
        for gene in V.dict_dataset[sect]:
            # draw tree only when query sequence exists
            if len(V.dict_dataset[sect][gene].list_qr_FI) > 0 or opt.queryonly is False:

                # Generating tree_interpretation opts for multithreading support
                tree_interpretation_opt.append(
                    (
                        f"{opt.runname}_{sect}_{gene}",
                        sect,
                        gene,
                        V,
                        path,
                        opt,
                    )
                )

    if opt.verbose < 3:
        p = mp.Pool(opt.thread)
        tree_interpretation_result = p.starmap(
            pipe_module_tree_interpretation, tree_interpretation_opt
        )
        p.close()
        p.join()

    else:
        # non-multithreading mode for debugging
        tree_interpretation_result = []
        for option in tree_interpretation_opt:
            tree_interpretation_result.append(pipe_module_tree_interpretation(*option))

    # put identified result to V
    for raw_result in tree_interpretation_result:
        if raw_result is not (None):
            for result in raw_result[2]:
                FI = V.dict_hash_FI[result.hash]
                if result.gene == "concatenated" or opt.concatenate is False:
                    # print("Updating final species")
                    FI.final_species = result.identifiedtaxon
                    FI.species_identifier = result.taxon_cnt
                else:
                    # print("Not updating final species")
                    FI.bygene_species[result.gene] = result.identifiedtaxon

    return V, path, opt

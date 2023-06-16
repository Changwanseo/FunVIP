from funid.src import ext, hasher
from funid.src.opt_generator import opt_generator
from copy import deepcopy
import multiprocessing as mp
import logging
import shutil
import os


def pipe_tree(V, path, opt, model_dict):
    # for tree, use hash dict with genus and species information
    tree_hash_dict = hasher.encode(V.list_FI, newick=True)

    # remove tree files already exists to prevent error
    for file in [f for f in os.listdir(path.sys_path) if f.endswith(".nwk")]:
        os.remove(f"{path.sys_path}/{file}")

    fasttree_opt = []  # for multiprocessing on fasttree

    tree_dataset = deepcopy(V.dict_dataset)

    # Before drawing tree, finalize datasets
    remove_dataset = []
    for group in tree_dataset:
        for gene in tree_dataset[group]:
            # draw tree only when query sequence exists if queryonly is True
            # draw tree only when outgroup sequence exists
            if len(tree_dataset[group][gene].list_qr_FI) == 0 and opt.queryonly is True:
                logging.warning(
                    f"Passing tree construction of {group} {gene} dataset because no query included"
                )
                remove_dataset.append((group, gene))
                # tree_dataset[group].pop(gene, None)
            elif tree_dataset[group][gene].list_og_FI == 0:
                logging.warning(
                    f"Passing tree construction of {group} {gene} dataset because no outgroup available"
                )
                # tree_dataset[group].pop(gene, None)
                remove_dataset.append((group, gene))

    for group, gene in remove_dataset:
        tree_dataset[group].pop(gene, None)

    # To indicate single genes for each group
    singlegene_dict = {}

    # Draw phylogenetic trees for each dataset
    for group in tree_dataset:
        singlegene_flag = 0
        singlegene = ""
        # For single gene, we don't have to draw tree multiple times
        if len(tree_dataset[group].keys()) == 2:
            singlegene_flag = 1
            singlegene = list(set(tree_dataset[group].keys()) - {"concatenated"})[0]
            singlegene_dict[group] = deepcopy(singlegene)

        for gene in tree_dataset[group]:
            # For single gene, skip concatenated
            if gene == "concatenated" and singlegene_flag == 1:
                logging.info(
                    f"Passing tree construction of {group} {gene} because single gene detected. Will copy from concatenated"
                )
            # Else for each gene
            elif gene != "concatenated":
                # If not trimming use this
                if opt.method.tree.lower() == "raxml":
                    ext.RAxML(
                        fasta=f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta",
                        out=f"{opt.runname}_{group}_{gene}.nwk",
                        hash_dict=tree_hash_dict,
                        path=path,
                        thread=opt.thread,
                        bootstrap=opt.bootstrap,
                        model=model_dict[group][gene],
                    )

                elif opt.method.tree.lower() == "iqtree":
                    ext.IQTREE(
                        fasta=f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta",
                        out=f"{opt.runname}_{group}_{gene}.nwk",
                        hash_dict=tree_hash_dict,
                        path=path,
                        memory=opt.memory,
                        thread=opt.thread,
                        bootstrap=opt.bootstrap,
                        model=model_dict[group][gene],
                    )

                # if fasttree, append to opt
                else:
                    if not (opt.method.tree.lower() == "fasttree"):
                        logging.warning(
                            "Tree construction method not selected, working for default opt, FastTree"
                        )
                    fasttree_opt.append(
                        (
                            f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta",
                            f"{opt.runname}_{group}_{gene}.nwk",
                            tree_hash_dict,
                            path,
                            model_dict[group][gene],
                        )
                    )
            # For concatenated datasets, use partition files
            else:
                if opt.method.tree.lower() == "raxml":
                    ext.RAxML(
                        fasta=f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta",
                        out=f"{opt.runname}_{group}_{gene}.nwk",
                        hash_dict=tree_hash_dict,
                        path=path,
                        thread=opt.thread,
                        bootstrap=opt.bootstrap,
                        partition=f"{path.out_alignment}/{opt.runname}_{group}.partition",
                        model=model_dict[group][gene],
                    )

                elif opt.method.tree.lower() == "iqtree":
                    ext.IQTREE(
                        fasta=f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta",
                        out=f"{opt.runname}_{group}_{gene}.nwk",
                        hash_dict=tree_hash_dict,
                        path=path,
                        memory=opt.memory,
                        thread=opt.thread,
                        bootstrap=opt.bootstrap,
                        partition=f"{path.out_alignment}/{opt.runname}_{group}.partition",
                        model=model_dict[group][gene],
                    )

                else:
                    if not (opt.method.tree.lower() == "fasttree"):
                        logging.warning(
                            "Tree method not selected, working for default opt, FastTree"
                        )
                    fasttree_opt.append(
                        (
                            f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta",
                            f"{opt.runname}_{group}_{gene}.nwk",
                            tree_hash_dict,
                            path,
                            model_dict[group][gene],
                        )
                    )

    # for fasttree, perform multiprocessing
    if opt.method.tree.lower() == "fasttree":
        # run multiprocessing start
        if opt.verbose < 3:
            p = mp.Pool(opt.thread)
            fasttree_result = p.starmap(ext.FastTree, fasttree_opt)
            p.close()
            p.join()

        else:
            # non-multithreading mode for debugging
            fasttree_result = []
            for option in fasttree_opt:
                fasttree_result.append(ext.FastTree(*option))

    for group in singlegene_dict:
        # copy concatenated result to single gene
        shutil.copy(
            f"{path.out_tree}/{opt.runname}_{group}_{singlegene_dict[group]}.nwk",
            f"{path.out_tree}/{opt.runname}_{group}_concatenated.nwk",
        )

        shutil.copy(
            f"{path.out_tree}/hash_{opt.runname}_{group}_{singlegene_dict[group]}.nwk",
            f"{path.out_tree}/hash_{opt.runname}_{group}_concatenated.nwk",
        )

    # decode alignments for each gene after building tree
    for group in tree_dataset:
        for gene in tree_dataset[group]:
            # draw tree only when query sequence exists
            if len(tree_dataset[group][gene].list_qr_FI) > 0 or opt.queryonly is False:
                # Concatenated MAFFT file does not exists
                if gene != "concatenated":
                    os.rename(
                        f"{path.out_alignment}/{opt.runname}_MAFFT_{group}_{gene}.fasta",
                        f"{path.out_alignment}/{opt.runname}_hash_MAFFT_{group}_{gene}.fasta",
                    )
                    hasher.decode(
                        tree_hash_dict,
                        f"{path.out_alignment}/{opt.runname}_hash_MAFFT_{group}_{gene}.fasta",
                        f"{path.out_alignment}/{opt.runname}_MAFFT_{group}_{gene}.fasta",
                    )

                os.rename(
                    f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta",
                    f"{path.out_alignment}/{opt.runname}_hash_trimmed_{group}_{gene}.fasta",
                )

                hasher.decode(
                    tree_hash_dict,
                    f"{path.out_alignment}/{opt.runname}_hash_trimmed_{group}_{gene}.fasta",
                    f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta",
                )

    return V, path, opt

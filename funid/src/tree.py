from funid.src import ext, hasher
from funid.src.opt_generator import opt_generator
import multiprocessing as mp
import shutil
import os


def pipe_tree(V, path, opt, model_dict):

    # for tree, use hash dict with genus and species information
    tree_hash_dict = hasher.encode(V.list_FI, newick=True)

    # remove tree files already exists to prevent error
    for file in [f for f in os.listdir(path.sys_path) if f.endswith(".nwk")]:
        os.remove(f"{path.sys_path}/{file}")

    # Single Gene
    fasttree_opt = []  # for multiprocessing on fasttree
    for group in V.dict_dataset:
        for gene in V.dict_dataset[group]:
            # draw tree only when query sequence exists
            if (
                len(V.dict_dataset[group][gene].list_qr_FI) > 0
                or opt.queryonly is False
            ):
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
                        thread=opt.thread,
                        bootstrap=opt.bootstrap,
                        model=model_dict[group][gene],
                    )

                # if fasttree, append to opt
                else:
                    if not (opt.method.tree.lower() == "fasttree"):
                        logging.warning(
                            "[Warning] Tree construction method not selected, working for default opt, FastTree"
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

    # MultiGene tree
    fasttree_concatenated_opt = []
    if opt.concatenate is True:
        for group in V.multigene_list:
            if opt.method.tree.lower() == "raxml":
                ext.RAxML(
                    f"{path.out_alignment}/{opt.runname}_{group}_concatenated.fasta",
                    f"{opt.runname}_{group}_concatenated.nwk",
                    tree_hash_dict,
                    path,
                    thread=opt.thread,
                    bootstrap=opt.bootstrap,
                    partition=f"{path.out_alignment}/{opt.runname}_{group}.partition",
                    model=model_dict[group][gene],
                )

            elif opt.method.tree.lower() == "iqtree":
                ext.IQTREE(
                    f"{path.out_alignment}/{opt.runname}_{group}_concatenated.fasta",
                    f"{opt.runname}_{group}_concatenated.nwk",
                    tree_hash_dict,
                    path,
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
                fasttree_concatenated_opt.append(
                    (
                        f"{path.out_alignment}/{opt.runname}_{group}_concatenated.fasta",
                        f"{opt.runname}_{group}_concatenated.nwk",
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
            fasttree_concatenated_result = p.starmap(
                ext.FastTree, fasttree_concatenated_opt
            )
            p.close()
            p.join()

        else:
            # non-multithreading mode for debugging
            fasttree_concatenated_result = []
            for option in fasttree_concatenated_opt:
                fasttree_concatenated_result.append(ext.FastTree(*option))

    # decode alignments for each gene after building tree
    for group in V.dict_dataset:
        for gene in V.dict_dataset[group]:
            # draw tree only when query sequence exists
            if (
                len(V.dict_dataset[group][gene].list_qr_FI) > 0
                or opt.queryonly is False
            ):
                os.rename(
                    f"{path.out_alignment}/{opt.runname}_MAFFT_{group}_{gene}.fasta",
                    f"{path.out_alignment}/{opt.runname}_hash_MAFFT_{group}_{gene}.fasta",
                )

                os.rename(
                    f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta",
                    f"{path.out_alignment}/{opt.runname}_hash_trimmed_{group}_{gene}.fasta",
                )

                # decoding
                hasher.decode(
                    tree_hash_dict,
                    f"{path.out_alignment}/{opt.runname}_hash_MAFFT_{group}_{gene}.fasta",
                    f"{path.out_alignment}/{opt.runname}_MAFFT_{group}_{gene}.fasta",
                )

                hasher.decode(
                    tree_hash_dict,
                    f"{path.out_alignment}/{opt.runname}_hash_trimmed_{group}_{gene}.fasta",
                    f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta",
                )

    # decode concatenated fasta after building tree
    if opt.concatenate is True:
        for group in V.dict_dataset:
            try:
                os.rename(
                    f"{path.out_alignment}/{opt.runname}_{group}_concatenated.fasta",
                    f"{path.out_alignment}/{opt.runname}_hash_{group}_concatenated.fasta",
                )

                hasher.decode(
                    tree_hash_dict,
                    f"{path.out_alignment}/{opt.runname}_hash_{group}_concatenated.fasta",
                    f"{path.out_alignment}/{opt.runname}_{group}_concatenated.fasta",
                )
            except:
                pass

    return V, path, opt

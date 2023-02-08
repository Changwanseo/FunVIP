def main():

    from funid.src import align
    from funid.src import tree_interpretation_pipe
    from funid.src import cluster
    from funid.src import dataset
    from funid.src import ext
    from funid.src import hasher
    from funid.src import initialize
    from funid.src import logger
    from funid.src import manage_input
    from funid.src import modeltest
    from funid.src import multigene
    from funid.src import patch
    from funid.src import reporter
    from funid.src import tool
    from funid.src import save
    from funid.src import search
    from funid.src import tree
    from funid.src import trim
    from funid.src.command import CommandParser
    from funid.src.tool import index_step
    from funid.src.opt_generator import opt_generator
    from time import time
    import pandas as pd
    import multiprocessing as mp
    import os
    import shelve
    import copy
    import argparse
    import shutil
    import sys
    import logging
    import PyQt5

    # Patch library errors, currently for ete3
    patch.patch()

    # For starting time stamp
    start_time = time()

    # Recursion limit exceeds for big trees
    sys.setrecursionlimit(1000000)

    # Argument parsing in commandline mode
    args = CommandParser().get_args()

    ############################################################################
    #     Start of initializing blocks should not be moved for function!!!     #
    ############################################################################

    # initialize working directories, and opts
    opt, path, list_info, list_warning, list_error = initialize.initialize(
        path_run=os.getcwd(), parser=args
    )
    tool.initialize_path(path)

    # setup logging
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        encoding="utf-8",
        level=tool.get_level(opt.verbose),
        handlers=[logging.FileHandler(path.log), logging.StreamHandler()],
    )

    # Delayed logging for option parsing
    for info in list_info:
        logging.info(info)

    for warning in list_warning:
        logging.warning(warning)

    for error in list_error:
        logging.error(error)

    # V contains all intermediate variables for FunID Run
    V = dataset.FunID_var()
    # R includes all reporting objects such as warnings, errors, statistics
    R = reporter.Report()

    ##########################################################################
    #     End of initializing blocks should not be moved for function!!!     #
    ##########################################################################

    ##########################################################################
    #                        Start of main Pipeline                          #
    ##########################################################################

    # IO stage
    if opt.continue_from_previous is False or index_step(opt.step) <= 0:
        step = "setup"
        logging.info("SETUP")

        # save metadata input to report
        # get input data
        V, R, opt = manage_input.data_input(V, R, opt, path)

        # get possible genus list for recognizing genus name
        V = initialize.get_genus_list(V, opt, path)

        R.update_report(V=V, path=path, opt=opt, step=step)
        save.save_session(opt=opt, path=path, global_var=globals(), var=vars())

    # Searching (BLAST or mmseqss)
    if opt.continue_from_previous is False or index_step(opt.step) <= 1:
        step = "search"
        logging.info("SEARCHING")

        # Generate BLAST or mmseqs matrices for further analysis
        V = search.search_df(V, path, opt)

        # Concatenate search matrix among genes
        V = multigene.concatenate_df(V, path, opt)

        R.update_report(V=V, path=path, opt=opt, step=step)
        save.save_session(opt=opt, path=path, global_var=globals(), var=vars())

    # Clustering
    if opt.continue_from_previous is False or index_step(opt.step) <= 2:
        step = "cluster"
        logging.info("CLUSTERING")

        # Assign group
        logging.info("ASSIGNING GROUPS")

        V, opt, path = cluster.pipe_cluster(V, opt, path)

        # Append query column to each of the df, and concatenated_df in df_dict as query assigned
        V = cluster.append_query_group(V)

        # Append query column to concatenated_df as query assigned
        # V = cluster.append_concatenated_query_group(V)

        # Both gene and group was assigned, dataset was confirmed
        V.generate_dataset(opt)

        # Appending outgroup
        logging.info("Appending outgroup")
        # For non-concatenated outgroup
        # ready for multiprocessing run
        # Pushing all v to multiprocessing requires too much memory
        V, path, opt = cluster.pipe_append_outgroup(V, path, opt)

        # remove invalid dataset from downstream analysis
        V.remove_invalid_dataset()
        # Save dataset to outgroups
        V.save_dataset(path, opt)

        R.update_report(V=V, path=path, opt=opt, step=step)
        save.save_session(opt=opt, path=path, global_var=globals(), var=vars())

    # Alignment
    if opt.continue_from_previous is False or index_step(opt.step) <= 3:
        step = "align"
        logging.info("ALIGNMENT")
        V, path, opt = align.pipe_alignment(V, path, opt)

        R.update_report(V=V, path=path, opt=opt, step=step)
        save.save_session(opt=opt, path=path, global_var=globals(), var=vars())

    # Trimming
    if opt.continue_from_previous is False or index_step(opt.step) <= 4:
        step = "trim"
        logging.info("TRIMMING")
        V, path, opt = trim.pipe_trimming(V, path, opt)

        R.update_report(V=V, path=path, opt=opt, step=step)
        save.save_session(opt=opt, path=path, global_var=globals(), var=vars())

    # Concatenation (multigene)
    if opt.continue_from_previous is False or index_step(opt.step) <= 5:
        step = "concatenate"
        logging.info("CONCATENATING MULTIGENE ALIGNMENTS")

        # Multigene
        V = multigene.combine_alignment(V, opt, path)

        R.update_report(V=V, path=path, opt=opt, step=step)
        save.save_session(opt=opt, path=path, global_var=globals(), var=vars())

    # Alignment validations - whether some of the sequences does not have overlapping regions
    V.validate_alignments(path, opt)

    # Modeltest
    if opt.continue_from_previous is False or index_step(opt.step) <= 6:
        step = "modeltest"
        logging.info("MODELTEST")

        # most of the modeltest algorithms are well working with multiprocessing
        model_dict = modeltest.modeltest(V, path, opt)

        R.update_report(V=V, path=path, opt=opt, step=step)
        save.save_session(opt=opt, path=path, global_var=globals(), var=vars())

    # Tree construction
    if opt.continue_from_previous is False or index_step(opt.step) <= 7:
        step = "tree"
        logging.info("TREE CONSTRUCTION")

        # Build phylogenetic trees
        V, path, opt = tree.pipe_tree(V, path, opt, model_dict)

        R.update_report(V=V, path=path, opt=opt, step=step)
        save.save_session(opt=opt, path=path, global_var=globals(), var=vars())

        # move tree files
        tool.cleanup_tree(path)

        R.update_report(V=V, path=path, opt=opt, step=step)
        save.save_session(opt=opt, path=path, global_var=globals(), var=vars())

    # Visualize
    if opt.continue_from_previous is False or index_step(opt.step) <= 8:
        step = "visualize"
        logging.info("VISUALIZE")

        # Running visualization
        V, path, opt = tree_interpretation_pipe.pipe_tree_interpretation(V, path, opt)

        # move tree image files
        tool.cleanup_tree_image(path)

        R.update_report(V=V, path=path, opt=opt, step=step)
        save.save_session(opt=opt, path=path, global_var=globals(), var=dir())

    # Report
    if opt.continue_from_previous is False or index_step(opt.step) <= 9:
        step = "report"
        logging.info("REPORT")

        # Reporting
        V.homogenize_dataset()

        R.update_report(V, path, opt, step=step)
        save.save_session(opt=opt, path=path, global_var=globals(), var=vars())

        end_time = time()
        logging.info(f"FunID ended in {end_time - start_time}")
        print(f"FunID ended in {end_time - start_time}")

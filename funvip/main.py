def main():
    from funvip.src import align
    from funvip.src import tree_interpretation_pipe
    from funvip.src import cluster
    from funvip.src import dataset
    from funvip.src import ext
    from funvip.src import hasher
    from funvip.src import initialize
    from funvip.src import logger
    from funvip.src import validate_input
    from funvip.src import modeltest
    from funvip.src import concatenate
    from funvip.src import patch
    from funvip.src import reporter
    from funvip.src import tool
    from funvip.src import save
    from funvip.src import search
    from funvip.src import tree
    from funvip.src import trim
    from funvip.src.command import CommandParser
    from funvip.src.tool import index_step
    from funvip.src.opt_generator import opt_generator
    from funvip.src.version import Version
    from time import time
    from time import sleep
    import pandas as pd
    import subprocess
    import multiprocessing as mp
    import os
    import shelve
    import copy
    import argparse
    import shutil
    import sys
    import logging
    import PyQt5

    # To initialize log print
    _ = subprocess.call("", shell=True)

    if "FunID" in sys.argv[0]:
        print(
            "\x1b[92m The name of the software changed to 'FunVIP'. Please use 'FunVIP' when citation. \x1b[0m"
        )
        sleep(3)

    # For starting time stamp
    time_start = time()

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

    # Initialize version management
    version = Version(opt, path)

    logger.setup_logging(list_info, list_warning, list_error, path, opt, tool)

    # V contains all intermediate variables for FunVIP Run
    V = dataset.FunVIP_var()
    # R includes all reporting objects such as warnings, errors, statistics
    R = reporter.Report()

    # Reload previous session from shelve if --continue selected
    if opt.continue_from_previous is True:
        var = save.load_session(opt, savefile=path.save)
        if "V" in var:
            V = var["V"]
        if "R" in var:
            R = var["R"]
        if "path" in var:
            path = var["path"]
        if "model_dict" in var:
            model_dict = var["model_dict"]

    # To deal with location changes when rerun
    path_root = f"{os.getcwd()}/{opt.runname}"
    path.init_workspace(path_root, opt)

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
        # R.initialize_metadata(opt)
        # get input data
        V, R, opt, GenMine_flag = validate_input.data_input(V, R, opt, path)

        # get possible genus list for recognizing genus name
        V = initialize.get_genus_list(V, opt, path)

        R.update_report(
            V=V,
            path=path,
            opt=opt,
            step=step,
            version=version,
            GenMine_flag=GenMine_flag,
        )

        # logging works well from here

        save.save_session(opt=opt, path=path, global_var=locals(), var=vars())

        time_setup = time()

        if not ("HS49HE" in V.dict_hash_FI):
            raise Exception

    # Searching (BLAST or mmseqs)
    if opt.continue_from_previous is False or index_step(opt.step) <= 1:
        step = "search"
        logging.info("SEARCHING")

        ## Generate BLAST or mmseqs matrices for further analysis
        # Also, gene informations are updated in this step
        V = search.search_df(V, path, opt)

        # Concatenate search matrix among genes
        V = concatenate.concatenate_df(V, path, opt)

        # Update gene assignment of query sequence inputs
        validate_input.update_queryfile(V, path, opt)

        R.update_report(
            V=V,
            path=path,
            opt=opt,
            step=step,
            version=version,
            GenMine_flag=GenMine_flag,
        )
        save.save_session(opt=opt, path=path, global_var=locals(), var=vars())

        time_search = time()

        if not ("HS49HE" in V.dict_hash_FI):
            raise Exception

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

        # Both gene and group was assigned, dataset are confirmed
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

        R.update_report(
            V=V,
            path=path,
            opt=opt,
            step=step,
            version=version,
            GenMine_flag=GenMine_flag,
        )
        save.save_session(opt=opt, path=path, global_var=locals(), var=vars())

        time_cluster = time()

        if not ("HS49HE" in V.dict_hash_FI):
            raise Exception

    # Alignment
    if opt.continue_from_previous is False or index_step(opt.step) <= 3:
        step = "align"
        logging.info("ALIGNMENT")
        V, path, opt = align.pipe_alignment(V, path, opt)

        R.update_report(
            V=V,
            path=path,
            opt=opt,
            step=step,
            version=version,
            GenMine_flag=GenMine_flag,
        )
        save.save_session(opt=opt, path=path, global_var=locals(), var=vars())

        time_align = time()

        if not ("HS49HE" in V.dict_hash_FI):
            raise Exception

    # Trimming
    if opt.continue_from_previous is False or index_step(opt.step) <= 4:
        step = "trim"
        logging.info("TRIMMING")
        V, path, opt = trim.pipe_trimming(V, path, opt)
        # Alignment validations - whether some of the sequences does not have overlapping regions
        V.validate_alignments(path=path, opt=opt)
        R.update_report(
            V=V,
            path=path,
            opt=opt,
            step=step,
            version=version,
            GenMine_flag=GenMine_flag,
        )
        save.save_session(opt=opt, path=path, global_var=locals(), var=vars())

        time_trim = time()

        if not ("HS49HE" in V.dict_hash_FI):
            raise Exception

    # Concatenation (multigene)
    if opt.continue_from_previous is False or index_step(opt.step) <= 5:
        step = "concatenate"
        logging.info("CONCATENATING MULTILOCI ALIGNMENTS")

        # Multigene
        V = concatenate.combine_alignment(V, opt, path)

        R.update_report(
            V=V,
            path=path,
            opt=opt,
            step=step,
            version=version,
            GenMine_flag=GenMine_flag,
        )
        save.save_session(opt=opt, path=path, global_var=locals(), var=vars())

        time_concatenate = time()

        if not ("HS49HE" in V.dict_hash_FI):
            raise Exception

    # Modeltest
    if opt.continue_from_previous is False or index_step(opt.step) <= 6:
        step = "modeltest"
        logging.info("MODELTEST")

        # most of the modeltest algorithms are well working with multiprocessing
        model_dict = modeltest.modeltest(V, path, opt)

        R.update_report(
            V=V,
            path=path,
            opt=opt,
            step=step,
            version=version,
            GenMine_flag=GenMine_flag,
        )
        save.save_session(opt=opt, path=path, global_var=locals(), var=vars())

        time_modeltest = time()

        if not ("HS49HE" in V.dict_hash_FI):
            raise Exception

    # Tree construction
    if opt.continue_from_previous is False or index_step(opt.step) <= 7:
        step = "tree"
        logging.info("TREE CONSTRUCTION")

        # Build phylogenetic trees
        V, path, opt = tree.pipe_tree(V, path, opt, model_dict)

        R.update_report(
            V=V,
            path=path,
            opt=opt,
            step=step,
            version=version,
            GenMine_flag=GenMine_flag,
        )
        save.save_session(opt=opt, path=path, global_var=locals(), var=vars())

        # move tree and image files
        tool.cleanup_tree(path)
        tool.cleanup_tree_image(path)

        R.update_report(
            V=V,
            path=path,
            opt=opt,
            step=step,
            version=version,
            GenMine_flag=GenMine_flag,
        )
        save.save_session(opt=opt, path=path, global_var=locals(), var=vars())

        time_tree = time()

        if not ("HS49HE" in V.dict_hash_FI):
            raise Exception

    # Visualize
    if opt.continue_from_previous is False or index_step(opt.step) <= 8:
        step = "visualize"
        logging.info("VISUALIZE")

        # The genus duplication occurs here
        # Running visualization
        V, path, opt = tree_interpretation_pipe.pipe_tree_interpretation(V, path, opt)

        # move tree image files
        tool.cleanup_tree(path)
        tool.cleanup_tree_image(path)

        R.update_report(
            V=V,
            path=path,
            opt=opt,
            step=step,
            version=version,
            GenMine_flag=GenMine_flag,
        )
        save.save_session(opt=opt, path=path, global_var=locals(), var=dir())

        time_visualize = time()

        if not ("HS49HE" in V.dict_hash_FI):
            raise Exception

    # Report
    if opt.continue_from_previous is False or index_step(opt.step) <= 9:
        step = "report"
        logging.info("REPORT")

        # Reporting
        V.homogenize_dataset()

        # Check inconsistent
        V.check_inconsistent()

        R.update_report(
            V, path, opt, step=step, version=version, GenMine_flag=GenMine_flag
        )
        save.save_session(opt=opt, path=path, global_var=locals(), var=vars())

        time_end = time()
        logging.info(f"FunVIP ended in {time_end - time_start}")

        if not ("HS49HE" in V.dict_hash_FI):
            raise Exception

        try:
            logging.info(f"Time setup : {time_setup-time_start}s")
        except:
            logging.warning(f"Failed logging initialization time")

        try:
            logging.info(f"Time search : {time_search-time_setup}s")
        except:
            logging.warning(f"Failed logging search time")

        try:
            logging.info(f"Time cluster : {time_cluster-time_search}s")
        except:
            logging.warning(f"Failed logging cluster time")

        try:
            logging.info(f"Time alignment : {time_align-time_cluster}s")
        except:
            logging.warning(f"Failed logging alignment time")

        try:
            logging.info(f"Time trimming : {time_trim-time_align}s")
        except:
            logging.warning(f"Failed logging trimming time")

        try:
            logging.info(f"Time concatenatation : {time_concatenate-time_trim}s")
        except:
            logging.warning(f"Failed logging concatenation time")

        try:
            logging.info(f"Time tree construction : {time_tree-time_concatenate}s")
        except:
            logging.warning(f"Failed logging tree construction time")

        try:
            logging.info(f"Time tree interpretation : {time_visualize-time_tree}s")
        except:
            logging.warning(f"Failed logging tree interpretation time")

        try:
            logging.info(f"Time report generation: {time_end-time_visualize}s")
        except:
            logging.warning(f"Failed logging reoprt generation time")

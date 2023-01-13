def main():

    from funid.src import (
        align,
        tree_interpretation_pipe,
        cluster,
        dataset,
        ext,
        hasher,
        initialize,
        logger,
        manage_input,
        modeltest,
        multigene,
        patch,
        reporter,
        tool,
        save,
        search,
        tree,
        trim,
    )

    from funid.src.command import CommandParser
    from funid.src.tool import initialize_path
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

    # Declaring available steps
    step = [
        "setup",
        "search",
        "cluster",
        "align",
        "trim",
        "concatenate",
        "modeltest",
        "tree",
        "visualize",
        "report",
    ]

    ############################################################################
    #     Start of initializing blocks should not be moved for function!!!     #
    ############################################################################

    # initialize working directories, and opts
    opt, path, list_info, list_warning, list_error = initialize.initialize(
        path_run=os.getcwd(), parser=args
    )
    initialize_path(path)

    # setup logging
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        encoding="utf-8",
        level=logging.DEBUG,
        handlers=[logging.FileHandler(path.log), logging.StreamHandler()],
    )

    # Delayed option parsing logging
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
    # Initialize logging process
    R.log.initialize_log(step)

    ##########################################################################
    #     End of initializing blocks should not be moved for function!!!     #
    ##########################################################################

    ##########################################################################
    #                        Start of main Pipeline                          #
    ##########################################################################

    # IO stage
    if opt.continue_from_previous is False or step.index(opt.step) <= 0:
        logging.info("Set up")

        # save metadata input to report
        R.initialize_metadata(opt)
        # input data
        V, R, opt = manage_input.data_input(V, R, opt, path)
        # get possible genus list for recognizing genus name
        V = initialize.get_genus_list(V, opt, path)

        save.save_session(opt=opt, path=path, global_var=globals(), var=vars())

    # Searching (BLAST or mmseqss)
    if opt.continue_from_previous is False or step.index(opt.step) <= 1:

        logging.info("SEARCHING")

        # Generate BLAST or mmseqs matrices for further analysis
        V = search.search_df(V, path, opt)
        # Concatenate search matrix among genes
        V = multigene.concatenate_df(V, path, opt)
        save.save_session(opt=opt, path=path, global_var=globals(), var=vars())

    # Clustering
    if opt.continue_from_previous is False or step.index(opt.step) <= 2:
        logging.info("CLUSTERING")

        # Assign section
        logging.info("ASSIGNING SECTION")

        V, opt, path = cluster.pipe_cluster(V, opt, path)

        # Save assigned section info in table format
        save.save(V.list_FI, path, opt)

        # Append query column to each of the df in df_dict as query assigned
        V = cluster.append_query_section(V)

        # Append query column to concatenated_df as query assigned
        if opt.concatenate is True:
            V = cluster.append_concatenated_query_section(V)

        # Both gene and section was assigned, dataset was confirmed
        V.generate_dataset(opt)

        # Appending outgroup
        logging.info("Appending outgroup")

        # For non-concatenated outgroup
        # ready for multiprocessing run
        # Pushing all v to multiprocessing requires too much memory
        V, path, opt = cluster.pipe_append_outgroup(V, path, opt)

        # Wrap up outgroup appending results
        # outgroup result : (gene, section, outgroup_list, section_list)
        V = cluster.outgroup_result_collector(V)
        # remove invalid dataset to be analyzed ()
        V.remove_invalid_dataset()
        # Save dataset to outgroups
        V.save_dataset(path, opt)

        save.save_session(opt=opt, path=path, global_var=globals(), var=vars())

    # Alignment
    if opt.continue_from_previous is False or step.index(opt.step) <= 3:

        logging.info("ALIGNMENT")
        V, path, opt = align.pipe_alignment(V, path, opt)
        save.save_session(opt=opt, path=path, global_var=globals(), var=vars())

    # Trimming
    if opt.continue_from_previous is False or step.index(opt.step) <= 4:
        logging.info("TRIMMING")
        V, path, opt = trim.pipe_trimming(V, path, opt)
        save.save_session(opt=opt, path=path, global_var=globals(), var=vars())

    # Concatenation (multigene)
    if opt.continue_from_previous is False or step.index(opt.step) <= 5:
        logging.info("CONCATENATING MULTIGENE ALIGNMENTS")
        if opt.concatenate is True:
            # Multigene
            V = multigene.combine_alignment(V, opt, path)
            save.save_session(opt=opt, path=path, global_var=globals(), var=vars())
        else:
            logging.info("Passing concatenated tree")

    # Alignment validations - whether some of the sequences does not have overlapping regions
    V.validate_alignments(path, opt)

    # Modeltest
    if opt.continue_from_previous is False or step.index(opt.step) <= 6:
        logging.info("MODELTEST")

        # most of the modeltest algorithms are well working with multiprocessing
        model_dict = modeltest.modeltest(V, path, opt)
        save.save_session(opt=opt, path=path, global_var=globals(), var=vars())

    # Tree construction
    if opt.continue_from_previous is False or step.index(opt.step) <= 7:
        logging.info("TREE CONSTRUCTION")

        # Build phylogenetic trees
        V, path, opt = tree.pipe_tree(V, path, opt, model_dict)
        save.save_session(opt=opt, path=path, global_var=globals(), var=vars())

        # move tree image files
        tool.cleanup_tree(path)
        save.save_session(opt=opt, path=path, global_var=globals(), var=vars())

    # Visualize
    if opt.continue_from_previous is False or step.index(opt.step) <= 8:
        logging.info("VISUALIZE")

        # Running visualization
        V, path, opt = tree_interpretation_pipe.pipe_tree_interpretation(V, path, opt)

        # move tree image files
        tool.cleanup_tree_image(path)

        # raise Exception
        save.save_session(opt=opt, path=path, global_var=globals(), var=dir())

    # Report
    if opt.continue_from_previous is False or step.index(opt.step) <= 9:
        logging.info("REPORT")

        writer = pd.ExcelWriter(
            f"{path.root}/{opt.runname}_report.xlsx", engine="xlsxwriter"
        )

        # raw_result: (out, df, report_list)
        # out: f"{opt.runname}_{sect}_{gene}"
        # df: db, query, others, total dataframe
        # report_list: report: SingleReport - accession, hash, section, gene, inputtaxon, identifiedtaxon

        # Reporting
        V.homogenize_dataset()
        R.statistics.update_statistics(V, opt)
        R.arrange(V, opt, f"{path.root}/{opt.runname}_Identification_result.xlsx")

        # Testing
        R.initialize_section_report(V, opt)
        R.render(opt, out=f"{path.root}/{opt.runname}_report.html")

        """
        try:
            R.initialize_section_report(V, opt)
            R.render(opt, out=f"{path.root}/{opt.runname}_report.html")
        except:
            logging.error(
                f"Data report is currently experimental feature. Failed in this analysis"
            )
        """

        save.save_session(opt=opt, path=path, global_var=globals(), var=vars())

        end_time = time()
        logging.info(f"FunID ended in {end_time - start_time}")

# Option generator for multiprocessing starmap
import logging
import os


def opt_generator(V, opt, path, step, thread=None):
    list_opt = []

    if thread is None:
        thread = opt.thread

    ## alignment options generator
    if step == "alignment":
        for group in V.dict_dataset:
            for gene in V.dict_dataset[group]:
                # For concatenated gene matrix, alignment should be done by each gene matrix and then concatenated
                if gene != "concatenated":
                    # double checking path
                    if os.path.isfile(
                        f"{path.out_adjusted}/{opt.runname}_Adjusted_{group}_{gene}.fasta"
                    ):
                        list_opt.append(
                            (
                                f"{path.out_adjusted}/{opt.runname}_Adjusted_{group}_{gene}.fasta",
                                f"{path.out_alignment}/{opt.runname}_MAFFT_{group}_{gene}.fasta",
                                path,
                                thread,
                                opt.mafft.algorithm,
                                "adjustdirection",
                                opt.mafft.op,
                                opt.mafft.ep,
                            )
                        )

    elif step == "trimming":
        # Generate trimming opts for multiprocessing
        for group in V.dict_dataset:
            for gene in V.dict_dataset[group]:
                if not (gene == "concatenated"):
                    if opt.method.trim.lower() == "gblocks":
                        list_opt.append(
                            (
                                f"{path.out_alignment}/{opt.runname}_MAFFT_{group}_{gene}.fasta",
                                f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta",
                                path,
                                opt,
                            )
                        )
                    elif opt.method.trim.lower() == "trimal":
                        list_opt.append(
                            (
                                f"{path.out_alignment}/{opt.runname}_MAFFT_{group}_{gene}.fasta",
                                f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta",
                                path,
                                opt,
                            )
                        )
                    else:  # for just copy
                        list_opt.append(
                            (
                                f"{path.out_alignment}/{opt.runname}_MAFFT_{group}_{gene}.fasta",
                                f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta",
                                path,
                                opt,
                            )
                        )

    # For tcs alignment validation
    elif step == "tcs":
        for group in V.dict_dataset:
            for gene in V.dict_dataset[group]:
                if gene != "concatenated":
                    tcs_out = (
                        f"{path.out_alignment}/tcs/{opt.runname}_{group}_{gene}.tcs"
                    )
                    list_opt.append(
                        (
                            f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta",
                            thread,
                            tcs_out,
                        )
                    )

    else:
        logging.error(f"[Error] Unexpected step {step} given for opt_generator")
        raise Exception

    return list_opt

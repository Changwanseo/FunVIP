# Option generator for multiprocessing starmap

# from .logger import Mes
import logging
import os


def opt_generator(V, opt, path, step):

    list_opt = []

    # alignment options generator
    if step == "alignment":
        for group in V.dict_dataset:
            for gene in V.dict_dataset[group]:
                # thread assignment
                if opt.verbose < 3:
                    if opt.thread in (1, 3):
                        thread = 1
                    else:
                        thread = 2
                else:
                    thread = opt.thread

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
                            )
                        )
                    elif opt.method.trim.lower() == "trimal":
                        list_opt.append(
                            (
                                f"{path.out_alignment}/{opt.runname}_MAFFT_{group}_{gene}.fasta",
                                f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta",
                                path,
                                opt.trimal.algorithm,
                                opt.trimal.gt,
                            )
                        )
                    else:  # for just copy
                        list_opt.append(
                            (
                                f"{path.out_alignment}/{opt.runname}_MAFFT_{group}_{gene}.fasta",
                                f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta",
                            )
                        )

    else:
        logging.error(f"[Error] Unexpected step {step} given for opt_generator")
        raise Exception

    return list_opt

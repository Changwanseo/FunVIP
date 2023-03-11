from funid.src import ext
from funid.src.opt_generator import opt_generator
from Bio import SeqIO
import multiprocessing as mp
import os
import logging
import shutil


def pipe_trimming(V, path, opt):
    trimming_opt = opt_generator(V, opt, path, step="trimming")

    # run multiprocessing start
    if opt.verbose < 3:
        p = mp.Pool(opt.thread)
        if opt.method.trim.lower() == "gblocks":
            trimming_result = p.starmap(ext.Gblocks, trimming_opt)
        elif opt.method.trim.lower() == "trimal":
            trimming_result = p.starmap(ext.Trimal, trimming_opt)
        else:
            trimming_result = p.starmap(shutil.copy, trimming_opt)
        p.close()
        p.join()

    else:
        # non-multithreading mode for debugging
        trimming_result = []
        for option in trimming_opt:
            if opt.method.trim.lower() == "gblocks":
                trimming_result.append(ext.Gblocks(*option))
            elif opt.method.trim.lower() == "trimal":
                trimming_result.append(ext.Trimal(*option))
            else:
                trimming_result.append(shutil.copy(*option))

    # Remove datasets if trimming results nothing
    trim_fail = []
    for group in V.dict_dataset:
        for gene in V.dict_dataset[group]:
            if not gene == "concatenated":
                if os.path.isfile(
                    f"{path.out_alignment}/{opt.runname}_MAFFT_{group}_{gene}.fasta"
                ):
                    if os.path.isfile(
                        f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta"
                    ):
                        seq_list = list(
                            SeqIO.parse(
                                f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta",
                                "fasta",
                            )
                        )
                        if len(seq_list[0].seq) == 0:
                            trim_fail.append((group, gene))
                    else:
                        pass
                else:
                    pass

    for fail in trim_fail:
        group = fail[0]
        gene = fail[1]
        V.dict_dataset[group].pop(gene)
        logging.warning(
            f"Dataset {group} {gene} has removed during trimming because no sequence left. Please check alignment if database contains bad sequences."
        )

        # If non of the genes left for group except for concatenate, remove group
        if len(V.dict_dataset[group]) == 1 and "concatenated" in V.dict_dataset[group]:
            logging.warning(
                f"Dataset {group} has removed because non of the genes are available"
            )
            V.dict_dataset.pop(group)

    return V, path, opt

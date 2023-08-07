from funid.src import ext
from funid.src.opt_generator import opt_generator
from Bio import AlignIO
import multiprocessing as mp
import os
import logging
import shutil


## Trimming function of FunID
# As most of the trimming functions affects inside of alignment
def trimming(alignment, out, path, opt):
    # Save alignment before running trimming
    original_msa = AlignIO.read(alignment, "fasta")

    # Running trimming
    if opt.method.trim.lower() == "gblocks":
        trimming_result = ext.Gblocks(fasta=alignment, out=out, path=path)
    elif opt.method.trim.lower() == "trimal":
        trimming_result = ext.Trimal(
            fasta=alignment,
            out=out,
            path=path,
            algorithm=opt.trimal.algorithm,
            threshold=opt.trimal.gt,
        )
    else:
        trimming_result = shutil.copy(fasta, out)

    # Repair mid alignment
    if not (opt.allow_innertrimming):
        # Repaired trimmend alignment
        trimmed_msa = AlignIO.read(out, "fasta")
        removed_columns = []
        original_length = len(original_msa[0])
        trimmed_length = len(trimmed_msa[0])

        if original_length > trimmed_length:
            for col_index in range(original_length):
                original_column = original_msa[:, col_index]
                trimmed_column = trimmed_msa[:, col_index]
                if original_column != trimmed_column:
                    removed_columns.append(col_index)

        survived = set(range(original_length)) - set(col_index)

        revived_msa = original_msa[:, min(survived) : max(survived) + 1]
        print(f"Conserved MSA: {min(survived)}:{max(survived) + 1}")

        AlignIO.write(revived_msa, out, "fasta")

    return trimming_result


# Trimming pipeline of FunID
def pipe_trimming(V, path, opt):
    trimming_opt = opt_generator(V, opt, path, step="trimming")
    # run multiprocessing start
    if opt.verbose < 3:
        p = mp.Pool(opt.thread)
        trimming_result = p.starmap(trimming, trimming_opt)
        p.close()
        p.join()

    else:
        # non-multithreading mode for debugging
        trimming_result = []
        for option in trimming_opt:
            trimming_result.append(trimming(*option))

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

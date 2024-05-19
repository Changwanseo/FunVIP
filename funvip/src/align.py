from funvip.src import ext
from funvip.src.opt_generator import opt_generator
import logging
from Bio import SeqIO
import multiprocessing as mp
import math


# Each module to be run in alignment multiprocessing
def module_alignment(
    in_fasta,
    out_fasta,
    path,
    thread,
    mafft_algorithm,
    adjustdirection,
    op,
    ep,
):
    # Running MAFFT
    ext.MAFFT(
        fasta=in_fasta,
        out=out_fasta,
        path=path,
        thread=thread,
        algorithm=mafft_algorithm,
        adjust=adjustdirection,
        maxiterate=1000,
        op=op,
        ep=ep,
    )

    # to prevent reversed sequence making error
    fasta_list = list(SeqIO.parse(out_fasta, "fasta"))
    for seq in fasta_list:
        if seq.description.startswith("_R_"):
            seq.id = ""
            seq.description = seq.description[3:]
    SeqIO.write(fasta_list, out_fasta, "fasta")

    # Fix unexpected spaces on mafft
    with open(out_fasta, "r") as fr:
        alignment = fr.read()

    with open(out_fasta, "w") as fw:
        fw.write(alignment.replace(" HS", "HS"))


# Alignment pipeline
def pipe_alignment(V, path, opt):
    alignment_opt = opt_generator(V, opt, path, step="alignment")

    # Multiprocessing start
    # Thread optimizations
    # Using O(N^2L^2) of MAFFT-G-ins-i
    # Use multithreading optimizations in non-verbose mode
    if opt.verbose < 3:
        time_consumptions = []
        ## Simulate time comsumption
        # Parsing probable time_consumption
        for aln_opt in alignment_opt:
            before_aln_file = aln_opt[0]
            before_seq_list = list(SeqIO.parse(before_aln_file, "fasta"))
            N = len(before_seq_list)
            L = max([len(seq.seq) for seq in before_seq_list])
            time_consumptions.append(N * N * L * L)

        time_consumptions = sorted(time_consumptions, reverse=True)
        best_distribute = 1
        best_thread = opt.thread
        best_time = 99999999999999999999
        # Optimize thread numbers
        for distribute_num in range(opt.thread):
            time_buffers = {}

            # Threads to be distributed
            each_thread_num = int(opt.thread / (distribute_num + 1))
            # Efficiency multipler estimated from
            # Rubio-Largo, A., Castelli, M., Vanneschi, L., & Vega-Rodríguez, M. A. (2018). A parallel multiobjective metaheuristic for multiple sequence alignment. Journal of Computational Biology, 25(9), 1009-1022.

            efficiency_multiplier = 1.0258 * math.log(each_thread_num) + 1

            for i in range(distribute_num + 1):
                time_buffers[i + 1] = 0

            for t in time_consumptions:
                time_buffers[min(time_buffers, key=time_buffers.get)] += int(
                    t / efficiency_multiplier
                )

            total_time_consumption = max(time_buffers.values())
            logging.debug(
                f"Time buffers for {distribute_num+1} workers: {time_buffers}"
            )
            logging.debug(total_time_consumption)

            if total_time_consumption < best_time:
                best_thread = each_thread_num
                best_distribute = distribute_num + 1
                best_time = total_time_consumption

        # logging.info(f"time_consumption_simulation: {time_buffers}")

        logging.info(
            f"Using {best_thread} threads for {best_distribute} workers in alignment by optimization"
        )

        # Update thread numbers by simulation
        alignment_opt = opt_generator(
            V, opt, path, step="alignment", thread=best_thread
        )
        p = mp.Pool(best_distribute)

        alignment_result = p.starmap(module_alignment, alignment_opt)
        p.close()
        p.join()

    else:
        # non-multithreading mode for debugging
        alignment_result = []
        for option in alignment_opt:
            alignment_result.append(module_alignment(*option))

    return V, path, opt

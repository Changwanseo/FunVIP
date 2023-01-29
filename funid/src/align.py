from funid.src import ext
from funid.src.opt_generator import opt_generator
from Bio import SeqIO
import multiprocessing as mp

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

    # run multiprocessing start
    # Thread optimizations
    if opt.verbose < 3:
        if opt.thread in (1, 3):
            p = mp.Pool(int(opt.thread))
        else:
            p = mp.Pool(int(opt.thread / 2))

        alignment_result = p.starmap(module_alignment, alignment_opt)
        p.close()
        p.join()

    else:
        # non-multithreading mode for debugging
        alignment_result = []
        for option in alignment_opt:
            alignment_result.append(module_alignment(*option))

    return V, path, opt

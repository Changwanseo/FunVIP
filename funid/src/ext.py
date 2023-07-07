# for running external programs
from sys import platform
from Bio import SeqIO
import logging
import os, subprocess
import shutil
import psutil
from copy import deepcopy
from pathlib import Path
from funid.src.save import save_tree
from funid.src.tool import mkdir


# Search methods
# BLAST
def blast(query, db, out, path, opt):
    path_blast = Path(f"{path.sys_path}/external/BLAST_Windows/bin/blastn.exe")

    # quotations make errors on windows platform
    if platform == "win32":
        CMD = f"{path_blast} -out {out} -query {query} -outfmt 6 -db {db} -word_size {opt.cluster.wordsize} -evalue {opt.cluster.evalue} -num_threads {opt.thread}"
    else:
        CMD = f"blastn -out '{out}' -query '{query}' -outfmt 6 -db '{db}' -word_size {opt.cluster.wordsize} -evalue {opt.cluster.evalue} -num_threads {opt.thread}"

    logging.info(CMD)
    Run = subprocess.call(CMD, shell=True)


# mmseqs
def mmseqs(query, db, out, tmp, path, opt):
    path_mmseqs = f"{path.sys_path}/external/mmseqs_Windows/mmseqs.bat"

    if platform == "win32":
        CMD = f"{path_mmseqs} easy-search {query} {db} {out} {tmp} --threads {opt.thread} -k {opt.cluster.wordsize} --search-type 3 -e {opt.cluster.evalue} --dbtype 2"
    else:
        CMD = f"mmseqs easy-search '{query}' '{db}' '{out}' '{tmp}' --threads {opt.thread} -k {opt.cluster.wordsize} --search-type 3 -e {opt.cluster.evalue} --dbtype 2"

    logging.info(CMD)
    Run = subprocess.call(CMD, shell=True)


# DB building methods
def makeblastdb(fasta, db, path):
    path_makeblastdb = f"{path.sys_path}/external/BLAST_Windows/bin/makeblastdb.exe"

    # To prevent makeblastdb error in windows, run it on temporary directory and move it
    if platform == "win32":
        # Save original path
        ori_path = deepcopy(os.getcwd())
        makeblastdb_path = f"{path.tmp}\\makeblastdb\\"
        # remove temporary path if exists
        if os.path.exists(makeblastdb_path):
            shutil.rmtree(makeblastdb_path)
        # make new temp makeblastdb directory
        mkdir(makeblastdb_path)
        os.chdir(makeblastdb_path)

        # Move makeblastdb.exe and destination file to temporate directory
        shutil.copy(fasta, makeblastdb_path)

        # Remove disk seperator to prevent error
        fasta_tmp = fasta.replace("\\", "/").split("/")[-1]
        db_tmp = db.replace("\\", "/").split("/")[-1]

        # run make blast db
        CMD = f"{path_makeblastdb} -in {fasta_tmp} -blastdb_version 4 -title {db_tmp} -dbtype nucl"
        logging.info(CMD)
        Run = subprocess.call(CMD, shell=True)
        # Change db names
        shutil.move(fasta_tmp + ".nsq", db + ".nsq")
        shutil.move(fasta_tmp + ".nin", db + ".nin")
        shutil.move(fasta_tmp + ".nhr", db + ".nhr")
        # return to original path
        os.chdir(ori_path)
        # remove temporary path
        if os.path.exists(makeblastdb_path):
            shutil.rmtree(makeblastdb_path)
    else:
        CMD = f"makeblastdb -in '{fasta}' -blastdb_version 4 -title '{db}' -dbtype nucl"
        logging.info(CMD)
        return_code = subprocess.call(CMD, shell=True)

        if return_code != 0:
            logging.error(f"[ERROR] Make blast_db failed!!")
            install_flag = 1

        # Change db names
        shutil.move(fasta + ".nsq", db + ".nsq")
        shutil.move(fasta + ".nin", db + ".nin")
        shutil.move(fasta + ".nhr", db + ".nhr")


def makemmseqsdb(fasta, db, path):
    path_makemmseqsdb = f"{path.sys_path}/external/mmseqs_Windows/mmseqs.bat"

    if platform == "win32":
        CMD = f"{path_makemmseqsdb} createdb {fasta} {db} --createdb-mode 0 --dbtype 2"
    else:
        CMD = f"mmseqs createdb '{fasta}' '{db}' --createdb-mode 0 --dbtype 2"
    logging.info(CMD)
    Run = subprocess.call(CMD, shell=True)


# Alignments
def MAFFT(
    fasta,
    out,
    path,
    thread=1,
    algorithm="localpair",
    adjust="adjustdirection",
    maxiterate=1000,
    op=1.3,
    ep=0.1,
):
    # validate if there are only 1 sequence
    seqlist = list(SeqIO.parse(fasta, "fasta"))
    if len(seqlist) == 1:
        logging.warning(
            f"{fasta} has only one sequence. Using original sequence as alignment"
        )
        shutil.copy(fasta, out)
    else:
        if platform == "win32":
            CMD = f"{path.sys_path}/external/MAFFT_Windows/mafft-win/mafft.bat --thread {thread} --{algorithm} --maxiterate {maxiterate} --{adjust} --op {op} --ep {ep} --quiet {fasta} > {out}"
        else:
            CMD = f"mafft --thread {thread} --{algorithm} --maxiterate {maxiterate} --{adjust} --op {op} --ep {ep} --quiet '{fasta}' > '{out}'"

        logging.info(CMD)
        try:
            Run = subprocess.call(CMD, shell=True)
        except:
            logging.error(f"Failed on {CMD}")
            raise Exception


# Trimming
def Gblocks(fasta, out, path):
    if platform == "win32":
        CMD = f"{path.sys_path}/external/Gblocks_Windows_0.91b/Gblocks_0.91b/Gblocks.exe {fasta} -t=d -b4=2 -b5=a -e=.gb"
    else:
        CMD = f"Gblocks '{fasta}' -t=d -b4=2 -b5=a -e=.gb"

    logging.info(CMD)
    Run = subprocess.call(CMD, shell=True)

    try:
        shutil.move(f"{fasta}.gb", out)
        shutil.move(f"{fasta}.gb.htm", f"{out}.html")
    except:  # when only one sequence and Gblocks failed
        shutil.move(fasta, out)


def Trimal(fasta, out, path, algorithm="gt", threshold=0.2):
    if algorithm == "gt":
        algorithm = f"{algorithm} {threshold}"

    if platform == "win32":
        CMD = f"{path.sys_path}/external/trimal.v1.2rev59/trimAl/bin/trimal.exe -in {fasta} -out {out} -{algorithm}"

    else:
        CMD = f"trimal -in {fasta} -out {out} -{algorithm}"

    logging.info(CMD)
    Run = subprocess.call(CMD, shell=True)

    # to remove unexpected hash
    fasta_list = list(SeqIO.parse(out, "fasta"))

    for seq in fasta_list:
        # if " " in seq.description:
        seq.id = seq.description.split(" ")[0]
        seq.description = ""

    SeqIO.write(fasta_list, out, "fasta")


# Modeltest
def Modeltest_ng(fasta, out, models, thread):
    if platform == "win32":
        logging.error("Modeltest-NG is not available in windows. Try IQTREE modeltest")
        raise Exception
    else:
        CMD = f"modeltest-ng -i '{fasta}' -o '{out}' -t ml -p {thread} --disable-checkpoint {models}"

    logging.info(CMD)
    Run = subprocess.call(CMD, shell=True)


# IQTREE ModelFinder
def ModelFinder(fasta, opt, path, thread):
    if opt.method.tree == "iqtree":
        model_term = "-m MFP"
    elif opt.method.tree == "raxml":
        model_term = "-m MF --mset raxml"
    elif opt.method.tree == "fasttree":
        model_term = "-m MF --mset JC,JC+G4,GTR,GTR+G4"
    else:
        logging.error(
            f"Modelterm cannot be selected to tree method {opt.method.tree} while running modelfinder"
        )
        raise Exception

    if platform == "win32":
        CMD = f"{path.sys_path}/external/iqtree-2.1.3-Windows/bin/iqtree2.exe --seqtype DNA -s {fasta} {model_term} -merit {opt.criterion} -T {thread} -mem {opt.memory}"
    else:
        # not final
        CMD = f"iqtree --seqtype DNA -s '{fasta}' {model_term} -merit {opt.criterion} -T {thread} -mem {opt.memory}"
    logging.info(CMD)
    Run = subprocess.call(CMD, shell=True)


# Tree building
def RAxML(
    fasta,
    out,
    hash_dict,
    path,
    thread=1,
    bootstrap=100,
    partition=None,
    model="-m GTRGAMMA",
):
    if model == "skip":
        model = ""

    # Because RAxML does not allows out location, change directory for running
    path_ori = os.getcwd()
    os.chdir(path.tmp)

    if platform == "win32":
        CMD = f"{path.sys_path}/external/RAxML_Windows/raxmlHPC-PTHREADS-AVX2.exe -s {fasta} -n {out} -p 1 -T {thread} -f a -# {bootstrap} -x 1 {model}"
    elif platform == "darwin":
        CMD = f"raxmlHPC-PTHREADS -s '{fasta}' -n '{out}' -p 1 -T {thread} -f a -# {bootstrap} -x 1 {model}"
    else:
        CMD = f"raxmlHPC-PTHREADS-AVX -s '{fasta}' -n '{out}' -p 1 -T {thread} -f a -# {bootstrap} -x 1 {model}"

    if not (partition is None):
        CMD += f" -q {partition}"

    logging.info(CMD)
    Run = subprocess.call(CMD, shell=True)

    # Return result to original directory
    os.chdir(path_ori)
    file = out.split("/")[-1]
    out = f"RAxML_bipartitions.{out}"
    save_tree(
        out=f"{path.tmp}/{out}",
        hash_dict=hash_dict,
        hash_file_path=f"{path.out_tree}/hash_{file}",
        decoded_file_path=f"{path.out_tree}/{file}",
    )


def FastTree(fasta, out, hash_dict, path, model=""):
    if model == "skip":
        model = ""
    if platform == "win32":
        CMD = f"{path.sys_path}/external/FastTree_Windows/FastTree.exe -quiet -nt {model} -log {path.tmp}/fasttreelog -seed 1 {fasta} > {path.tmp}/{out}"
    else:
        CMD = f"FastTree -quiet -nt {model} -log {path.tmp}/fasttreelog -seed 1 '{fasta}' > {path.tmp}/{out}"

    logging.info(CMD)
    Run = subprocess.call(CMD, shell=True)
    file = out.split("/")[-1]
    save_tree(
        out=f"{path.tmp}/{out}",
        hash_dict=hash_dict,
        hash_file_path=f"{path.out_tree}/hash_{file}",
        decoded_file_path=f"{path.out_tree}/{file}",
        fix=True,
    )


def IQTREE(
    fasta,
    out,
    hash_dict,
    path,
    memory=f"{max(2, int(psutil.virtual_memory().total / (1024**3)))}G",
    thread=1,
    bootstrap=1000,
    partition=None,
    model="",
):
    if model == "skip":
        model = ""

    if bootstrap < 1000:
        logging.warning("IQTREE requires at least 1000 bootstrap, setting to 1000")
        bootstrap = 1000

    if platform == "win32":
        CMD = f"{path.sys_path}/external/iqtree-2.1.3-Windows/bin/iqtree2.exe -s {fasta} -B {bootstrap} -T {thread} {model}"
    else:
        CMD = f"iqtree -s {fasta} -B {bootstrap} -T {thread} {model}"

    # Partitioned analysis cannot be used with memory
    if not (partition is None):
        CMD += f" -q {partition}"
    else:
        CMD += f" -mem {memory}"

    logging.info(CMD)
    Run = subprocess.call(CMD, shell=True)
    try:
        if partition is None:
            shutil.move(f"{fasta}.treefile", f"{path.tmp}/{out}")
        else:
            shutil.move(f"{partition}.contree", f"{path.tmp}/{out}")
    except:
        logging.error(
            "IQTREE FAILED. Maybe due to memory problem if partitioned analysis included."
        )

    file = out.split("/")[-1]
    save_tree(
        out=f"{path.tmp}/{out}",
        hash_dict=hash_dict,
        hash_file_path=f"{path.out_tree}/hash_{file}",
        decoded_file_path=f"{path.out_tree}/{file}",
    )

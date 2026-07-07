# for running external programs
from sys import platform
from Bio import SeqIO
import logging
import os, subprocess
import shutil
import psutil
from copy import deepcopy
from pathlib import Path
from funvip.src.save import save_tree
from funvip.src.tool import mkdir
from funvip.src.exceptions import ExternalToolError, SearchError


# Search methods
# BLAST
def blast(query, db, out, path, opt):
    path_blast = Path(f"{path.sys_path}/external/BLAST_Windows/bin/blastn.exe")

    # quotations make errors on windows platform when space does not exists
    if platform == "win32":
        if " " in out:
            out = f'"{out}"'
        if " " in query:
            query = f'"{query}"'
        if " " in db:
            db = f'"{db}"'

        CMD = f"{path_blast} -out {out} -query {query} -outfmt 6 -db {db} -word_size {opt.cluster.wordsize} -evalue {opt.cluster.evalue} -num_threads {opt.thread} -max_target_seqs {opt.cluster.max_target_seqs}"
    else:
        CMD = f"blastn -out '{out}' -query '{query}' -outfmt 6 -db '{db}' -word_size {opt.cluster.wordsize} -evalue {opt.cluster.evalue} -num_threads {opt.thread} -max_target_seqs {opt.cluster.max_target_seqs}"

    logging.info(CMD)
    Run = subprocess.call(CMD, shell=True)
    if Run != 0:
        raise SearchError(f"blastn failed (exit code {Run}) for query {query} vs db {db}")


# mmseqs
def mmseqs(query, db, out, tmp, path, opt):
    path_mmseqs = f"{path.sys_path}/external/mmseqs_Windows/mmseqs.bat"

    if platform == "win32":
        if " " in out:
            out = f'"{out}"'
        if " " in query:
            query = f'"{query}"'
        if " " in db:
            db = f'"{db}"'
        if " " in tmp:
            tmp = f'"{tmp}"'
        CMD = f"{path_mmseqs} easy-search {query} {db} {out} {tmp} --threads {opt.thread} -k {opt.cluster.wordsize} --search-type 3 -e {opt.cluster.evalue} --dbtype 2 --max-seqs {opt.cluster.max_target_seqs}"
    else:
        CMD = f"mmseqs easy-search '{query}' '{db}' '{out}' '{tmp}' --threads {opt.thread} -k {opt.cluster.wordsize} --search-type 3 -e {opt.cluster.evalue} --dbtype 2  --max-seqs {opt.cluster.max_target_seqs}"

    logging.info(CMD)
    Run = subprocess.call(CMD, shell=True)
    if Run != 0:
        raise SearchError(f"mmseqs easy-search failed (exit code {Run}) for query {query} vs db {db}")


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
        # I cannot find any "quiet" options for makeblastdb
        Run = subprocess.call(CMD, stdout=subprocess.DEVNULL, shell=True)
        if Run != 0:
            raise ExternalToolError(f"makeblastdb failed (exit code {Run}) for {fasta}")
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
        # I cannot find any "quiet" options for makeblastdb
        return_code = subprocess.call(CMD, stdout=subprocess.DEVNULL, shell=True)

        if return_code != 0:
            raise ExternalToolError(f"makeblastdb failed (exit code {return_code}) for {fasta}")

        # Change db names
        shutil.move(fasta + ".nsq", db + ".nsq")
        shutil.move(fasta + ".nin", db + ".nin")
        shutil.move(fasta + ".nhr", db + ".nhr")


def makemmseqsdb(fasta, db, path):
    path_makemmseqsdb = f"{path.sys_path}/external/mmseqs_Windows/mmseqs.bat"

    if " " in fasta:
        fasta = f'"{fasta}"'

    if " " in db:
        db = f'"{db}"'

    if platform == "win32":
        CMD = f"{path_makemmseqsdb} createdb {fasta} {db} --createdb-mode 0 --dbtype 2"
    else:
        CMD = f"mmseqs createdb '{fasta}' '{db}' --createdb-mode 0 --dbtype 2"
    logging.info(CMD)
    Run = subprocess.call(CMD, shell=True)
    if Run != 0:
        raise ExternalToolError(f"mmseqs createdb failed (exit code {Run}) for {fasta}")


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
            if " " in out:
                out = f'"{out}"'
            if " " in fasta:
                fasta = f'"{fasta}"'

            CMD = f"{path.sys_path}/external/MAFFT_Windows/mafft-win/mafft.bat --thread {thread} --{algorithm} --maxiterate {maxiterate} --{adjust} --op {op} --ep {ep} --quiet {fasta} > {out}"
        else:
            CMD = f"mafft --thread {thread} --{algorithm} --maxiterate {maxiterate} --{adjust} --op {op} --ep {ep} --quiet '{fasta}' > '{out}'"

        logging.info(CMD)
        Run = subprocess.call(CMD, shell=True)
        if Run != 0:
            raise ExternalToolError(f"MAFFT failed (exit code {Run}) for {fasta}")


# Trimming
def Gblocks(fasta, out, path):
    if platform == "win32":
        if " " in fasta:
            fasta = f'"{fasta}"'
        CMD = f"{path.sys_path}/external/Gblocks_Windows_0.91b/Gblocks_0.91b/Gblocks.exe {fasta} -t=d -b4=2 -b5=a -e=.gb -p=t"
    else:
        CMD = f"Gblocks '{fasta}' -t=d -b4=2 -b5=a -e=.gb -p=t"

    logging.info(CMD)
    Run = subprocess.call(CMD, shell=True)

    try:
        shutil.move(f"{fasta}.gb", out)
    except FileNotFoundError:  # when only one sequence and Gblocks failed
        shutil.move(fasta, out)

    # Parse and return column statistics. (Gblocks returns a nonzero exit code even
    # on success, so its exit code is intentionally not checked.) Default to the
    # (-1, -1) failure sentinel so a missing/format-drifted .gb.txt cannot NameError.
    start_pos = -1
    end_pos = -1
    if os.path.exists(f"{fasta}.gb.txt"):
        with open(f"{fasta}.gb.txt", "r") as f:
            for line in f:
                if line.startswith("Flanks:"):
                    flank_log = (
                        line.replace("Flanks:", "")
                        .replace("  ", " ")
                        .replace("[", "")
                        .replace("]", "")
                        .strip()
                        .split(" ")
                    )
                    logging.debug(f"Gblocks flanks for {fasta}: {flank_log}")
                    try:
                        flank_log = [int(x) for x in flank_log]
                        start_pos = flank_log[0]
                        end_pos = flank_log[-1] - 1
                    except (ValueError, IndexError):
                        start_pos = -1
                        end_pos = -1

    logging.debug(f"Gblocks start_pos: {start_pos}, end_pos: {end_pos}")

    try:
        shutil.move(f"{fasta}.gb.txt", path.extlog)
    except FileNotFoundError:
        pass

    return (start_pos, end_pos)


def Trimal(fasta, out, path, algorithm="gt", threshold=0.2):
    if algorithm == "gt":
        algorithm = f"{algorithm} {threshold}"

    if platform == "win32":
        if " " in fasta:
            fasta = f'"{fasta}"'
        if " " in out:
            out_dir = f'"{out}"'
            out_colnumbering = f'"{out}.colnumbering"'
        else:
            out_dir = out
            out_colnumbering = f"{out}.colnumbering"

        CMD = f"{path.sys_path}/external/trimal.v1.4/trimAl/bin/trimal.exe -in {fasta} -out {out_dir} -{algorithm} -terminalonly -colnumbering > {out_colnumbering}"

    else:
        CMD = f"trimal -in '{fasta}' -out '{out}' -{algorithm} -terminalonly -colnumbering > '{out}.colnumbering'"

    logging.info(CMD)
    Run = subprocess.call(CMD, shell=True)
    if Run != 0 or not os.path.exists(out):
        raise ExternalToolError(f"trimal failed (exit code {Run}) for {fasta}")

    # to remove unexpected hash included - maybe not needed after stabilization
    fasta_list = list(SeqIO.parse(out, "fasta"))

    for seq in fasta_list:
        # if " " in seq.description:
        seq.id = seq.description.split(" ")[0]
        seq.description = ""

    SeqIO.write(fasta_list, out, "fasta")

    # Parse and return column statistics
    with open(f"{out}.colnumbering", "r") as f:
        line = f.read()
        cols = line.replace("#ColumnsMap", "").strip().split(", ")
        try:
            cols = [int(x) for x in cols]
            start_pos = cols[0]
            end_pos = cols[-1]
        except (ValueError, IndexError):
            start_pos = -2
            end_pos = -2

    try:
        shutil.move(f"{out}.colnumbering", path.extlog)
    except FileNotFoundError:
        pass

    # Trimal uses 0 based position, return with +1
    return (start_pos + 1, end_pos + 1)


# Modeltest
def Modeltest_ng(fasta, out, path, models, thread):
    path_modeltestng = f"{path.sys_path}/external/modeltest-ng_Windows/modeltest-ng.exe"
    if platform == "win32":
        CMD = f"{path_modeltestng} -i '{fasta}' -o '{out}' -t ml -p {thread} --disable-checkpoint {models}"
        """
        logging.error("Modeltest-NG is not available in windows. Try IQTREE modeltest")
        raise Exception
        """

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
        raise ExternalToolError(
            f"cannot select a model term for tree method {opt.method.tree!r} in ModelFinder"
        )

    if platform == "win32":
        if " " in fasta:
            fasta = f'"{fasta}"'
        CMD = f"{path.sys_path}/external/iqtree/bin/iqtree2.exe --seqtype DNA -s {fasta} {model_term} -merit {opt.criterion} -nt AUTO -ntmax {thread} -mem {opt.memory} --quiet"
    else:
        # not final
        CMD = f"iqtree --seqtype DNA -s '{fasta}' {model_term} -merit {opt.criterion} -nt AUTO -ntmax {thread} -mem {opt.memory} --quiet"
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
    version="old",
):
    if model == "skip":
        model = ""

    # Because RAxML does not allow an out location, change directory for running.
    # Wrap in try/finally so the original cwd is always restored, even if RAxML
    # (or CMD selection) fails -- otherwise later relative-path operations break.
    path_ori = os.getcwd()
    os.chdir(path.tmp)
    try:
        if platform == "win32":
            if " " in fasta:
                fasta = f'"{fasta}"'
            if " " in out:
                out = f'"{out}"'

            CMD = f"{path.sys_path}/external/RAxML_Windows/raxmlHPC-PTHREADS-AVX2.exe -s {fasta} -n {out} -p 1 -T {thread} -f a -# {bootstrap} -x 1 {model} --silent"
        elif platform == "darwin":
            # For Rosetta
            if (
                subprocess.run(
                    "raxmlHPC-PTHREADS -v",
                    shell=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.STDOUT,
                ).returncode
                == 0
            ):
                CMD = f"raxmlHPC-PTHREADS -s '{fasta}' -n '{out}' -p 1 -T {thread} -f a -# {bootstrap} -x 1 {model}"
            # For arm native
            elif (
                subprocess.run(
                    "raxmlHPC -v",
                    shell=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.STDOUT,
                ).returncode
                == 0
            ):
                CMD = f"raxmlHPC -s '{fasta}' -n '{out}' -p 1 -T {thread} -f a -# {bootstrap} -x 1 {model}"
            else:
                raise ExternalToolError(
                    "cannot find a working RAxML for this Apple Silicon system"
                )

        else:
            if version == "old":
                CMD = f"raxmlHPC-PTHREADS-AVX -s '{fasta}' -n '{out}' -p 1 -T {thread} -f a -# {bootstrap} -x 1 {model} --silent"
            elif version == "new":
                CMD = f"raxmlHPC-PTHREADS-AVX2 -s '{fasta}' -n '{out}' -p 1 -T {thread} -f a -# {bootstrap} -x 1 {model} --silent"
            else:
                raise ExternalToolError(f"unexpected RAxML version {version!r}")

        if not (partition is None):
            CMD += f" -q {partition}"

        logging.info(CMD)
        Run = subprocess.call(CMD, shell=True)
        if Run != 0:
            raise ExternalToolError(f"RAxML failed (exit code {Run}) for {fasta}")
    finally:
        # Always restore the original working directory
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
        """
        if " " in model:
            model = f"{model}"
        """
        if " " in fasta:
            fasta = f'"{fasta}"'
        if " " in path.tmp:
            path_tmp = f'"{path.tmp}/fasttreelog"'
        else:
            path_tmp = f"{path.tmp}/fasttreelog"
        if " " in path.tmp or " " in out:
            path_out = f'"{path.tmp}/{out}"'
        else:
            path_out = f"{path.tmp}/{out}"
        CMD = f"{path.sys_path}/external/FastTree_Windows/FastTree.exe -quiet -nt {model} -log {path_tmp} -seed 1 {fasta} > {path_out}"
    else:
        CMD = f"FastTree -quiet -nt {model} -log {path.tmp}/fasttreelog -seed 1 '{fasta}' > {path.tmp}/{out}"

    logging.info(CMD)
    Run = subprocess.call(CMD, shell=True)
    tree_file = f"{path.tmp}/{out}"
    if Run != 0 or not os.path.exists(tree_file) or os.path.getsize(tree_file) == 0:
        raise ExternalToolError(f"FastTree failed (exit code {Run}) for {fasta}")
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
        # For working with space
        if " " in fasta:
            tmp_fasta = f'"{fasta}"'
        else:
            tmp_fasta = fasta

        CMD = f"{path.sys_path}/external/iqtree/bin/iqtree2.exe -s {tmp_fasta} -B {bootstrap} -nt AUTO -ntmax {thread} {model}"
    else:
        CMD = f"iqtree -s {fasta} -B {bootstrap} -nt AUTO -ntmax {thread} {model}"

    logging.info(f"partition: {partition}")
    # Partitioned analysis cannot be used with memory option
    if not (partition is None):
        if " " in partition:
            tmp_partition = f'"{partition}"'
        else:
            tmp_partition = partition
        CMD += f" -q {tmp_partition}"
    else:
        CMD += f" -mem {memory}"

    logging.info(CMD)
    Run = subprocess.call(CMD, shell=True)
    try:
        if partition is None:
            shutil.move(f"{fasta}.contree", f"{path.tmp}/{out}")
        else:
            shutil.move(f"{partition}.contree", f"{path.tmp}/{out}")
    except FileNotFoundError as e:
        raise ExternalToolError(
            f"IQTREE failed (exit code {Run}, no .contree produced) for {fasta}; "
            "possibly a memory problem for partitioned analysis"
        ) from e

    file = out.split("/")[-1]
    save_tree(
        out=f"{path.tmp}/{out}",
        hash_dict=hash_dict,
        hash_file_path=f"{path.out_tree}/hash_{file}",
        decoded_file_path=f"{path.out_tree}/{file}",
    )


# TCS calculation from T-COFFEE
def TCS(fasta, thread, out):
    if platform == "win32":
        raise ExternalToolError("TCS (from T-COFFEE) is only available on Linux")
    else:
        # T_COFFEE env variable MAX_N_PID_4_TCOFFEE should be changed for 64bit machine
        # should be already done in installation check process
        # subprocess.call("export MAX_N_PID_4_TCOFFEE=4194304", shell=True)
        CMD = f"t_coffee -infile {fasta} -cpu {thread} -method fast_pair -type DNA -evaluate -output score_ascii -outfile {out} -quiet"

    logging.info(CMD)
    # Even though "quiet" option exists, TCS show some blank lines
    # Run = subprocess.call(CMD, stdout=open(os.devnull, "wb"), shell=True)
    Run = subprocess.run(
        CMD, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, shell=True
    ).returncode

    if Run != 0:
        raise ExternalToolError(f"TCS (T-COFFEE) failed (exit code {Run})")

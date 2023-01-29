from funid.src import cluster, tool, hasher, manage_input
from funid.src.ext import blast, makeblastdb, mmseqs, makemmseqsdb
from funid.src.manage_input import save_df
import copy
import pandas as pd
import numpy as np
import os
import logging
import subprocess
import shutil
import hashlib


def cleanblastdb(db):
    os.remove(f"{db}.nsq")
    os.remove(f"{db}.nin")
    os.remove(f"{db}.nhr")


def cleanmmseqsdb(db):
    os.remove(f"{db}.dbtype")
    os.remove(f"{db}_h.dbtype")
    os.remove(f"{db}.index")
    os.remove(f"{db}_h.index")
    os.remove(f"{db}.lookup")
    os.remove(f"{db}.source")
    os.remove(f"{db}")
    os.remove(f"{db}_h")


# Make blast db file or mmseqs db file
def create_search_db(opt, db_fasta, db, path) -> None:
    # make blastdb
    if opt.method.search.lower() in ("blast", "blastn"):
        makeblastdb(fasta=db_fasta, db=db, path=path)
    # make mmseqssdb
    elif opt.method.search.lower() in (
        "mmseqss",
        "mmseqs",
        "mmseqs2",
        "mmseqss2",
    ):
        makemmseqsdb(fasta=db_fasta, db=db, path=path)
    else:
        logging.error("DEVELOPMENTAL ERROR on building search DB!")
        raise Exception


# Merge fragmented search matches from given blast or mmseqss results
def merge_fragments(df) -> pd.DataFrame():

    # Check if list (or pandas series or one-column DataFrame) has only one value
    def get_unique(series) -> str:
        if len(set(series)) == 1:
            return list(series)[0]
        elif len(set(series)) == 0:
            logging.error(f"Found 0 values in {series} while merging search results")
            raise Exception
        else:
            logging.error(
                f"Found {len(set(series))} values in {series} while merging search results"
            )
            raise Exception

    # Calculate overall percent identity of fragments
    def calculate_pident(df):
        pident = df["pident"]
        length = df["length"]
        if len(pident) != len(length):
            logging.error(
                f"During merge blast fragments, found pident {pident} and len {length} are different"
            )
            raise Exception
        elif len(pident) == 0:
            logging.error(f"No percent identitiy {pident} found")
            raise Exception
        else:
            # Calculate overall percent identity
            return np.sum(np.array(pident) * np.array(length)) / np.sum(length)

    # Return empty df because if causes error
    if len(df.index) == 0:  # faster way for df.empty
        return df

    # pident needs access to other columns, will be calculated in next line
    # qseqid and sseqid will return itself, after checking if they are unique
    # mismatch, gaps, and bitscore were calculated as sum of fragments
    # evalues were calculated by multiplying them, because they are probability

    df = df.groupby(["qseqid", "sseqid"], dropna=True, as_index=False).aggregate(
        {
            "qseqid": lambda x: set(x),
            "sseqid": lambda x: set(x),
            "pident": lambda x: tuple(x),
            "length": lambda x: tuple(x),
            "mismatch": np.sum,
            "gaps": np.sum,
            "qstart": lambda x: tuple(x),
            "qend": lambda x: tuple(x),
            "sstart": lambda x: tuple(x),
            "send": lambda x: tuple(x),
            "evalue": np.prod,
            "bitscore": np.sum,
        }
    )

    # Calculate pident
    df["pident"] = df.apply(calculate_pident, axis=1)

    # Merge length
    df["length"] = df["length"].apply(sum).astype(int)

    # Update qseqid and sseqid
    df["qseqid"] = df["qseqid"].apply(lambda x: tuple(x)[0])
    df["sseqid"] = df["sseqid"].apply(lambda x: tuple(x)[0])

    return df


# Core search : Running blast or mmseqss
def search(query_fasta, db_fasta, path, opt) -> pd.DataFrame():

    # Working on saved database
    # find key for which db to use
    if opt.cachedb is True or opt.usecache is True:

        # get hash number of the given db to compare
        hash_ = hashlib.md5(open(db_fasta, "rb").read()).hexdigest()

        if hash_ is None:
            logging.error(f"Database file {db_fasta} missing")
            raise Exception

        # Try to parse DB
        if opt.usecache is True:

            # When succesfully parsed DB
            if (
                os.path.isdir(f"{path.in_db}/{opt.method.search.lower()}/{hash_}")
                is True
            ):
                logging.info("[INFO] Found existing database! Skipping database build")
                db = f"{path.in_db}/{opt.method.search.lower()}/{hash_}/{hash_}"

            # When parsing existing DB failed
            else:

                logging.info("No existing database found")

                # Try save DB
                if opt.cachedb is True:
                    logging.info(
                        f"--cachedb selected, {opt.method.search.lower()} database will be saved"
                    )

                    # Create DB saving directory
                    os.mkdir(f"{path.in_db}/{opt.method.search.lower()}/{hash_}")
                    db = f"{path.in_db}/{opt.method.search.lower()}/{hash_}/{hash_}"

                    # DB saving starts
                    logging.info("The database is in first run, caching database")

                    # Create search database
                    create_search_db(opt, db_fasta, db, path)

                # Passing Save DB
                else:
                    logging.info(
                        "--cachedb not selected, saving database will be passed"
                    )
                    db = f"{path.tmp}/{opt.runname}"

                    # Create search database
                    create_search_db(opt, db_fasta, db, path)

    else:
        db = f"{path.tmp}/{opt.runname}"
        # Create search database
        create_search_db(opt, db_fasta, db, path)

    # run searching
    # blast
    if opt.method.search.lower() in ("blast", "blastn"):
        blast_result = blast(
            query=query_fasta,
            db=db,
            out=f"{path.tmp}/{opt.runname}.m8",
            path=path,
            opt=opt,
        )
        # remove db when saving not enabled
        if opt.cachedb is False:
            cleanblastdb(db)
    # mmseqs
    elif opt.method.search.lower() in ("mmseqs", "mmseq", "mmseq2", "mmseqs2"):
        mmseqs(
            query=query_fasta,
            db=db,
            out=f"{path.tmp}/{opt.runname}.m8",
            tmp=path.tmp,
            path=path,
            opt=opt,
        )
        # remove db when saving not enabled
        if opt.cachedb is False:
            cleanmmseqsdb(db)
        # remove temporary file
        # shutil.rmtree(f"{path.tmp}/{opt.runname}")
    else:
        logging.error("DEVELOPMENTAL ERROR on searching!")
        raise Exception

    # Parse out
    df = pd.read_csv(
        f"{path.tmp}/{opt.runname}.m8",
        sep="\t",
        header=None,
        names=[
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "mismatch",
            "gaps",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
        ],
    )

    print(df)

    # Remove temporary files
    os.remove(f"{path.tmp}/{opt.runname}.m8")

    # merge fragmented matches
    df = merge_fragments(df)

    return df


# Main search: Get blast or mmseqss dataframe results from given dataset
def search_df(V, path, opt):

    # Initialize variable for search result
    dict_search = {}

    # Ready for by sseqid hash, which group to append
    # generating group dict to assign group column next to the output dataframe
    group_dict = {}
    for FI in V.list_FI:
        group_dict[FI.hash] = FI.group

    # genes available for db and query
    V.list_db_gene = copy.copy(opt.gene)
    V.list_qr_gene = copy.copy(opt.gene)

    # make blastdb for this run
    list_db_FI = tool.select(V.list_FI, datatype="db")

    # Make database by genes
    for gene in opt.gene:

        db_state = manage_input.save_fasta(
            list_db_FI, gene, f"{path.tmp}/{opt.runname}_DB_{gene}.fasta", by="hash"
        )

        if db_state == 0:
            logging.warning(
                f"No data for {gene} found in database. Excluding from analysis"
            )
            V.list_db_gene.remove(gene)

    # get query fasta from funinfo_list
    list_qr_FI = tool.select(V.list_FI, datatype="query")

    # hash_dict format - hash : id
    V.dict_id_hash = hasher.encode(V.list_FI)

    # if no query exists and only database sequence exists
    if len(list_qr_FI) == 0:
        logging.warning("No query selected. Changing to database validation mode")
        # db by db search analysis
        for gene in V.list_db_gene:
            df_search = search(
                query_fasta=f"{path.tmp}/{opt.runname}_DB_{gene}.fasta",
                db_fasta=f"{path.tmp}/{opt.runname}_DB_{gene}.fasta",
                path=path,
                opt=opt,
            )

            # append group column
            df_search["subject_group"] = df_search["sseqid"].apply(
                lambda x: group_dict.get(x)
            )
            # add to search dict
            # dict_search[gene] = df_search
            V.dict_gene_SR[gene] = df_search

            # Save dataframe
            if opt.savesearchresult is True:
                save_df(
                    hasher.decode_df(V.dict_id_hash, df_search),
                    f"{path.out_matrix}/{opt.runname}_BLAST_result_{gene}.{opt.matrixformat}",
                    fmt=opt.matrixformat,
                )

    # if query sequence exists
    else:
        query_state = manage_input.save_fasta(
            list_qr_FI,
            "unclassified",
            f"{path.tmp}/{opt.runname}_Query_unclassified.fasta",
            by="hash",
        )

        # first, run for unclassified query for assign gene
        dict_unclassified = {}

        if query_state == 1:  # if unclassified sequence exists
            for gene in V.list_db_gene:
                df_search = search(
                    query_fasta=f"{path.tmp}/{opt.runname}_Query_unclassified.fasta",
                    db_fasta=f"{path.tmp}/{opt.runname}_DB_{gene}.fasta",
                    path=path,
                    opt=opt,
                )

                # if dataframe exists
                if isinstance(df_search, pd.DataFrame):
                    if not df_search.empty:
                        dict_unclassified[gene] = df_search

            # assign gene by search result to unassigned sequences
            V = cluster.assign_gene(dict_unclassified, V)

        # then, gene by gene BLAST for outgroup
        for gene in opt.gene:
            # if database is very confident, so no more validation is required
            if opt.confident is True:
                query_state = manage_input.save_fasta(
                    [FI for FI in V.list_FI if FI.datatype == "query"],
                    gene,
                    f"{path.tmp}/{opt.runname}_Query_{gene}.fasta",
                    by="hash",
                )

            # for most of the cases
            else:
                query_state = manage_input.save_fasta(
                    V.list_FI,
                    gene,
                    f"{path.tmp}/{opt.runname}_Query_{gene}.fasta",
                    by="hash",
                )

            # if no sequence exists for gene, remove it
            if query_state == 0:
                V.list_qr_gene.remove(gene)
            else:
                # If db column and query column does not matches -> should be moved to manage_input
                if not (gene in V.list_db_gene):
                    logging.error(
                        f"Gene {gene} found in query, but not found in database. Please add {gene} column to database"
                    )
                    raise Exception

        # BLAST or mmseqs search
        # This part should be changed by using former search result for faster performance

        for gene in V.list_qr_gene:

            df_search = search(
                query_fasta=f"{path.tmp}/{opt.runname}_Query_{gene}.fasta",
                db_fasta=f"{path.tmp}/{opt.runname}_DB_{gene}.fasta",
                path=path,
                opt=opt,
            )

            # append group column
            df_search["subject_group"] = df_search["sseqid"].apply(
                lambda x: group_dict.get(x)
            )
            V.dict_gene_SR[gene] = df_search

            # Save dataframe
            if opt.savesearchresult is True:
                save_df(
                    hasher.decode_df(V.dict_id_hash, df_search),
                    f"{path.out_matrix}/{opt.runname}_BLAST_result_{gene}.{opt.matrixformat}",
                    opt.matrixformat,
                )

        # remove tmp file after search
        os.remove(f"{path.tmp}/{opt.runname}_Query_unclassified.fasta")

    # remove temporary files
    for gene in opt.gene:
        try:
            os.remove(f"{path.tmp}/{opt.runname}_Query_{gene}.fasta")
        except:
            pass
        try:
            os.remove(f"{path.tmp}/{opt.runname}_DB_{gene}.fasta")
        except:
            pass

    return V

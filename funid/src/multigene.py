from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import logging
from functools import reduce
import pandas as pd
from funid.src import search, hasher


def combine_alignment(V, opt, path):

    multigene_list = []
    for group in V.dict_dataset:
        if "concatenated" in V.dict_dataset[group]:
            # get alignment length
            len_dict = {}  # length of each of the alignments
            seq_dict = {}
            hash_set = set()
            gene_list = []

            for gene in V.dict_dataset[group]:
                if not (gene == "concatenated"):

                    fasta_list = list(
                        SeqIO.parse(
                            f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta",
                            "fasta",
                        )
                    )

                    len_dict[gene] = len(fasta_list[0].seq)
                    seq_dict[gene] = {}

                    total_dataset = (
                        [
                            FI.hash
                            for FI in V.dict_dataset[group]["concatenated"].list_qr_FI
                        ]
                        + [
                            FI.hash
                            for FI in V.dict_dataset[group]["concatenated"].list_db_FI
                        ]
                        + [
                            FI.hash
                            for FI in V.dict_dataset[group]["concatenated"].list_og_FI
                        ]
                    )

                    for seq in fasta_list:
                        if seq.description in total_dataset:  # if available hash
                            seq_dict[gene][seq.description] = seq
                            hash_set.add(seq.description)

            # Generate partition file
            with open(
                f"{path.out_alignment}/{opt.runname}_{group}.partition", "w"
            ) as fw:
                tot_len = 0
                for gene in sorted(len_dict.keys()):

                    gene_list.append(gene)
                    fw.write(f"DNA, {gene}= {tot_len+1}-{tot_len+len_dict[gene]}\n")
                    tot_len += len_dict[gene]

            # Generate concatenated alignment file
            concatenate_list = []
            for hash_id in hash_set:
                tmp_seq = ""
                for gene in gene_list:
                    if hash_id in seq_dict[gene]:
                        tmp_seq += str(seq_dict[gene][hash_id].seq)
                    else:
                        tmp_seq += "-" * len_dict[gene]  # add gaps for ids without gene

                concatenate_list.append(
                    SeqRecord(id=hash_id, description="", seq=Seq(tmp_seq))
                )

            SeqIO.write(
                concatenate_list,
                f"{path.out_alignment}/{opt.runname}_trimmed_{group}_concatenated.fasta",
                "fasta",
            )

            SeqIO.write(
                concatenate_list,
                f"{path.out_alignment}/{opt.runname}_MAFFT_{group}_concatenated.fasta",
                "fasta",
            )

        else:
            logging.warning(
                f"For {group}, no more than genes were detected. Passing combining alignment"
            )

    V.multigene_list = multigene_list

    return V


# for concatenating blast result to get single bitscore among results
def concatenate_df(V, path, opt):

    if opt.concatenate is True:

        df_list = [V.dict_gene_SR[gene] for gene in V.dict_gene_SR]

        if len(df_list) <= 0:
            logging.warning(
                f"Stop concatenating because same or less than 0 gene exists"
            )
            return None

        else:
            # drop unused columns
            cnt = 0

            # initializing because search result for specific gene may not exists
            while 1:
                if isinstance(df_list[cnt], pd.DataFrame):
                    if not (df_list[cnt].empty):
                        V.cSR = df_list[cnt]
                        V.cSR.drop(
                            columns=[
                                "pident",
                                "length",
                                "mismatch",
                                "gaps",
                                "qstart",
                                "qend",
                                "sstart",
                                "send",
                                "evalue",
                            ],
                            inplace=True,
                        )
                        cnt += 1
                        break

                cnt += 1

            # generate concatenated search result
            for df in df_list[cnt:]:
                if isinstance(df_list[cnt], pd.DataFrame):
                    if not (df_list[cnt].empty):
                        df.drop(
                            columns=[
                                "pident",
                                "length",
                                "mismatch",
                                "gaps",
                                "qstart",
                                "qend",
                                "sstart",
                                "send",
                                "evalue",
                            ],
                            inplace=True,
                        )
                        V.cSR = pd.merge(
                            V.cSR,
                            df,
                            how="outer",
                            on=["sseqid", "qseqid"],
                            suffixes=("1", "2"),
                        )

                        drop_list = []
                        rename_dict = {}

                        # concatenate bitscore columns
                        if "bitscore1" in V.cSR.columns:
                            V.cSR["bitscore1"].fillna(value=0, inplace=True)
                            V.cSR["bitscore2"].fillna(value=0, inplace=True)
                            V.cSR["bitscore"] = V.cSR["bitscore1"] + V.cSR["bitscore2"]
                            drop_list.append("bitscore1")
                            drop_list.append("bitscore2")

                        # concatenate query group column
                        if "query_group1" in V.cSR.columns:
                            V.cSR["query_group1"].fillna(
                                V.cSR["query_group2"], inplace=True
                            )
                            rename_dict["query_group1"] = "query_group"
                            drop_list.append("query_group2")

                        # concatenate subject group column
                        if "subject_group1" in V.cSR.columns:

                            V.cSR["subject_group1"].fillna(
                                V.cSR["subject_group2"], inplace=True
                            )
                            rename_dict["subject_group1"] = "subject_group"
                            drop_list.append("subject_group2")

                        V.cSR.rename(columns=rename_dict, inplace=True)

                        V.cSR.drop(
                            columns=drop_list,
                            inplace=True,
                        )

        V.cSR = V.cSR

        # Save it
        search.save_df(
            hasher.decode_df(V.dict_hash_name, V.cSR),
            f"{path.out_matrix}/{opt.runname}_BLAST_result_concatenated.{opt.matrixformat}",
            fmt=opt.matrixformat,
        )

    else:
        logging.info("Concatenation not selected, passing concatenation")

    return V

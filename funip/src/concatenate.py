from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import itertools
import logging
from functools import reduce
import pandas as pd
from copy import deepcopy
from funip.src import search, hasher
from scipy.optimize import minimize
import numpy as np
import shutil


# Combine trimmed alignment of each gene to make concatenated matrix, and generate partition file
def combine_alignment(V, opt, path):
    for group in V.dict_dataset:
        if "concatenated" in V.dict_dataset[group]:
            # If more than one locus exists
            if len(V.dict_dataset[group]) > 2:
                # get alignment length
                # length of each of the alignments
                len_dict = {}
                seq_dict = {}
                hash_set = set()
                gene_list = []

                for gene in V.dict_dataset[group]:
                    if not (gene == "concatenated"):
                        gene_list.append(gene)

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
                                for FI in V.dict_dataset[group][
                                    "concatenated"
                                ].list_qr_FI
                            ]
                            + [
                                FI.hash
                                for FI in V.dict_dataset[group][
                                    "concatenated"
                                ].list_db_FI
                            ]
                            + [
                                FI.hash
                                for FI in V.dict_dataset[group][
                                    "concatenated"
                                ].list_og_FI
                            ]
                        )
                        for seq in fasta_list:
                            # if available hash
                            if seq.description in total_dataset:
                                seq_dict[gene][seq.description] = seq
                                hash_set.add(seq.description)

                # Save partition information
                V.partition[group] = {"len": len_dict, "order": gene_list}

                # Generate partition file
                with open(
                    f"{path.out_alignment}/{opt.runname}_{group}.partition", "w"
                ) as fw:
                    tot_len = 0
                    print(gene_list)
                    for gene in gene_list:
                        # gene_list.append(gene)
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
                            # add gaps for ids without gene
                            tmp_seq += "-" * len_dict[gene]

                    concatenate_list.append(
                        SeqRecord(id=hash_id, description="", seq=Seq(tmp_seq))
                    )

                SeqIO.write(
                    concatenate_list,
                    f"{path.out_alignment}/{opt.runname}_trimmed_{group}_concatenated.fasta",
                    "fasta",
                )

            # For single gene, use same alignment from dataset of corresponding gene
            elif len(V.dict_dataset[group]) == 2:
                singlegene = list(set(V.dict_dataset[group].keys()) - {"concatenated"})[
                    0
                ]
                # Copy single gene alignment
                shutil.copy(
                    f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{singlegene}.fasta",
                    f"{path.out_alignment}/{opt.runname}_trimmed_{group}_concatenated.fasta",
                )
                # Generate partition file
                gene_length = len(
                    str(
                        list(
                            SeqIO.parse(
                                f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{singlegene}.fasta",
                                "fasta",
                            )
                        )[0].seq
                    )
                )
                with open(
                    f"{path.out_alignment}/{opt.runname}_{group}.partition", "w"
                ) as fw:
                    fw.write(f"DNA, {singlegene}= 1-{gene_length}\n")

                # Save partition information
                V.partition[group] = {
                    "len": {singlegene: gene_length},
                    "order": [singlegene],
                }

            else:
                logging.error(f"V.dict_dataset {group}: {V.dict_dataset[group]}")
                logging.error(
                    f"[DEVELOPMENTAL ERROR] Failed constructing concatenated alignment for {group}"
                )
                raise Exception

        else:
            logging.warning(
                f"For {group}, no more than genes were detected. Passing combining alignment"
            )

    return V


# for concatenating blast result to get single bitscore among results
def concatenate_df(V, path, opt):
    logging.info("Concatenating search results")

    gene_list = []
    df_list = []
    for gene in V.dict_gene_SR.keys():
        # Leave non-empty dataframes
        if isinstance(V.dict_gene_SR[gene], pd.DataFrame):
            gene_list.append(gene)
            df = deepcopy(V.dict_gene_SR[gene].set_index(["qseqid", "sseqid"]))
            df_list.append(df)

    if len(df_list) <= 0:
        logging.warning(f"Stop concatenating because same or less than 0 gene exists")
        return V

    else:
        # For multigene analysis, some genes can be missing
        # Directly adding those values can make troubles if some of the genes are missing
        # Therefore, add predicted values to fill them with linear regression
        # drop unused columns
        cnt = 0

        # Change biscore names to delimitate genes
        for n, df in enumerate(df_list):
            gene = gene_list[n]
            df[f"{gene}_bitscore"] = df["bitscore"]
            df[f"{gene}_subject_group"] = df["subject_group"]

        # Concatenate dataframes from multiple genes
        df_multigene_regression_ori = pd.concat(df_list, axis=1)

        # Drop unnecessary columns for processing
        df_multigene_regression_ori.drop(
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
                "bitscore",
                "subject_group",
            ],
            inplace=True,
        )

        # Column name managing on subject_group
        def same_merge(x, list_col):
            values = x[list_col].dropna()
            if values.empty:
                raise Exception
            return values.iloc[0]

        df_multigene_regression_ori[
            "subject_group"
        ] = df_multigene_regression_ori.apply(
            lambda x: same_merge(x, [f"{gene}_subject_group" for gene in gene_list]),
            axis=1,
        )

        # For regression, leave anchor points with all genes existing
        df_multigene_regression = df_multigene_regression_ori

        for gene in gene_list:
            df_multigene_regression = df_multigene_regression[
                df_multigene_regression[f"{gene}_bitscore"].notna()
            ]

        # Perform regression
        # Get regression line
        def distance_to_line(line, pts, l0=None, p0=None):
            """
            In three genes situation
            (C0 - X0) / K0 = (C1 -X1) / K1 = (C2 - X2) / K2 = K

            line = (C0, C1, C2)
            pts = [(X0, X1, X2), (Y0, Y1, Y2) ... ]
            p0 = (K0, K1, K2)

            line defined between l0 and line
            points defined between p0 and pts

            This function calcuates following vector distance
            D = (P - (P.dot.u) * u).length
            """
            # line origin other than (0,0,0,..)
            if l0 is not None:
                line = line - l0
            # points origin other than (0,0,0,..)
            if p0 is not None:
                pts = pts - p0

            # dot product
            dp = np.dot(pts, line)
            # dot product value divided by normalized vector of line
            # np.linalg.norm(line) : size of the line vector
            # pp should be orthographic projected length of the dot
            pp = dp / np.linalg.norm(line)
            # norm value of point
            # length from p0 to point
            pn = np.linalg.norm(pts, axis=1)

            return np.sqrt(np.clip(pn**2 - pp**2, a_min=1e-10, a_max=None))

        # Optimization function
        def optimize_regression_line(points):
            n = points.shape[1]  # Dimensionality of the points

            # Define the objective function to minimize (R-squared)
            def objective(x):
                K = x[:n]
                C = x[n:]
                distances = distance_to_line(p0=C, line=K, pts=points)
                mean_squared = np.mean(distances**2)
                logging.debug(f"Mean_squared: {mean_squared}")
                return mean_squared

            # Initial guess for C and K values
            # If C and K are same, it causes initialization error
            x0 = np.full(2 * n, 1)
            x0[:n] = 1

            result = minimize(objective, x0)

            # Extract the optimized C and K values
            C_optimized = result.x[n:]
            K_optimized = result.x[:n]

            return (
                C_optimized,
                K_optimized,
            )

        # Reset df before filling it
        def calculate_prediction(row, gene_list, coeff, grad):
            # Calculate linear_constant of the strain
            linear_constant = []
            for k, gene in enumerate(gene_list):
                if not np.isnan(row[f"{gene}_bitscore"]):
                    #  (coeff - value) / gradient = linear constant
                    linear_constant.append(
                        (coeff[k] - row[f"{gene}_bitscore"]) / grad[k]
                    )

            # Predict unknown BLAST value
            for k, gene in enumerate(gene_list):
                if np.isnan(row[f"{gene}_bitscore"]):
                    #  prediction value  = coeff - linear constant * gradient
                    prediction = coeff[k] - np.mean(linear_constant) * grad[k]
                    row[f"{gene}_bitscore"] = prediction
            return row

        def apply_prediction(row, gene_list, coeff, grad):
            row = calculate_prediction(row, gene_list, coeff, grad)
            return row

        # Change to numpy for faster cazlculation
        np_bitscore = df_multigene_regression[
            [f"{gene}_bitscore" for gene in gene_list]
        ].to_numpy()

        # get coefficient and gradient with regression
        coeff, grad = optimize_regression_line(np_bitscore)

        # Inform users about linear regression result
        # (C0 - X0) / K0 = (C1 -X1) / X1 = (C2 - X2) / X2 = K
        trend_line_string = ""
        for n, gene in enumerate(gene_list):
            trend_line_string += f" {coeff[n]} - {gene} / {grad[n]} ="

        trend_line_string += " K"

        logging.info(
            f"Trend line for Linear regression, \n {trend_line_string} \n Calculated to fill blank genes"
        )

        # fill empty blast results for each gene with regression
        # This might be accelerated by using "apply", but coded manually initially because of logical complexity
        df_multigene_regression = df_multigene_regression_ori.copy()
        df_multigene_regression = df_multigene_regression.apply(
            apply_prediction, args=(gene_list, coeff, grad), axis=1
        )

        # Also update each gene bitscore matrix
        # This part is needed, for multigene analysis, for example ITS, CaM and RPB2
        # If query only exists for ITS and CaM, no blast result for RPB2 were generated
        # So we need to fill it out with linear regression

        # Get summation and save to concatenated search result
        df_multigene_regression["bitscore"] = df_multigene_regression[
            [f"{gene}_bitscore" for gene in gene_list]
        ].mean(axis=1)
        V.cSR = df_multigene_regression.reset_index()

    # Save it
    # decode df is not working well here
    if opt.nosearchresult is False:
        search.save_df(
            hasher.decode_df(hash_dict=V.dict_id_hash, df=V.cSR),
            f"{path.out_matrix}/{opt.runname}_BLAST_result_concatenated.{opt.tableformat}",
            fmt=opt.tableformat,
        )

    return V

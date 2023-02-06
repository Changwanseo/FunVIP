import datetime
import pandas as pd
import matplotlib.pyplot as plt
import io
import logging
import plotly.express as px
from tabulate import tabulate
from funid.src.tool import index_step
from funid.src.save import save_df

# For version reporting
__version__ = "0.3.0"

### Temporary report for tree_interpretation_pipe
# To collect result from multiprocessing on tree_interpretation
# All abnormalities should be listed here
class Singlereport:
    def __init__(self):
        self.id = ""
        self.hash = ""
        self.group = ""
        self.gene = ""
        self.species_original = ""
        self.species_assigned = ""

        # abnormalities
        # for if ambiguous clade exists e.g. Amanita subglobosa "1", Amanita subglobosa "2"
        # 0 if none of them exists
        self.ambiguous = ""
        self.flat = False

    def update_group(self, group):
        self.group = group

    def update_gene(self, gene):
        self.gene = gene

    # Update original species string from (genus, species) tuple
    def update_species_original(self, taxon):
        self.species_original = " ".join(list(taxon))

    # Update assigned species string
    def update_species_assigned(self, taxon_str):
        self.species_assigned = taxon_str

    def __repr__(self):
        return f"{self.hash} {self.gene}"


#############################################################

#############################################################
### Final reporting data
### All reports should be generated from here
class Report:

    # Try to make report in dictionary form, most of them should be printed in tabular form
    def __init__(self):
        # Should be finished after release
        # self.metadata = Metadata()

        # Dataset - group by gene
        # Example)
        #  GROUP        ITS     BenA        CaM        Concatenated
        #  Penicillium   O       O       Fail outgroup      O
        #  Exilicaluis   O     no query      O              O
        # Whether to use just simple O/X or detailed description is not determined

        self.dataset = {"GROUP": []}

        # Identification result statistics : how many of them were well identified
        # Example)
        #

        self.statistics = {
            "GROUP": [],
            "IDENTIFIED": [],
            "AMBIGUOUS": [],
            "NEW SPECIES CANDIDATE": [],
            "MISIDENTIFIED": [],
            "ERROR": [],
            "TOTAL": [],
        }

        # For linking in html, should be added later (on GUI)
        # self.statistics_link = {}

        # By gene results, and sequence should be added
        self.result = {
            "ID": [],
            "HASH": [],
            "DATATYPE": [],
            "GROUP_ORIGINAL": [],
            "GROUP_ASSIGNED": [],
            "SPECIES_ORIGINAL": [],
            "SPECIES_ASSIGNED": [],
            "FLAT_BRANCH": [],
            "INCONSISTENT": [],
            "AMBIGUOUS": [],
            "STATUS": [],
        }

        self.query_result = None

    def update_dataset(self, V, opt):

        for gene in V.list_qr_gene:
            self.dataset[gene.upper()] = []

        self.dataset["CONCATENATE"] = []

        for group in V.list_group:
            self.dataset["GROUP"].append(group)
            if group in V.dict_dataset:
                for gene in V.list_qr_gene:
                    # if group - gene dataset exists
                    if gene in V.dict_dataset[group]:
                        self.dataset[gene.upper()].append("O")
                    # if gene exists but group does not exists
                    else:
                        self.dataset[gene.upper()].append("-")

                if "concatenated" in V.dict_dataset[group]:
                    self.dataset["CONCATENATE"].append("O")
                else:
                    self.dataset["CONCATENATE"].append("-")

            else:
                # If group does not exists
                for gene in V.list_qr_gene:
                    self.dataset[gene.upper()].append("-")

                self.dataset["CONCATENATE"].append("-")

    # Combining result into dictionary form
    def update_result(self, V, opt):
        # get all gene list -> check if V.list_qr_gene is enough for queryonly : False analysis
        set_gene = set(V.list_db_gene + V.list_qr_gene)

        # Generate each gene identification report
        for gene in sorted(list(set_gene)):
            if gene != "concatenated":
                self.result[f"{gene.upper()}_ASSIGNED"] = []

        # Concatenated analysis will be operated seperately
        # set_gene.discard("concatenated") try if without discarding works

        # Collect result from each hash
        for _hash in V.dict_hash_FI:

            FI = V.dict_hash_FI[_hash]

            # Collect result from only used sequences
            if FI.adjusted_group in V.dict_dataset:
                # If any of appropriate gene used
                if any(
                    key in V.dict_dataset[FI.adjusted_group] for key in FI.seq.keys()
                ):

                    # Default results
                    self.result["ID"].append(FI.id)
                    self.result["HASH"].append(FI.hash)
                    self.result["DATATYPE"].append(FI.datatype)

                    if str(FI.group).strip() == "":
                        self.result["GROUP_ORIGINAL"].append("-")
                    else:
                        self.result["GROUP_ORIGINAL"].append(FI.group)

                    self.result["GROUP_ASSIGNED"].append(FI.adjusted_group)

                    if f"{FI.ori_genus} {FI.ori_species}".strip() == "":
                        self.result["SPECIES_ORIGINAL"].append("-")
                    else:
                        self.result["SPECIES_ORIGINAL"].append(
                            f"{FI.ori_genus} {FI.ori_species}"
                        )

                    # Collect identification result for each gene analysis
                    possible_taxon = set()
                    for gene in set_gene:
                        # Check if data analysis had performed for specific FI, group, gene combination
                        if (
                            gene in FI.bygene_species
                            and gene in V.dict_dataset[FI.adjusted_group]
                        ):
                            self.result[f"{gene.upper()}_ASSIGNED"].append(
                                FI.bygene_species[gene]
                            )
                            possible_taxon.add(FI.bygene_species[gene])
                        else:
                            self.result[f"{gene.upper()}_ASSIGNED"].append("-")

                    # Add final identification result
                    if FI.final_species != "":
                        self.result["SPECIES_ASSIGNED"].append(f"{FI.final_species}")
                        possible_taxon.add(FI.final_species)
                    else:
                        self.result["SPECIES_ASSIGNED"].append("UNDETERMINED")

                    # Add abnormalities
                    self.result["FLAT_BRANCH"].append("/".join(FI.flat))
                    if len(possible_taxon) > 1:
                        self.result["INCONSISTENT"].append("inconsistent")
                    else:
                        self.result["INCONSISTENT"].append("-")

                    self.result["AMBIGUOUS"].append(FI.species_identifier)

                    if FI.final_species.strip() == "":
                        self.result["STATUS"].append("ERROR")
                    elif (
                        f"{FI.ori_genus} {FI.ori_species}".strip() == ""
                        and "sp." in FI.final_species.strip()
                    ):
                        self.result["STATUS"].append("new species candidate")
                    elif f"{FI.ori_genus} {FI.ori_species}".strip() == "":
                        self.result["STATUS"].append("undetermined")
                    elif FI.final_species == f"{FI.ori_genus} {FI.ori_species}":
                        self.result["STATUS"].append("correctly identified")
                    else:
                        self.result["STATUS"].append("misidentified")
            else:
                # print(f"DEBUGGING {self.dataset['GROUP']} {FI.adjusted_group}")
                # logging.info(f"{FI} removed from report because not used in analysis")
                pass

        ## update query only result
        self.query_result = pd.DataFrame(self.result)
        # Filter if queryonly is True
        if opt.queryonly is True:
            self.query_result = self.query_result[
                self.query_result["DATATYPE"] == "query"
            ]

        ### Update statistics by result

        # Groupby group
        df_result_group = self.query_result.groupby(["GROUP_ASSIGNED"])

        # Count groups
        for group in sorted(list(set(self.query_result["GROUP_ASSIGNED"]))):
            df_group = df_result_group.get_group(group)

            # Collect statistics
            cnt_correctly_identified = list(df_group["STATUS"]).count(
                "correctly identified"
            )
            cnt_undetermined = list(df_group["STATUS"]).count("undetermined")
            cnt_new_species_candidate = list(df_group["STATUS"]).count(
                "new species candidate"
            )
            cnt_misidentified = list(df_group["STATUS"]).count("misidentified")
            cnt_error = list(df_group["STATUS"]).count("ERROR")
            cnt_total = sum(
                (
                    cnt_correctly_identified,
                    cnt_undetermined,
                    cnt_new_species_candidate,
                    cnt_misidentified,
                    cnt_error,
                )
            )

            # Write into dictionary
            self.statistics["GROUP"].append(group)
            self.statistics["IDENTIFIED"].append(cnt_correctly_identified)
            self.statistics["AMBIGUOUS"].append(cnt_undetermined)
            self.statistics["NEW SPECIES CANDIDATE"].append(cnt_new_species_candidate)
            self.statistics["MISIDENTIFIED"].append(cnt_misidentified)
            self.statistics["ERROR"].append(cnt_error)
            self.statistics["TOTAL"].append(cnt_total)

        # Add final summations
        self.statistics["GROUP"].append("TOTAL")
        self.statistics["IDENTIFIED"].append(sum(self.statistics["IDENTIFIED"]))
        self.statistics["AMBIGUOUS"].append(sum(self.statistics["AMBIGUOUS"]))
        self.statistics["NEW SPECIES CANDIDATE"].append(
            sum(self.statistics["NEW SPECIES CANDIDATE"])
        )
        self.statistics["MISIDENTIFIED"].append(sum(self.statistics["MISIDENTIFIED"]))
        self.statistics["ERROR"].append(sum(self.statistics["ERROR"]))
        self.statistics["TOTAL"].append(sum(self.statistics["TOTAL"]))

    ### Main report runner
    # Update report by pipeline step
    def update_report(self, V, path, opt, step):

        if step == "setup":
            self.report_text(V, path, opt, step)
        elif step == "search":
            self.report_table(V, path, opt, step)
            self.report_text(V, path, opt, step)
        elif step == "cluster":
            self.update_dataset(V, opt)
            self.report_table(V, path, opt, step)
            self.report_text(V, path, opt, step)
            # Datasets are made after clustering
        elif step == "align":
            self.report_text(V, path, opt, step)
        elif step == "trim":
            self.report_text(V, path, opt, step)
        elif step == "concatenate":
            self.report_text(V, path, opt, step)
        elif step == "modeltest":
            self.report_text(V, path, opt, step)
        elif step == "tree":
            self.report_text(V, path, opt, step)
        elif step == "visualize":
            self.report_text(V, path, opt, step)
        elif step == "report":
            self.update_result(V, opt)
            self.report_table(V, path, opt, step)
            self.report_text(V, path, opt, step)
        else:
            logging.error(
                f"DEVELOPMENTAL ERROR : BAD STEP INPUT {step} WHILE UPDATE REPORT"
            )
            raise Exception

    # Generate report in text file
    def report_text(self, V, path, opt, step):

        # Rewrite everytime when called, the io step won't be that much
        with open(f"{path.root}/{opt.runname}.report.txt", "wt", encoding="UTF8") as f:
            if index_step(step) >= 0:
                f.write("FunID Report\n\n")
                ## Write runinfo
                f.write("[INFO]\n")
                # Runname
                f.write(f"RUNNAME      : {opt.runname}\n")
                # Running date
                f.write(
                    f"DATE         : {datetime.datetime.now().strftime('%Y%m%d-%H%M%S')}\n"
                )
                # Running time
                # f.write(f"TIME CONSUMED: {}\n")
                # Number of warnings
                # f.write(f"WARNINGS     : {}\n")
                # Number of errors
                # f.write(f"ERRORS       : {}\n")
                f.write(f"FINAL STEP   : {step}")

            f.write("\n")

            ## Write locations to result files
            if index_step(step) >= 0:
                # DB and Query locations can be written after input step
                f.write(f"* DB files are saved in {path.out_db}\n")
                f.write(f"* Query files are saved in {path.out_query}\n")
                # Options location can be written after option validation step
                f.write(f"* Options are saved in {path.root}/Options.yaml\n")

            if index_step(step) >= 2:
                # Datset table can be written after cluster and outgroup selection step
                f.write(
                    f"* Dataset table can be found in {path.root}/{opt.runname}_Dataset.{opt.matrixformat}\n"
                )
                # Identification table can be written after cluster step
                f.write(
                    f"* Identification result table can be found in {path.root}/{opt.runname}_Identification.{opt.matrixformat}\n"
                )
                f.write(f"* Log files can be found in {path.log}\n")
                f.write("\n")

            ## Write options used (in concise form)
            if index_step(step) >= 0:
                f.write(f"[OPTION]\n")
                f.write(f"DB:                     {opt.query}\n")
                f.write(f"QUERY:                  {opt.db}\n")
                f.write(f"GENE:                   {opt.gene}\n")
                f.write(f"EMAIL:                  {opt.email}\n")
                f.write(f"API:                    {opt.api}\n")
                f.write(f"TEST:                   {opt.test}\n")
                f.write(f"THREAD:                 {opt.thread}\n")
                f.write(f"OUTDIR:                 {opt.outdir}\n")
                f.write(f"RUNNAME:                {opt.runname}\n")
                f.write(f"MODE:                   {opt.mode}\n")
                f.write(f"CONTINUE_FROM_PREVIOUS: {opt.continue_from_previous}\n")
                f.write(f"CRITERION:              {opt.criterion}\n")
                f.write(f"STEP:                   {opt.step}\n")
                f.write(f"LEVEL:                  {opt.level}\n")
                f.write(f"QUERYONLY:              {opt.queryonly}\n")
                f.write(f"CONFIDENT:              {opt.confident}\n")
                f.write(f"VERBOSE:                {opt.verbose}\n")
                f.write(f"MAXOUTGROUP:            {opt.maxoutgroup}\n")
                f.write(f"COLLAPSEDISTCUTOFF:     {opt.collapsedistcutoff}\n")
                f.write(f"COLLAPSEBSCUTOFF:       {opt.collapsebscutoff}\n")
                f.write(f"BOOTSTRAP:              {opt.bootstrap}\n")
                f.write(f"SOLVEFLAT:              {opt.solveflat}\n")
                f.write(f"REGEX:                  {opt.regex}\n")
                f.write(f"AVX:                    {opt.avx}\n")
                f.write(f"CACHEDB:                {opt.cachedb}\n")
                f.write(f"USECACHE:               {opt.usecache}\n")
                f.write(f"MATRIXFORMAT:           {opt.matrixformat}\n")
                f.write(f"NOSEARCHMATRIX:         {opt.nosearchmatrix}\n")
                f.write(f"NOSEARCHRESULT:         {opt.nosearchresult}\n")
                f.write(f"METHOD:                 \n")
                f.write(f" - SEARCH:              {opt.method.search}\n")
                f.write(f" - ALIGNMENT:           {opt.method.alignment}\n")
                f.write(f" - TRIM:                {opt.method.trim}\n")
                f.write(f" - MODELTEST:           {opt.method.modeltest}\n")
                f.write(f" - TREE:                {opt.method.tree}\n")
                f.write(f"VISUALIZE:              \n")
                f.write(f" - BSCUTOFF:            {opt.visualize.bscutoff}\n")
                f.write(f" - FULLGENUS:           {opt.visualize.fullgenus}\n")
                f.write(f" - HIGHLIGHT:           {opt.visualize.highlight}\n")
                f.write(f" - HEIGHTMULTIPLIER:    {opt.visualize.heightmultiplier}\n")
                f.write(f" - MAXWORDLENGTH:       {opt.visualize.maxwordlength}\n")
                f.write(f" - BGCOlOR:             {opt.visualize.backgroundcolor}\n")
                f.write(f" - OUTGROUPCOLOR:       {opt.visualize.outgroupcolor}\n")
                f.write(f" - FTYPE:               {opt.visualize.ftype}\n")
                f.write(f" - FSIZE:               {opt.visualize.fsize}\n")
                f.write(f" - FSIZE_BOOTSTRAP:     {opt.visualize.fsize_bootstrap}\n")
                f.write(f"CLUSTER:                \n")
                f.write(f" - CUTOFF:              {opt.cluster.cutoff}\n")
                f.write(f" - EVALUE:              {opt.cluster.evalue}\n")
                f.write(f" - WORDSIZE:            {opt.cluster.wordsize}\n")
                f.write(f"MAFFT:                  \n")
                f.write(f" - ALGORITHM:           {opt.mafft.algorithm}\n")
                f.write(f" - OP:                  {opt.mafft.op}\n")
                f.write(f" - EP:                  {opt.mafft.ep}\n")
                f.write(f"TRIMAL:                 \n")
                f.write(f" - ALGORITHM:           {opt.trimal.algorithm}\n")
                f.write(f" - GT:                  {opt.trimal.gt}\n")
                f.write("\n")

            ## Write datasets
            # Should be written after clustering
            if index_step(step) >= 2:
                f.write(f"[DATASET]\n")
                f.write(tabulate(self.dataset, headers=self.dataset.keys()))
                f.write("\n\n")
                f.write("O : Group - Gene dataset analysis done\n")
                f.write(
                    "- : Group - Gene dataset analysis cannot be performed (absence of query sequence, no appropriate outgroup etc..)\n"
                )
                f.write("\n\n")

            ## Write identification statistics
            # Should be written after tree_interpretation
            if index_step(step) >= 9:
                f.write(f"[STATISTICS]\n")
                f.write(tabulate(self.statistics, headers=self.statistics.keys()))
                f.write("\n\n")
                f.write(
                    "IDENTIFIED : Number of well-identified strains without any concerns. \n"
                )
                f.write(
                    "AMBIGUOUS : Multiple clades with same taxon name. Your database may contain misidentified sequences. \n"
                )
                f.write(
                    "NEW SPECIES CANDIDATE : New species candidate strains found by topology, phylogenetic distance and bootstrap criteria\n"
                )
                f.write(
                    "MISIDENTIFIED : Strains that shows different identification result from original annotation\n"
                )
                f.write(
                    "ERROR : Strains that cannot be analyzed. Please check if appropriate database sequence / outgroup sequences are given\n"
                )
                f.write("\n\n")

            ## Write identification result
            # Should be written after tree_interpretation
            if index_step(step) >= 9:
                f.write(f"[IDENTIFICATION]\n")
                if opt.queryonly is True:
                    f.write(
                        tabulate(
                            self.query_result,
                            headers=self.query_result.columns,
                            showindex=False,
                        )
                    )
                else:
                    f.write(tabulate(self.result, headers=self.result.keys()))
                f.write("\n\n")
                f.write("ID : Name of the strain\n")
                f.write(
                    "HASH : Temporary name of the strain to prevent unexpected error during run. Use this when manually edit intermediate step data and run from middle, or debugging unexpectively terminated run\n"
                )
                f.write("DATATYPE : query or database\n")
                f.write("GROUP_ORIGINAL : group name given by user\n")
                f.write("GROUP_ASSIGNED : group assigned by FunID clustering\n")
                f.write("SPECIES_ORIGINAL : species name given by user\n")
                f.write(
                    "SPECIES_ASSIGNED : final species name (usually result from concatenated analysis) assigned by FunID tree_interpretation\n"
                )
                f.write(
                    "FLAT_BRANCH : Strains with flat_branch in phylogenetic analysis. If checked, please check your barcode region have enough taxonomic resolution\n"
                )
                f.write(
                    "INCONSISTENT : Strains that show different identification result across genes. Please check sequences were contaminated or misused. \n"
                )
                f.write(
                    "AMBIGUOUS : Multiple clades with same taxon name. Your database may contain misidentified sequences. \n"
                )

                f.write("\n\n")

            ## Write identification methods
            # Can be written after clustering step
            if index_step(step) >= 1:
                f.write(f"[METHOD]\n")
                f.write(f"Sequences were identified with FunID {__version__}\n")
                f.write("\n")

                cnt_db = len([FI for FI in V.list_FI if FI.datatype == "db"])
                cnt_query = len([FI for FI in V.list_FI if FI.datatype == "query"])

                f.write("- Sequence validation -\n")
                f.write(
                    f"Total of {cnt_db} database strains and {cnt_query} strains "
                    f"were used for analysis. "
                )

                f.write(
                    f"Sequences used for analysis were first adjusted to prevent errors during analysis, "
                    f"such as containing invalid bases, non-ascii unicodes, or empty database.\n"
                    # f"During sequences validation step, {cnt_warning} warnings and {cnt_error} errors occured"
                )
                f.write("\n")
                """
                for warning in manage_input_warning:
                    f.write(f"\n")
                for error in manage_input_error:
                    f.write(f"\n")
                """
                f.write(
                    "* Most of the warnings in this step are usually typo problems "
                    "(blanks, tabs, foreign languages that cannot be used in certain programs - like german umlauts) "
                    "and can be automatically fixed by FunID. So you don't have to consider about it that much if you are not going to directly publish.\n"
                )
                f.write("\n")

            if index_step(step) >= 2:
                f.write("- Sequence type identification -\n")
                f.write(
                    f"Sequence type (loci or gene) of query sequences were examined through {opt.method.search} search "
                    f"by using gene of the closeset match. However, ambiguous match with multiple number of genes were revised. "
                    # f"During gene clustering, {cnt_warning} warnings and {cnt_error} errors occured "
                )
                f.write("\n")
                """
                for warning in gene_clustering_warning:
                    f.write(f"\n")
                for error in gene_clustering_error:
                    f.write(f"\n")
                """
                f.write("\n")

                f.write("- Group assignment -\n")
                f.write(
                    f"Sequences were grouped to each datasets by {opt.level} level for phylogenetic analysis through {opt.method.search} search. "
                    f"Therefore, {opt.level} of query sequences are either assigned (when {opt.level} is not given) or validated for clustering "
                    # f"During clustering, {cnt_warning} warnings and {cnt_error}  errors occured "
                )
                f.write("\n")
                """
                for warning in group_clustering_warning:
                    f.write(f"\n")
                for error in group_clustering_error:
                    f.write(f"\n")
                """
                f.write("\n")

                f.write("- Outgroup selection -\n")
                f.write(
                    f"Outgroup sequences were appended to each datasets. "
                    f"{opt.maxoutgroup} sequences closest, but apperantly in distinct group were found by {opt.method.search} search. "
                    # f"During outgroup selection, {cnt_warning} warnings and {cnt_error} errors occured"
                )
                f.write("\n")
                """
                for warning in append_outgroup_warning:
                    f.write(f"\n")
                for error in append_outgroup_error:
                    f.write(f"\n")
                """
                f.write("\n")

                f.write("- Multiple sequence alignment -\n")
                f.write(
                    f"Multiple sequence alignment to each datasets were performed with {opt.method.alignment}. "
                    f"{opt.mafft.algorithm} algorithm was selected with --op {opt.mafft.op} and --ep {opt.mafft.ep} options. "
                    # f"During alignment, {cnt_warning} warnings and {cnt_error} errors occured"
                )
                f.write("\n")
                """
                for warning in mafft_warning:
                    f.write(f"\n")
                for error in mafft_error:
                    f.write(f"\n")
                    """
                f.write("\n")

                f.write("- Alignment trimming -\n")

                if not opt.method.trim == "none":
                    f.write(
                        f"Trimming to each datasets were performed with {opt.method.trim}. "
                        # f"{somewhat} options were used"
                        # f"During trimming, {cnt_warning} warnings and {cnt_error} errors occured"
                    )
                else:
                    f.write(f"No trimming performed on sequence alignments. ")
                f.write("\n")
                """
                for warning in trim_warning:
                    f.write(f"\n")
                for error in trim_error:
                    f.write(f"\n")
                    """
                f.write("\n")

                f.write("- Phylogenetic tree construction -\n")
                f.write(
                    f"Maximum Likelihood tree analysis to each datasets were performed with {opt.method.tree}. "
                    # f"{somewhat} options were used"
                    # f"During trimming, {cnt_warning} warnings and {cnt_error} errors occured"
                )
                f.write("\n")
                """
                for warning in tree_warning:
                    f.write(f"\n")
                for error in tree_error:
                    f.write(f"\n")
                    """
                f.write("\n")

                f.write("- Phylogenetic tree interpretation and identification -\n")
                f.write(
                    f"Tree interpretation were performed to each datasets for species delimitation and tree visualization. "
                    f"Trees were rerooted by outgroup clades. "
                    f"Leaves in ambiguous tree topology, with over {opt.collapsedistcutoff} tree distance "
                )

                if opt.collapsebscutoff < 100:
                    f.write(
                        f"or over {opt.collapsebscutoff} support"
                        # f"During tree interpretation, {cnt_warning} warnings and {cnt_error} errors occured"
                    )

                f.write(
                    f"in common ancestor diverge were considered as distinct species. "
                )

                f.write("\n")
                """
                for warning in tree_interpretation_warning:
                    f.write(f"\n")
                for error in tree_interpretation_error:
                    f.write(f"\n")
                    """
                f.write("\n")

                """
                f.write(
                    f"As a result of FunID analysis, a total of {cnt_query} strains were identified"
                    f"{cnt_query} sequences constists of"
                    f"{cnt_consistent} well identified strains,"
                    f"{cnt_ambiguous} strains that showed different results by genes,"
                    f"{cnt_new_species_candidates} new species candidates,"
                    f"{cnt_error} failed on analysis due to error"
                    f"During tree interpreation, {cnt_warning} warnings and {cnt_error} errors occured"
                )
                for warning in tree_interpretation_warning:
                    f.write(f"\n")
                for error in tree_interpretation_error:
                    f.write(f"\n")
                """

            f.write("\n")
            f.write(
                f"* For precise identification and publications, please double check all warnings and errors\n"
            )
            f.write("\n")

            # Generate software list by options
            software_list = ["FunID"]
            if index_step(step) >= 1:
                if opt.method.search == "blast":
                    software_list.append("BLASTn")
                elif opt.method.search == "mmseqs":
                    software_list.append("MMseqs2")

            if index_step(step) >= 3:
                if opt.method.alignment == "mafft":
                    software_list.append("MAFFT")

            if index_step(step) >= 4:
                if opt.method.trim == "gblocks":
                    software_list.append("Gblocks")
                elif opt.method.trim == "trimal":
                    software_list.append("TrimAl")

            if index_step(step) >= 6:
                if opt.method.modeltest == "modeltest-ng":
                    software_list.append("Modeltest-NG")
                elif opt.method.modeltest == "iqtree":
                    software_list.append("Partitionfinder")

            if index_step(step) >= 7:
                if opt.method.tree == "fasttree":
                    software_list.append("FastTree")
                elif opt.method.tree == "iqtree":
                    software_list.append("IQTREE2")
                elif opt.method.tree == "raxml":
                    software_list.append("RAxML")

            # append software list by step
            dict_version = {
                "FunID": f"{__version__}",
                "BLASTn": "2.13.0",
                "MMseqs2": "14.7e284",
                "MAFFT": "7.453",
                "Gblocks": "0.91b",
                "TrimAl": "1.4",
                "Modeltest-NG": "0.1.7",
                "Partitionfinder": "2.2.0.3",
                "FastTree": "2.1.11",
                "IQTREE2": "2.2.0.3",
                "RAxML": "8.2.12",
            }

            ### Version for each software notation
            f.write(f"[VERSIONS]\n")

            for n, software in enumerate(software_list):
                f.write(f"{'{:<15}'.format(software)} : {dict_version[software]}\n")

            f.write(f"\n")

            ### Write citations
            f.write(f"[CITATION]\n")
            dict_citation = {
                "FunID": "https://github.com/Changwanseo/FunID",
                "BLASTn": "Altschul, S. F., Gish, W., Miller, W., Myers, E. W., & Lipman, D. J. (1990). Basic local alignment search tool. Journal of molecular biology, 215(3), 403-410.",
                "MMseqs2": "Steinegger, M., & Söding, J. (2017). MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature biotechnology, 35(11), 1026-1028.",
                "MAFFT": "Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular biology and evolution, 30(4), 772-780.",
                "Gblocks": "Talavera, G., & Castresana, J. (2007). Improvement of phylogenies after removing divergent and ambiguously aligned blocks from protein sequence alignments. Systematic biology, 56(4), 564-577.",
                "TrimAl": "Capella-Gutiérrez, S., Silla-Martínez, J. M., & Gabaldón, T. (2009). trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. Bioinformatics, 25(15), 1972-1973.",
                "Modeltest-NG": "Darriba, D., Posada, D., Kozlov, A. M., Stamatakis, A., Morel, B., & Flouri, T. (2020). ModelTest-NG: a new and scalable tool for the selection of DNA and protein evolutionary models. Molecular biology and evolution, 37(1), 291-294.",
                "Partitionfinder": "Lanfear, R., Frandsen, P. B., Wright, A. M., Senfeld, T., & Calcott, B. (2017). PartitionFinder 2: new methods for selecting partitioned models of evolution for molecular and morphological phylogenetic analyses. Molecular biology and evolution, 34(3), 772-773.",
                "FastTree": "Price, M. N., Dehal, P. S., & Arkin, A. P. (2010). FastTree 2–approximately maximum-likelihood trees for large alignments. PloS one, 5(3), e9490.",
                "IQTREE2": "Minh, B. Q., Schmidt, H. A., Chernomor, O., Schrempf, D., Woodhams, M. D., Von Haeseler, A., & Lanfear, R. (2020). IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. Molecular biology and evolution, 37(5), 1530-1534.",
                "RAxML": "Stamatakis, A. (2014). RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics, 30(9), 1312-1313.",
            }

            for n, software in enumerate(software_list):
                f.write(
                    f"[{n+1}] {'{:<15}'.format(software)} : {dict_citation[software]}\n"
                )

    def report_table(self, V, path, opt, step):

        # Generate list of table to be written
        list_table = []

        if index_step(step) >= 2:
            list_table.append("dataset")

        if index_step(step) >= 9:
            list_table.append("identification")
            list_table.append("statistics")

        # Write tables
        for table in list_table:
            if table == "dataset":
                save_df(
                    df=pd.DataFrame(self.dataset),
                    out=f"{path.root}/{opt.runname}.dataset.{opt.matrixformat}",
                    fmt=opt.matrixformat,
                )
            elif table == "identification":
                save_df(
                    df=pd.DataFrame(self.result),
                    out=f"{path.root}/{opt.runname}.result.{opt.matrixformat}",
                    fmt=opt.matrixformat,
                )

            elif table == "statistics":
                save_df(
                    df=pd.DataFrame(self.statistics),
                    out=f"{path.root}/{opt.runname}.statistics.{opt.matrixformat}",
                    fmt=opt.matrixformat,
                )
            else:
                logging.error(
                    f"DEVELOPMENTAL ERROR : WRONG TABLE TYPE {table} ATTEMPTED TO BE WRITTEN"
                )
                raise Exception

    def report_html(V, path, opt):
        pass

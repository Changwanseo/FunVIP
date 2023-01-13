from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
import io
import logging
import plotly.express as px


try:
    # from upsetplot import from_indicators, UpSet
    import datapane as dp

    # import altair as alt
except:
    logging.error(
        "Please update your conda environment with following commands:\n Linux: conda env update --file FunID.yaml --prune\n Windows: conda env update --file FunID_Windows.yaml"
    )
    raise Exception


# Temporary report for tree_interpretation_pipe
class Singlereport:
    def __init__(self):
        self.accession = ""
        self.hash = ""
        self.section = ""
        self.gene = ""
        self.inputtaxon = ""
        self.identifiedtaxon = ""
        # for if ambiguous clade exists e.g. Amanita subglobosa 1, Amanita subglobosa 2
        # 0 if none of them exists
        self.taxon_cnt = ""

    def update_genussection(self, string, gene):

        if "concatenated" in string:
            self.section = "_".join(string.split("_")[1:])
        else:
            self.section = string.split("_")[-2]
        self.gene = gene

    def update_inputtaxon(self, taxon):
        self.inputtaxon = " ".join(list(taxon))

    def update_identifiedtaxon(self, taxon):
        self.identifiedtaxon = " ".join(list(taxon))

    def __repr__(self):
        return f"{self.hash} {self.gene}"


#############################################################

#############################################################


# Metadata for the run
class Metadata:
    def __init__(self):
        self.runname = ""
        self.datetime = ""
        self.time_consumption = 0
        self.mode = ""


class Statistics:
    def __init__(self):
        # list of genes designated by option
        self.list_gene = []
        # list of genes found in db
        self.list_gene_db = []
        # list of genes found in query
        self.list_gene_query = []
        # list of genes used for analysis
        self.list_gene_used = []

        # number of total database accessions
        self.cnt_total_db_acc = 0
        # number of used database accessions
        self.cnt_used_db_acc = 0
        # number of downloaded seqs
        self.cnt_db_downloaded_seqs = 0
        # number of failed seqs
        self.cnt_db_failed_seqs = 0

        # number of total query accessions
        self.cnt_query_acc = 0
        # number of query accessions with sections assigned
        self.cnt_query_section_assigned = 0
        # number of query accessions finally assigned to species level
        self.cnt_query_identified = 0

        # number of total sections used for analysis
        self.list_section_used = []

    def update_statistics(self, V, opt):

        self.list_gene = opt.gene
        self.cnt_query_acc = len([FI for FI in V.list_FI if FI.datatype == "Query"])
        self.cnt_query_section_assigned = len(
            [
                FI
                for FI in V.list_FI
                if (FI.datatype == "Query" and FI.adjusted_section != "")
            ]
        )
        self.cnt_query_identified = len(
            [
                FI
                for FI in V.list_FI
                if (FI.datatype == "Query" and FI.final_species != "")
            ]
        )


class Section_Report:
    def __init__(self):
        self.section = ""
        self.genes = []  # genes used for section analysis
        self.gene_statistics = None
        self.upsetplot = None

        # identification statistics
        # number of accessions correctly identified in all possible genes compared to original identifications
        self.cnt_consistent_correctly_identified = 0
        # number of accessions correctly identified in concatenated analysis compared to original identifications
        self.cnt_inconsistent_correctly_identified = 0
        # number of accessions identified to new species candidates
        self.cnt_new_species_candidates = 0
        # number of accessions misidentified
        self.cnt_misidentified = 0
        # number of accessions not included to tree analysis
        self.cnt_errors = 0
        # query identification plot
        self.identificationplot = None
        # clade ambiguity plot
        self.ambiguityplot = None

        self.cnt_ambiguous_species = 0
        self.cnt_unambiguous_species = 0

    def initialize_data(self, V, sect):

        #  for each section
        if sect != "all":

            self.section = sect
            self.genes = list(V.dict_dataset[sect].keys())
        else:

            self.section = sect
            self.genes = list(V.dict_gene_SR.keys()) + ["concatenated"]  # hard coding

    def update_upsetplot(self, V, sect):

        #  for each section
        if sect != "all":

            dataset = V.dict_dataset[sect]

            # gene availaty upsetplot reports
            # initialize upsetplot - columns
            dict_upsetplot = {}
            for gene in self.genes:
                dict_upsetplot[gene] = {}
            dict_upsetplot["Type"] = {}

            # gather all FIs
            set_all_FI = set()
            set_qr_FI = set()
            set_db_FI = set()
            set_og_FI = set()

            for gene in dataset:
                set_qr_FI = set_qr_FI.union(set(dataset[gene].list_qr_FI))
                set_db_FI = set_db_FI.union(set(dataset[gene].list_db_FI))
                set_og_FI = set_og_FI.union(set(dataset[gene].list_og_FI))

            set_all_FI = set_all_FI.union(set_qr_FI)
            set_all_FI = set_all_FI.union(set_db_FI)
            set_all_FI = set_all_FI.union(set_og_FI)

            # initialize upsetplot - index
            for FI in set_all_FI:
                for gene in self.genes:
                    dict_upsetplot[gene][FI.hash] = False

            for FI in set_qr_FI:
                dict_upsetplot["Type"][FI.hash] = "qr"

            for FI in set_db_FI:
                dict_upsetplot["Type"][FI.hash] = "db"

            for FI in set_og_FI:
                dict_upsetplot["Type"][FI.hash] = "og"

            for gene in dataset:
                for FI in dataset[gene].list_qr_FI:
                    dict_upsetplot[gene][FI.hash] = True

                for FI in dataset[gene].list_db_FI:
                    dict_upsetplot[gene][FI.hash] = True

                for FI in dataset[gene].list_og_FI:
                    dict_upsetplot[gene][FI.hash] = True

            df_upsetplot = pd.DataFrame(dict_upsetplot)

            for datatype in ("qr", "db", "og", "all"):
                if datatype == "all":
                    df_upsetplot_input = df_upsetplot
                else:
                    df_upsetplot_input = df_upsetplot[df_upsetplot["Type"] == datatype]

            input_upsetplot = from_indicators(
                indicators=self.genes, data=df_upsetplot_input
            )

            upset = UpSet(input_upsetplot, intersection_plot_elements=0)
            upset.add_stacked_bars(by="Type", elements=10)
            upset.plot()
            plt.title(f"Gene availability for section {sect}")

            tmp = io.BytesIO()
            plt.savefig(tmp, format="svg")
            plt.close()

            self.upsetplot = tmp

        else:

            # gene availaty upsetplot reports
            # initialize upsetplot - columns
            dict_upsetplot = {}
            for gene in self.genes:
                dict_upsetplot[gene] = {}
            dict_upsetplot["Type"] = {}

            for FI in V.list_FI:
                for gene in self.genes:
                    dict_upsetplot[gene][FI.hash] = False

            for FI in V.list_FI:
                if FI.datatype == "Query":
                    dict_upsetplot["Type"][FI.hash] = "qr"
                elif FI.datatype == "DB":
                    dict_upsetplot["Type"][FI.hash] = "db"
                else:
                    raise Exception

                for gene in FI.seq.keys():
                    dict_upsetplot[gene][FI.hash] = True

            df_upsetplot = pd.DataFrame(dict_upsetplot)

            df_upsetplot_input = df_upsetplot

            input_upsetplot = from_indicators(
                indicators=self.genes, data=df_upsetplot_input
            )

            upset = UpSet(input_upsetplot, intersection_plot_elements=0)
            upset.add_stacked_bars(by="Type", elements=10)
            upset.plot()
            plt.title(f"Overall gene availability for inputs")

            tmp = io.BytesIO()
            plt.savefig(tmp, format="svg")
            plt.close()

            self.upsetplot = tmp

            # WIP for all data
            pass

    def update_cntdata(self, V, overall=False):

        # Identification report xlsx file generation
        dict_tmp = {
            "accession": [],
            "hash": [],
            "section": [],
            "original identification": [],
        }

        def counter(self, FI, sect, set_gene):

            # making output excel
            # 0: consistent identification over all genes, 1: not
            inconsistency_flag = 0
            # 1 if correctly identified in final
            correct_flag = 0
            # 1 if new species
            wrong_flag = 0
            newspecies_flag = 0
            # 1 if no concatenated analysis result exists
            error_flag = 0

            if FI.datatype == "Query" and (
                FI.adjusted_section == sect or sect == "all"
            ):
                for gene in set_gene:
                    if gene in FI.bygene_species:
                        if FI.bygene_species[gene] != f"{FI.final_species}":
                            inconsistency_flag = 1
                    else:
                        pass

                # if final species designated by concatenated analysis
                if FI.final_species != "":
                    if "sp." in f"{FI.final_species}":
                        newspecies_flag = 1
                    elif f"{FI.final_species}" == f"{FI.ori_genus} {FI.ori_species}":
                        correct_flag = 1

                # if not final species designated and only one gene analysis exists, take it
                elif len(set_gene) == 1:

                    if "sp." in f"{FI.genus} {FI.bygene_species[list(set_gene)[0]]}":
                        newspecies_flag = 1
                    elif (
                        f"{FI.genus} {FI.bygene_species[list(set_gene)[0]]}"
                        == f"{FI.ori_genus} {FI.ori_species}"
                    ):
                        correct_flag = 1
                else:
                    error_flag = 1

                if error_flag == 1:
                    self.cnt_errors += 1
                elif newspecies_flag == 1:
                    self.cnt_new_species_candidates += 1
                elif correct_flag == 1 and inconsistency_flag == 1:
                    self.cnt_inconsistent_correctly_identified += 1
                elif correct_flag == 1 and inconsistency_flag == 0:
                    self.cnt_consistent_correctly_identified += 1
                else:
                    self.cnt_misidentified += 1

            else:
                pass

        set_gene = set(self.genes)  ## problematic

        #  for each section
        if overall == False:
            sect = self.section
            for h in V.dict_hash_FI:
                FI = V.dict_hash_FI[h]
                counter(self, FI, sect, set_gene)

        # for all sections
        else:
            for FI in V.list_FI:
                counter(self, FI, "all", set_gene)

    def update_ambiguity_plot(self, V, overall=False):

        set_ambiguous = set()
        set_unambiguous = set()

        if overall == False:
            # print(self.section)
            for FI in V.list_FI:
                # if FI is in the right section
                if FI.adjusted_section == self.section or FI.section == self.section:
                    # print(f"Found {FI} for section {self.section}")
                    # FI.final_identification = ""
                    if not (FI.species_identifier == 0):
                        # print("Found non-zero identifier")
                        # debugging
                        # print(FI.final_species, end=" | ")
                        # print(FI.ori_species, end=" | ")
                        # print(FI.species_identifier)
                        # remove from unambiguous list if exists
                        set_unambiguous.discard(FI.final_species)
                        set_ambiguous.add(FI.final_species)
                    else:
                        # print(f"Found zero identifier")
                        if not (FI.final_species in set_ambiguous):
                            # print(f"Final species: {FI.final_species}")
                            set_unambiguous.add(FI.final_species)
                        else:
                            print(f"Exceptional cases: {FI.final_species}")

            # raise Exception

        else:
            pass

        # print(set_ambiguous)
        # print(set_unambiguous)
        self.cnt_ambiguous_species = len(set_ambiguous)
        # print(self.cnt_ambiguous_species)
        self.cnt_unambiguous_species = len(set_unambiguous)
        # print(self.cnt_unambiguous_species)


# logs for each steps
class Singlelog:
    def __init__(self):
        self.info = {}
        self.warning = {}
        self.error = {}


class Log:
    def __init__(self):
        self.log = {}

    def initialize_log(self, list_step):
        for step in list_step:
            self.log[step] = Singlelog()


# Yes, real report
class Report:
    def __init__(self):

        # metadata of the current run
        self.metadata = Metadata()
        # Statistics of the current run
        self.statistics = Statistics()
        # All results - section - report pair dictionary
        self.section_report = {}
        # total report
        self.total_report = Section_Report()
        # Logs - message, warnings, errors
        self.log = Log()

    def initialize_metadata(self, opt):
        self.metadata.runname = opt.runname
        self.metadata.datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        # self.metadata.time_consumption = opt.time_consumption
        self.metadata.mode = opt.mode

        ###### WIP ######

    def initialize_section_report(self, V, opt):

        sections = tuple(V.dict_dataset.keys())

        for sect in sections:
            self.section_report[sect] = Section_Report()
            self.section_report[sect].initialize_data(V, sect)
            if opt.concatenate is True:
                self.section_report[sect].update_upsetplot(V, sect)
            self.section_report[sect].update_cntdata(V)
            self.section_report[sect].update_ambiguity_plot(V)

        self.section_report["all"] = Section_Report()
        self.section_report["all"].initialize_data(V, "all")
        if opt.concatenate is True:
            self.section_report["all"].update_upsetplot(V, "all")
        self.section_report["all"].update_cntdata(V, overall=True)

    # arrange bunch of Singlereport into single dataframe and save
    def arrange(self, V, opt, out):

        # Identification report xlsx file generation
        dict_tmp = {
            "accession": [],
            "hash": [],
            "section": [],
            "original identification": [],
            "final identification": [],
            "compare": [],
        }

        set_gene = set(V.list_db_gene + V.list_qr_gene)

        for gene in sorted(list(set_gene)):
            if gene != "concatenated":
                # dict_tmp["original identification"] = []
                dict_tmp[f"{gene} identified"] = []

        set_gene.discard("concatenated")

        for h in V.dict_hash_FI:  # h for hash
            FI = V.dict_hash_FI[h]
            if opt.queryonly is False or FI.datatype == "Query":
                dict_tmp["accession"].append(FI.accession)
                dict_tmp["hash"].append(FI.hash)
                dict_tmp["section"].append(FI.adjusted_section)
                dict_tmp["original identification"].append(
                    f"{FI.ori_genus} {FI.ori_species}"
                )

                for gene in set_gene:
                    if gene in FI.bygene_species:
                        dict_tmp[f"{gene} identified"].append(FI.bygene_species[gene])
                    else:
                        dict_tmp[f"{gene} identified"].append("-")

                if FI.final_species != "":
                    dict_tmp["final identification"].append(f"{FI.final_species}")

                # elif len(set_gene) == 1:
                #    dict_tmp[f"final identification"].append(
                #        f"{FI.genus} {FI.bygene_species[list(set_gene)[0]]}"
                #    )

                else:
                    dict_tmp[f"final identification"].append("undetermined")

                if FI.final_species.strip() == "":
                    dict_tmp["compare"].append("error")
                elif (
                    f"{FI.ori_genus} {FI.ori_species}".strip() == ""
                    and "sp." in FI.final_species.strip()
                ):
                    dict_tmp["compare"].append("new species candidates")
                elif f"{FI.ori_genus} {FI.ori_species}".strip() == "":
                    dict_tmp["compare"].append("identified")
                elif FI.final_species == f"{FI.ori_genus} {FI.ori_species}":
                    dict_tmp["compare"].append("correctly identified")
                else:
                    dict_tmp["compare"].append("misidentified")

        pd.DataFrame(dict_tmp).to_excel(out, index=None)

    def render(self, opt, out="test.html"):  # render to report

        dp_list = []

        # overall plotting
        dp_list.append(dp.Text("# FunID Report"))
        dp_list.append(
            dp.Text(
                "FunID is a automated fungal identification pipeline with phylogenetic analysis"
            )
        )
        dp_list.append(dp.Text("## Overall statistics"))
        dp_list.append(dp.Text("### Metadata"))
        dp_list.append(dp.Text(f"__Run Name__ : {self.metadata.runname}"))
        dp_list.append(dp.Text(f"__Date__ : {self.metadata.datetime}"))
        dp_list.append(dp.Text(f"__Mode__ : {self.metadata.mode}"))

        dict_df_all = {
            "section": [],
            "consistent_correctly_identified": [],
            "inconsistent_correctly_identified": [],
            "new_species_candidates": [],
            "misidentified": [],
            "errors": [],
        }

        # stacked bar graph of identification status for overall results

        dict_df_all["section"].append("all")
        dict_df_all["consistent_correctly_identified"].append(
            self.section_report["all"].cnt_consistent_correctly_identified
        )
        dict_df_all["inconsistent_correctly_identified"].append(
            self.section_report["all"].cnt_inconsistent_correctly_identified
        )
        dict_df_all["new_species_candidates"].append(
            self.section_report["all"].cnt_new_species_candidates
        )
        dict_df_all["misidentified"].append(
            self.section_report["all"].cnt_misidentified
        )
        dict_df_all["errors"].append(self.section_report["all"].cnt_errors)

        dp_list.append(
            dp.Group(
                dp.BigNumber(
                    heading="Total Query Input", value=self.statistics.cnt_query_acc
                ),
                dp.BigNumber(
                    heading="Total Query Identified",
                    value=self.statistics.cnt_query_identified,
                ),
                dp.Plot(
                    px.bar(
                        pd.DataFrame(dict_df_all),
                        x="section",
                        y=[
                            "consistent_correctly_identified",
                            "inconsistent_correctly_identified",
                            "new_species_candidates",
                            "misidentified",
                            "errors",
                        ],
                    )
                ),
                columns=3,
            )
        )

        # all_plots = []

        if not (self.section_report["all"].identificationplot is None):
            dp_list.append(
                dp.Plot(
                    self.section_report["all"].identificationplot,
                    caption=f"Overall identification result compared to original annotations",
                )
            )

        """
        if not (self.section_report["all"].upsetplot is None):

            all_plots.append(
                dp.HTML(self.section_report["all"].upsetplot.getvalue().decode("utf-8"))
            )

        dp_list.append(
            dp.Group(
                *all_plots,
                columns=2,
            )
        )
        """

        # stacked bar of identification status for each of the section

        list_index = []

        dict_df = {
            "section": [],
            "consistent_correctly_identified": [],
            "inconsistent_correctly_identified": [],
            "new_species_candidates": [],
            "misidentified": [],
            "errors": [],
        }

        for sect in self.section_report:
            if not (sect == "all"):
                dict_df["section"].append(sect)
                dict_df["consistent_correctly_identified"].append(
                    self.section_report[sect].cnt_consistent_correctly_identified
                )
                dict_df["inconsistent_correctly_identified"].append(
                    self.section_report[sect].cnt_inconsistent_correctly_identified
                )
                dict_df["new_species_candidates"].append(
                    self.section_report[sect].cnt_new_species_candidates
                )
                dict_df["misidentified"].append(
                    self.section_report[sect].cnt_misidentified
                )
                dict_df["errors"].append(self.section_report[sect].cnt_errors)

        dp_list.append(
            dp.Plot(
                px.bar(
                    pd.DataFrame(dict_df),
                    x="section",
                    y=[
                        "consistent_correctly_identified",
                        "inconsistent_correctly_identified",
                        "new_species_candidates",
                        "misidentified",
                        "errors",
                    ],
                )
            )
        )

        # same stacked bar graph, for percentage
        dp_list.append(
            dp.Plot(
                px.histogram(
                    pd.DataFrame(dict_df),
                    barnorm="percent",
                    x="section",
                    y=[
                        "consistent_correctly_identified",
                        "inconsistent_correctly_identified",
                        "new_species_candidates",
                        "misidentified",
                        "errors",
                    ],
                )
            )
        )

        # stacked bar graph for ambiguous species
        dict_df = {"section": [], "ambiguous": [], "unambiguous": []}
        for sect in self.section_report:
            if not (sect == "all"):
                dict_df["section"].append(sect)
                dict_df["ambiguous"].append(
                    self.section_report[sect].cnt_ambiguous_species
                )
                dict_df["unambiguous"].append(
                    self.section_report[sect].cnt_unambiguous_species
                )

        dp_list.append(
            dp.Plot(
                px.bar(
                    pd.DataFrame(dict_df),
                    x="section",
                    y=["ambiguous", "unambiguous"],
                )
            )
        )

        dp_list.append(
            dp.Plot(
                px.histogram(
                    pd.DataFrame(dict_df),
                    barnorm="percent",
                    x="section",
                    y=["ambiguous", "unambiguous"],
                )
            )
        )
        """
        dp_list.append(dp.Text("## Sectional statistics"))
        # per section plotting
        section_list = sorted(list(self.section_report.keys()))


        tab_list = []
        for sect in section_list:
            if sect != "all":

                section_plots = []

                # pie chart for section
                if not (self.section_report[sect].identificationplot is None):
                    section_plots.append(
                        dp.Plot(
                            self.section_report[sect].identificationplot,
                            caption=f"Identification result of section {sect} compared to original annotations",
                        )
                    )

                # upset plot for section
                if not (self.section_report[sect].upsetplot is None):

                    section_plots.append(
                        dp.HTML(
                            self.section_report[sect]
                            .upsetplot.getvalue()
                            .decode("utf-8")
                        )
                    )

                # merge all to tab
                # tab_list.append(dp.Text(f"### Statistics for Section {sect}"))
                tab_list.append(
                    dp.Group(
                        dp.Text(f"### Statistics for Section {sect}"),
                        dp.Group(
                            *section_plots,
                            columns=2,
                            label=f"Statistics for Section {sect}",
                        ),
                        columns=1,
                    )
                )

        print(tab_list)
        dp_list.append(dp.Select(blocks=tab_list))
        """

        dp.Report(*tuple(dp_list)).save(path=out)

# FunIP_dev/FunIP/src/command.py

import argparse


class CommandParser:
    def __init__(self) -> None:
        self.parser = argparse.ArgumentParser(
            description="Fungal Identification Pipeline", prog="FunIP"
        )

    def get_args(self) -> argparse.Namespace:
        # Mandatory options
        group_required = self.parser.add_argument_group(
            title="required", description="Main options"
        )
        group_required.add_argument(
            "--query", "-q", nargs="*", help="Query fasta or table files", type=str
        )
        group_required.add_argument(
            "--db",
            "-d",
            nargs="+",
            help="Database table files",
            type=str,
        )
        group_required.add_argument(
            "--gene", "-g", nargs="*", help="Gene names to be analyzed", type=str
        )

        # Mandatory options - ncbi notations
        group_ncbi = self.parser.add_mutually_exclusive_group(required=False)
        group_ncbi.add_argument(
            "--email",
            "-e",
            help="E-mail notation to download data from GenBank",
            type=str,
        )
        group_ncbi.add_argument(
            "--api",
            "-a",
            help="NCBI API strings to download data from GenBank",
            type=str,
        )

        # Test settings
        group_test = self.parser.add_argument_group(
            title="test", description="Test run setups"
        )
        group_test.add_argument(
            "--test",
            help="Use test dataset, [Penicillium]",
            type=str,
        )

        # Running options
        group_run = self.parser.add_argument_group(
            title="run", description="Running options"
        )
        group_run.add_argument(
            "--thread",
            "-t",
            help="Threads to be used for pipeline, default : system maximum",
            type=int,
        )
        group_run.add_argument(
            "--memory",
            "-m",
            help="Max memory limit in 'nG' form, ex: '16G', should be more than 4G, default : system maximum",
            type=str,
        )
        group_run.add_argument(
            "--outdir", help="Out file location, default : current directory", type=str
        )
        group_run.add_argument(
            "--runname",
            help="Name prefix to current run : default : current timestamp",
            type=str,
        )
        group_run.add_argument(
            "--mode",
            help="Mode setup in species identification, see documents for detailed explanations, [validation, identification] default : validation",
            type=str,
        )
        group_run.add_argument(
            "--continue",
            dest="continue_from_previous",
            action="store_true",
            help="Continue from previous run",
        )
        group_run.add_argument(
            "--step",
            help="[WIP] Steps to continue from previous run, will be ignored if invalid --continue option [setup, search, cluster, align, trim, concatenate, modeltest, tree, visualize, report]",
            type=str,
        )
        group_run.add_argument(
            "--level",
            help="Taxonomic level for each phylogenetic tree. Should be one of [subseries, series, subsection, section, subtribe, tribe, subfamily, family, suborder, order, subclass, class, subphylum, phylum, subdivision, division, subkingdom, kingdom]",
            type=str,
        )
        group_run.add_argument(
            "--all",
            action="store_true",
            help="Run FunIP for all database sequences, regardless of corrresponding sequences exists in query, default : False",
        )
        group_run.add_argument(
            "--confident",
            help="Skip blast analysis among database sequences, use it when your database sequences contains large number of misidentified sequences, default : False",
        )

        # Method options
        group_method = self.parser.add_argument_group(
            title="method", description="Methods for each step of pipeline"
        )
        group_method.add_argument(
            "--search",
            help="Search methods to be used in selecting genes, groups and outgroups, [blast, mmseqs], default : mmseqs",
            type=str,
        )
        group_method.add_argument(
            "--alignment",
            help="Multiple sequence alignment methods, [mafft], default : mafft",
            type=str,
        )
        group_method.add_argument(
            "--trim",
            help="Trimming methods, [trimal, gblocks, none], default : trimal",
            type=str,
        )
        group_method.add_argument(
            "--modeltest",
            help="Model test methods, [iqtree, modeltestng, none], default : none",
            type=str,
        )
        group_method.add_argument(
            "--tree",
            help="Tree methods to build phylogenetic tree, [fasttree, iqtree, raxml], default : fasttree",
            type=str,
        )

        # Visualize
        group_visualize = self.parser.add_argument_group(
            title="visualize",
            description="Visualization options for drawing phylogenetic tree",
        )
        group_visualize.add_argument(
            "--bscutoff",
            help="Bootstrap cutoff for visualize, default : 70",
            type=int,
        )
        group_visualize.add_argument(
            "--highlight",
            help="Color to highlight query sequences in tree visualization. Either in html svg recognizable string or hex code, default: #AA0000",
            type=str,
        )
        group_visualize.add_argument(
            "--heightmultiplier",
            help="Height multiplier in drawing collapsing nodes. Change it if you want to show collapse node more or less expanded. Default: 6",
            type=float,
        )
        group_visualize.add_argument(
            "--maxwordlength",
            help="Maximum letters to be shown in single line of tree annotation. Default: 48",
            type=int,
        )

        group_visualize.add_argument(
            "--backgroundcolor",
            help='List of background colors to be shown in tree, default: #f4f4f4, #c6c6c6. Input should be used with quotes, delimit with spaces and recommended to be used as hex codes. To remove background, use --backgroundcolor "#FFFFFF" "#FFFFFF" ',
            nargs="*",
            type=str,
        )
        group_visualize.add_argument(
            "--outgroupcolor",
            help="Background colors to indicate outgroup, default: #999999",
            type=str,
        )

        group_visualize.add_argument(
            "--ftype",
            help="Font to use for phylogenetic tree, default: Arial",
            type=str,
        )
        group_visualize.add_argument(
            "--fsize",
            help="Font size to use for phylogenetic tree, default: 10",
            type=float,
        )
        group_visualize.add_argument(
            "--fsize_bootstrap",
            help="Font size to use for bootstrap support in phylogenetic tree, default: 9",
            type=float,
        )

        # Advanced
        group_advanced = self.parser.add_argument_group(
            title="advanced", description="Advanced options for minor controls"
        )
        group_advanced.add_argument(
            "--verbose",
            "-v",
            help="Verbosity level, 0: quiet, 1: info, 2: warning, 3: debug, default : 2",
            type=int,
        )
        group_run.add_argument(
            "--maxoutgroup",
            help="Maximum outgroup numbers to include in phylogenetic analysis, default : 1",
        )
        group_advanced.add_argument(
            "--collapsedistcutoff",
            help="Maximum tree distance to be considered as same species, default : 0.01",
            type=float,
        )
        group_advanced.add_argument(
            "--collapsebscutoff",
            help="Minimum bootstrap to be considered as same species, default : 100",
            type=float,
        )
        group_advanced.add_argument(
            "--bootstrap",
            help="Boostrap number for tree analysis, will be ignored if fasttree is selected for tree method, default : 1000",
            type=int,
        )
        group_advanced.add_argument(
            "--solveflat",
            action="store_true",
            help="Whether to automatically detect 0 length branch and automatically solve them, default : True",
        )
        group_advanced.add_argument(
            "--regex",
            nargs="*",
            help="Regex groups to parse strain numbers from your input. Maybe useful if your sequence descriptions are dirty. See documentation",
            type=str,
        )

        group_advanced.add_argument(
            "--cluster-cutoff",
            dest="cluster_cutoff",
            help="Minimum percent identity to be considered as same group in clustering analysis. Should be between 0 and 1, default : 0.97",
            type=float,
        )
        group_advanced.add_argument(
            "--cluster-evalue",
            dest="cluster_evalue",
            help="E-value cutoffs for blast/mmseqs search, default : 0.0000001",
            type=float,
        )
        group_advanced.add_argument(
            "--cluster-wordsize",
            dest="cluster_wordsize",
            help="Word size for blast/mmseqs search, default : 7",
            type=int,
        )
        group_advanced.add_argument(
            "--mafft-algorithm",
            dest="mafft_algorithm",
            help="MAFFT algorithm for alignment, see mafft documents, will be ignored if mafft not selected for alignment option, default : auto",
        )
        group_advanced.add_argument(
            "--mafft-op",
            dest="mafft_op",
            help="MAFFT op (gap opening penalty) value, default : 1.3",
            type=float,
        )
        group_advanced.add_argument(
            "--mafft-ep",
            dest="mafft_ep",
            help="MAFFT ep value, default : 0.1",
            type=float,
        )
        group_advanced.add_argument(
            "--trimal-algorithm",
            dest="trimal_algorithm",
            help="Trimal algorithm for trimming, see trimal documents, will be ignored if trimal not selected for trimming option, default : gt",
        )
        group_advanced.add_argument(
            "--trimal-gt", dest="trimal_gt", help="gt value for trimal", type=float
        )
        group_advanced.add_argument(
            "--allow-innertrimming",
            dest="allow_innertrimming",
            help="Turn off FunIP adjustment to not to trim inner alignment columns",
            action="store_true",
        )
        group_advanced.add_argument(
            "--criterion",
            help="Modeltest criterion to use, either AIC, AICc or BIC",
            type=str,
        )

        group_advanced.add_argument(
            "--noavx",
            action="store_true",
            help="do not use AVX for RAxML, default: False",
        )
        group_advanced.add_argument(
            "--outgroupoffset",
            help="outgroupoffset value. Highering this value may select more distant outgroup, default : 20",
            type=int,
        )

        # Cache
        group_cache = self.parser.add_argument_group(
            title="cache",
            description="Save search database for faster run in next time",
        )
        group_cache.add_argument(
            "--cachedb",
            action="store_true",
            help="Cache current search database, turn off if your database is too big for system directory, default : True",
        )
        group_cache.add_argument(
            "--usecache",
            action="store_true",
            help="Use cached search database, turn off if your cached database makes error, default : True",
        )

        # Save
        group_save = self.parser.add_argument_group(
            title="save", description="Run saving options"
        )
        group_save.add_argument(
            "--matrixformat",
            help="Default format for search matrix files, [csv, xlsx, parquet, feather], default : csv",
        )
        group_save.add_argument(
            "--nosearchresult",
            action="store_true",
            help="Do not save blast/mmseqs search matrix, use when dataset gets too big and generates IO bottleneck",
        )

        # Preset
        group_setting = self.parser.add_argument_group(
            title="setting", description="Presets for one-step settings"
        )
        group_setting.add_argument(
            "--preset",
            help="[fast, accurate], or json formatted option config file. Check documentation for each preset, default : fast",
            type=str,
        )

        return self.parser.parse_args()

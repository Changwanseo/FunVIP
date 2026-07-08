# FunVIP_dev/FunVIP/src/command.py

import argparse
from importlib.metadata import version


class CommandParser:
    def __init__(self) -> None:
        self.parser = argparse.ArgumentParser(
            description="Fungal Validation & Identification Pipeline", prog="FunVIP"
        )

    def get_args(self) -> argparse.Namespace:
        # Mandatory options
        group_required = self.parser.add_argument_group(
            title="required", description="Main options"
        )
        group_required.add_argument(
            "--query",
            "-q",
            nargs="*",
            help="Query fasta or table files, delimitated by space",
            type=str,
        )
        group_required.add_argument(
            "--db",
            "-d",
            nargs="+",
            help="Database table files, delimitated by space",
            type=str,
        )
        group_required.add_argument(
            "--gene", "-g", nargs="*", help="Gene names to be analyzed", type=str
        )

        # Mandatory when NCBI download - raise Exception when manage input
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
            help="Run on a bundled test dataset, one of [penicillium, terrei, sanghuangporus]. Requires --email for the bundled GenBank accessions.",
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
            help="Analysis mode, one of [identification, validation]. default : identification",
            type=str,
        )
        group_run.add_argument(
            "--continue",
            dest="continue_from_previous",
            action="store_true",
            help="Continue from previous run, default: False",
        )
        group_run.add_argument(
            "--step",
            help="With --continue, resume from this pipeline step, one of [setup, search, cluster, align, trim, concatenate, modeltest, tree, visualize, report]. Ignored without a valid --continue.",
            type=str,
        )
        group_run.add_argument(
            "--level",
            help="Taxonomic level at which each phylogenetic tree is built, one of [genus, subseries, series, subsection, section, subtribe, tribe, subfamily, family, suborder, order, subclass, class, subphylum, phylum, subdivision, division, subkingdom, kingdom]. default : genus",
            type=str,
        )
        group_run.add_argument(
            "--all",
            action="store_true",
            help="Run FunVIP for all database sequences, regardless of corresponding sequences exists in query, default : False",
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
            help="Search method for selecting genes, groups and outgroups, one of [blast, mmseqs]. default : blast",
            type=str,
        )
        group_method.add_argument(
            "--alignment",
            help="Multiple sequence alignment methods, [mafft], default : mafft",
            type=str,
        )
        group_method.add_argument(
            "--notcs",
            action="store_true",
            help="Skip T-COFFEE TCS(Transitive Consistency Score) for alignment validation. default: False",
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
            help='Alternating background band colors in the tree, space-delimited hex codes in quotes. default: "#ffe0e0" "#ffefef". To remove the background use --backgroundcolor "#FFFFFF" "#FFFFFF".',
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
            help="Font size for phylogenetic tree labels, default: 14",
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
            help="Maximum number of outgroup sequences to include in each phylogenetic tree, default : 3",
            type=int,
        )
        group_advanced.add_argument(
            "--collapsedistcutoff",
            help="Maximum tree distance to be considered as same species, default : 0.01",
            type=float,
        )
        group_advanced.add_argument(
            "--collapsebscutoff",
            help="Minimum bootstrap support to collapse a clade as one species; the default 101 (above the 100 maximum) disables bootstrap-based collapsing. default : 101",
            type=float,
        )
        group_advanced.add_argument(
            "--bootstrap",
            help="Bootstrap replicates for tree analysis (ignored when tree method is fasttree). default : 100 (the accurate preset uses 1000)",
            type=int,
        )
        group_advanced.add_argument(
            "--nosolveflat",
            action="store_true",
            help="Do not detect 0 length branch and automatically solve them, default : False",
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
            help="E-value cutoff for the blast/mmseqs clustering search; an explicit value overrides the preset. default : 0.0001",
            type=float,
        )
        group_advanced.add_argument(
            "--cluster-wordsize",
            dest="cluster_wordsize",
            help="Word size for blast/mmseqs search, default : 7",
            type=int,
        )

        group_advanced.add_argument(
            "--cluster-max_target_seqs",
            dest="cluster_max_target_seqs",
            help="Max_target_seqs for blast/mmseqs search, increase this when expected search match is not found. default : 100",
            type=int,
        )
        group_advanced.add_argument(
            "--mafft-algorithm",
            dest="mafft_algorithm",
            help="MAFFT algorithm for alignment (see MAFFT docs). default : auto",
            type=str,
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
            help="trimAl algorithm for trimming (see trimAl docs). default : gt",
            type=str,
        )
        group_advanced.add_argument(
            "--trimal-gt", dest="trimal_gt", help="trimAl -gt (gap threshold), between 0 and 1. default : 0.2", type=float
        )
        group_advanced.add_argument(
            "--allow-innertrimming",
            dest="allow_innertrimming",
            help="Turn off FunVIP adjustment to not to trim inner alignment columns, default: False",
            action="store_true",
        )
        group_advanced.add_argument(
            "--criterion",
            help="Model-selection criterion for modeltest, one of [AIC, AICc, BIC]. default : BIC",
            type=str,
        )

        group_advanced.add_argument(
            "--noavx",
            action="store_true",
            help="Do not use AVX for RAxML, default: False",
        )
        group_advanced.add_argument(
            "--outgroupoffset",
            help="outgroupoffset value. Highering this value may select more distant outgroup, default : 20",
            type=int,
        )
        group_advanced.add_argument(
            "--nosuspicious",
            action="store_true",
            help="Do not include suspicious samples for sequence-set. Mostly for metabarcoding analysis. May deduce inaccurate result with problematic database.",
        )
        group_advanced.add_argument(
            "--terminate",
            action="store_true",
            help="Terminate FunVIP run when critical error detected",
        )

        # Cache
        group_cache = self.parser.add_argument_group(
            title="cache",
            description="Save search database for faster run in next time",
        )
        group_cache.add_argument(
            "--nocachedb",
            action="store_true",
            help="Disable caching current search database. Use it if your database is too big for system directory, default : True",
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
        # Compatibility test required for this one
        """
        group_save.add_argument(
            "--tableformat",
            help="Default format for table files, [csv, xlsx, parquet, feather], default : csv",
        )
        """
        group_save.add_argument(
            "--tableformat",
            help="Default format for table files, [csv, xlsx], default : csv",
        )

        group_save.add_argument(
            "--nosearchresult",
            action="store_true",
            help="Do not save blast/mmseqs search matrix, use when dataset gets too big and generates IO bottleneck, default: False",
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

        # version
        self.parser.add_argument(
            "--version", action="version", version=f"FunVIP {version('FunVIP')}"
        )

        return self.parser.parse_args()

FunID documentation

usage: FunID [-h] [--query [QUERY ...]] [--db DB [DB ...]] [--gene [GENE ...]] [--email EMAIL | --api API] [--test TEST] [--thread TH
             [--memory MEMORY] [--outdir OUTDIR] [--runname RUNNAME] [--mode MODE] [--continue] [--step STEP] [--level LEVEL] [--all]
             [--confident CONFIDENT] [--search SEARCH] [--alignment ALIGNMENT] [--trim TRIM] [--modeltest MODELTEST] [--tree TREE] [-
             [--highlight HIGHLIGHT] [--heightmultiplier HEIGHTMULTIPLIER] [--maxwordlength MAXWORDLENGTH] [--backgroundcolor [BACKGR
             [--outgroupcolor OUTGROUPCOLOR] [--ftype FTYPE] [--fsize FSIZE] [--fsize_bootstrap FSIZE_BOOTSTRAP] [--verbose VERBOSE]
             [--maxoutgroup MAXOUTGROUP] [--collapsedistcutoff COLLAPSEDISTCUTOFF] [--collapsebscutoff COLLAPSEBSCUTOFF] [--bootstrap
             [--solveflat] [--regex [REGEX ...]] [--cluster-cutoff CLUSTER_CUTOFF] [--cluster-evalue CLUSTER_EVALUE]
             [--cluster-wordsize CLUSTER_WORDSIZE] [--mafft-algorithm MAFFT_ALGORITHM] [--mafft-op MAFFT_OP] [--mafft-ep MAFFT_EP]
             [--trimal-algorithm TRIMAL_ALGORITHM] [--trimal-gt TRIMAL_GT] [--allow-innertrimming] [--criterion CRITERION] [--noavx]
             [--outgroupoffset OUTGROUPOFFSET] [--cachedb] [--usecache] [--matrixformat MATRIXFORMAT] [--nosearchresult] [--preset PR

Fungal Identification Pipeline

options:
  -h, --help            show this help message and exit
  --email EMAIL, -e EMAIL
                        E-mail notation to download data from GenBank
                        (Without email, NCBI may ban your IP. We will remove requirement of email when connection to GenBank is not needed, but currently, email is mandatory)


  --api API, -a API     Dummy option - NCBI API strings to download data from GenBank
                        (Will be developed for faster download from GenBank, but currently may not work)

  Main options

  --query [QUERY ...], -q [QUERY ...]
                        Query fasta or table files. Data that you want to analyze. You can either write file names (if they are in current directory), or file path (either full path or relative path is okay). See figure (LINK) to how to generate query files.

  --db DB [DB ...], -d DB [DB ...] [REQUIRED!!!!]
                        Database table files. You can either write file names (if they are in current directory), or file path (either full path or relative path is okay). See figure (LINK) to how to generate database files.

  --gene [GENE ...], -g [GENE ...] [REQUIRED!!!!]
                        Gene names to be analyzed. Should be same to corresponding column name (Capitzliation differences are okay)

test:
  Test run setups

  --test TEST           Use test dataset, [Penicillium] or [Terrei]
                        Penicillium dataset compromises data from GenMine paper. It will take 10+ minutes in 12c system.
                        
                        Terrei dataset compromises data from FunID paper. It will take 3+ minutes in 12c system. This test requires additional email

run:
  Running options

  --thread THREAD, -t THREAD
                        Threads to be used for pipeline, default : system maximum
                        If FunID uses too much memory, decreasing thread may help

  --memory MEMORY, -m MEMORY
                        Max memory limit in 'nG' form, ex: '16G', should be more than 4G, default : system maximum
                        This limit only applies for mmseqs and iqtree. Therefore, if FunID exceeds memory in other steps, please decrease thread number

  --outdir OUTDIR       Out file location, default : current directory
                        If you run FunID with --continue, this is mandatory

  --runname RUNNAME     Name prefix to current run : default : current timestamp
                        If you run FunID with --continue, this is mandatory

  --mode MODE           Mode setup in species identification, see documents for detailed explanations, [validation, identification]
                        The mode decides species assignment of zero-length polytomy branch, if species name of queries are given.
                        For example, two database sample, Penicillium A, and Penicillium B shows zero-length polytomy branch.
                        In validation mode, if query is given as Penicillium A, FunID decides query as Penicillium A (Believes query information)
                        In identification mode, if query is given as Penicillium A, FunID decides query as Penicillium sp. 1 (Cannot completely decide if Penicillium A or Penicillium B are right)
                    
  --continue            Continue from previous run. Should be used with --outdir, --runname, --step to designate previous run directory and which steps to be continued

  --step STEP           Steps to continue from previous run, will be ignored without --continue option [setup, search, clust
                        concatenate, modeltest, tree, visualize, report]

  --level LEVEL         Taxonomic level for each phylogenetic tree. Should be one of [subseries, series, subsection, section, subtrib
                        subfamily, family, suborder, order, subclass, class, subphylum, phylum, subdivision, division, subkingdom, kingdom

  --all                 Run FunID for all database sequences, regardless of corrresponding sequences exists in query, default : False
                        If you don't have any query sequences, --all flag will be automatically on.

  --confident CONFIDENT
                        Skip blast analysis among database sequences, use it when your database sequences contains large number of mi
                        sequences, and confirm with higher level assignment default : False

  --maxoutgroup MAXOUTGROUP
                        Maximum outgroup numbers to include in phylogenetic analysis, default : 3

method:
  Methods for each step of pipeline

  --search SEARCH       Search methods to be used in selecting genes, groups and outgroups, [blast, mmseqs], default : mmseqs
                        BLAST is more sensitive, while mmseqs is faster in large number of sequences. MMseqs uses large space of memory, so more than 16GB or RAM is recommended

  --alignment ALIGNMENT
                        Multiple sequence alignment methods, [mafft], default : mafft

  --trim TRIM           Trimming methods, [trimal, gblocks, none], default : trimal
                        Default trimming methods are curated to not to remove internal columns. If you want to use original behavior of trimal or gblocks, use with --allow-innertrimming

  --modeltest MODELTEST
                        Model test methods, [iqtree, modeltestng, none], default : none
                        Modeltest-ng is currently not available with Windows, but will be supported soon!

  --tree TREE           Tree methods to build phylogenetic tree, [fasttree, iqtree, raxml], default : fasttree

visualize:
  Visualization options for drawing phylogenetic tree

  --bscutoff BSCUTOFF   Bootstrap cutoff for visualize, default : 70
                        If 70, bootstrap values which are equals or above 70 are represented in tree

  --highlight HIGHLIGHT
                        Color to highlight query sequences in tree visualization. Either in html svg recognizable string or hex code, 
                        i.e. red, blue, magenta, #03030A

  --heightmultiplier HEIGHTMULTIPLIER
                        Height multiplier in drawing collapsing nodes. Change it if you want to show collapse node more or less expan
                        
  --maxwordlength MAXWORDLENGTH
                        Maximum letters to be shown in single line of tree annotation. Default: 48
  --backgroundcolor [BACKGROUNDCOLOR ...]
                        List of background colors to be shown in tree, default: #f4f4f4, #c6c6c6. Input should be used with quotes, d
                        and recommended to be used as hex codes. To remove background, use --backgroundcolor "#FFFFFF" "#FFFFFF"
  --outgroupcolor OUTGROUPCOLOR
                        Background colors to indicate outgroup, default: #999999
  --ftype FTYPE         Font to use for phylogenetic tree, default: Arial
  --fsize FSIZE         Font size to use for phylogenetic tree, default: 10
  --fsize_bootstrap FSIZE_BOOTSTRAP
                        Font size to use for bootstrap support in phylogenetic tree, default: 9

advanced:
  Advanced options for minor controls

  --verbose VERBOSE, -v VERBOSE
                        Verbosity level, 0: quiet, 1: info, 2: warning, 3: debug, default : 2
  --collapsedistcutoff COLLAPSEDISTCUTOFF
                        Maximum tree distance to be considered as same species, default : 0.01
  --collapsebscutoff COLLAPSEBSCUTOFF
                        Minimum bootstrap to be considered as same species, default : 100
  --bootstrap BOOTSTRAP
                        Boostrap number for tree analysis, will be ignored if fasttree is selected for tree method, default : 1000
  --solveflat           Whether to automatically detect 0 length branch and automatically solve them, default : True
  --regex [REGEX ...]   Regex groups to parse strain numbers from your input. Maybe useful if your sequence descriptions are dirty. S
  --cluster-cutoff CLUSTER_CUTOFF
                        Minimum percent identity to be considered as same group in clustering analysis. Should be between 0 and 1, de
  --cluster-evalue CLUSTER_EVALUE
                        E-value cutoffs for blast/mmseqs search, default : 0.0000001
  --cluster-wordsize CLUSTER_WORDSIZE
                        Word size for blast/mmseqs search, default : 7
  --mafft-algorithm MAFFT_ALGORITHM
                        MAFFT algorithm for alignment, see mafft documents, will be ignored if mafft not selected for alignment optio
  --mafft-op MAFFT_OP   MAFFT op (gap opening penalty) value, default : 1.3
  --mafft-ep MAFFT_EP   MAFFT ep value, default : 0.1
  --trimal-algorithm TRIMAL_ALGORITHM
                        Trimal algorithm for trimming, see trimal documents, will be ignored if trimal not selected for trimming opti
  --trimal-gt TRIMAL_GT
                        gt value for trimal
  --allow-innertrimming
                        Turn off FunID adjustment to not to trim inner alignment columns
  --criterion CRITERION
                        Modeltest criterion to use, either AIC, AICc or BIC
  --noavx               do not use AVX for RAxML, default: False
  --outgroupoffset OUTGROUPOFFSET
                        outgroupoffset value. Highering this value may select more distant outgroup, default : 20

cache:
  Save search database for faster run in next time

  --cachedb             Cache current search database, turn off if your database is too big for system directory, default : True
  --usecache            Use cached search database, turn off if your cached database makes error, default : True

save:
  Run saving options

  --matrixformat MATRIXFORMAT
                        Default format for search matrix files, [csv, xlsx, parquet, feather], default : csv
  --nosearchresult      Do not save blast/mmseqs search matrix, use when dataset gets too big and generates IO bottleneck

setting:
  Presets for one-step settings

  --preset PRESET       [fast, accurate], or json formatted option config file. Check documentation for each preset, default : fast
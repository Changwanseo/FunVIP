# FunVIP command-line parameters

_Auto-generated from `funvip/src/command.py` by `docs/gen_param_reference.py`;
do not edit by hand. Regenerate after changing the CLI._

## options

| Option | Type | Description |
|---|---|---|
| `-h`, `--help` | flag | show this help message and exit |
| `--email`, `-e` | str | E-mail notation to download data from GenBank |
| `--api`, `-a` | str | NCBI API strings to download data from GenBank |
| `--version` | flag | show program's version number and exit |

## required

| Option | Type | Description |
|---|---|---|
| `--query`, `-q` | str | Query fasta or table files, delimitated by space |
| `--db`, `-d` | str | Database table files, delimitated by space |
| `--gene`, `-g` | str | Gene names to be analyzed |

## test

| Option | Type | Description |
|---|---|---|
| `--test` | str | Run on a bundled test dataset, one of [penicillium, terrei, sanghuangporus]. Requires --email for the bundled GenBank accessions. |

## run

| Option | Type | Description |
|---|---|---|
| `--thread`, `-t` | int | Threads to be used for pipeline, default : system maximum |
| `--memory`, `-m` | str | Max memory limit in 'nG' form, ex: '16G', should be more than 4G, default : system maximum |
| `--outdir` | str | Out file location, default : current directory |
| `--runname` | str | Name prefix to current run : default : current timestamp |
| `--mode` | str | Analysis mode, one of [identification, validation]. default : identification |
| `--continue` | flag | Continue from previous run, default: False |
| `--step` | str | With --continue, resume from this pipeline step, one of [setup, search, cluster, align, trim, concatenate, modeltest, tree, visualize, report]. Ignored without a valid --continue. |
| `--level` | str | Taxonomic level at which each phylogenetic tree is built, one of [genus, subseries, series, subsection, section, subtribe, tribe, subfamily, family, suborder, order, subclass, class, subphylum, phylum, subdivision, division, subkingdom, kingdom]. default : genus |
| `--all` | flag | Run FunVIP for all database sequences, regardless of corresponding sequences exists in query, default : False |
| `--confident` | str | Skip blast analysis among database sequences, use it when your database sequences contains large number of misidentified sequences, default : False |
| `--maxoutgroup` | int | Maximum number of outgroup sequences to include in each phylogenetic tree, default : 3 |

## method

| Option | Type | Description |
|---|---|---|
| `--search` | str | Search method for selecting genes, groups and outgroups, one of [blast, mmseqs]. default : blast |
| `--alignment` | str | Multiple sequence alignment methods, [mafft], default : mafft |
| `--notcs` | flag | Skip T-COFFEE TCS(Transitive Consistency Score) for alignment validation. default: False |
| `--trim` | str | Trimming methods, [trimal, gblocks, none], default : trimal |
| `--modeltest` | str | Model test methods, [iqtree, modeltestng, none], default : none |
| `--tree` | str | Tree methods to build phylogenetic tree, [fasttree, iqtree, raxml], default : fasttree |

## visualize

| Option | Type | Description |
|---|---|---|
| `--bscutoff` | int | Bootstrap cutoff for visualize, default : 70 |
| `--highlight` | str | Color to highlight query sequences in tree visualization. Either in html svg recognizable string or hex code, default: #AA0000 |
| `--heightmultiplier` | float | Height multiplier in drawing collapsing nodes. Change it if you want to show collapse node more or less expanded. Default: 6 |
| `--maxwordlength` | int | Maximum letters to be shown in single line of tree annotation. Default: 48 |
| `--backgroundcolor` | str | Alternating background band colors in the tree, space-delimited hex codes in quotes. default: "#ffe0e0" "#ffefef". To remove the background use --backgroundcolor "#FFFFFF" "#FFFFFF". |
| `--outgroupcolor` | str | Background colors to indicate outgroup, default: #999999 |
| `--ftype` | str | Font to use for phylogenetic tree, default: Arial |
| `--fsize` | float | Font size for phylogenetic tree labels, default: 14 |
| `--fsize_bootstrap` | float | Font size to use for bootstrap support in phylogenetic tree, default: 9 |

## advanced

| Option | Type | Description |
|---|---|---|
| `--verbose`, `-v` | int | Verbosity level, 0: quiet, 1: info, 2: warning, 3: debug, default : 2 |
| `--collapsedistcutoff` | float | Maximum tree distance to be considered as same species, default : 0.01 |
| `--collapsebscutoff` | float | Minimum bootstrap support to collapse a clade as one species; the default 101 (above the 100 maximum) disables bootstrap-based collapsing. default : 101 |
| `--bootstrap` | int | Bootstrap replicates for tree analysis (ignored when tree method is fasttree). default : 100 (the accurate preset uses 1000) |
| `--nosolveflat` | flag | Do not detect 0 length branch and automatically solve them, default : False |
| `--regex` | str | Regex groups to parse strain numbers from your input. Maybe useful if your sequence descriptions are dirty. See documentation |
| `--cluster-cutoff` | float | Minimum percent identity to be considered as same group in clustering analysis. Should be between 0 and 1, default : 0.97 |
| `--cluster-evalue` | float | E-value cutoff for the blast/mmseqs clustering search; an explicit value overrides the preset. default : 0.0001 |
| `--cluster-wordsize` | int | Word size for blast/mmseqs search, default : 7 |
| `--cluster-max_target_seqs` | int | Max_target_seqs for blast/mmseqs search, increase this when expected search match is not found. default : 100 |
| `--mafft-algorithm` | str | MAFFT algorithm for alignment (see MAFFT docs). default : auto |
| `--mafft-op` | float | MAFFT op (gap opening penalty) value, default : 1.3 |
| `--mafft-ep` | float | MAFFT ep value, default : 0.1 |
| `--trimal-algorithm` | str | trimAl algorithm for trimming (see trimAl docs). default : gt |
| `--trimal-gt` | float | trimAl -gt (gap threshold), between 0 and 1. default : 0.2 |
| `--allow-innertrimming` | flag | Turn off FunVIP adjustment to not to trim inner alignment columns, default: False |
| `--criterion` | str | Model-selection criterion for modeltest, one of [AIC, AICc, BIC]. default : BIC |
| `--noavx` | flag | Do not use AVX for RAxML, default: False |
| `--outgroupoffset` | int | outgroupoffset value. Highering this value may select more distant outgroup, default : 20 |
| `--nosuspicious` | flag | Do not include suspicious samples for sequence-set. Mostly for metabarcoding analysis. May deduce inaccurate result with problematic database. |
| `--terminate` | flag | Terminate FunVIP run when critical error detected |

## cache

| Option | Type | Description |
|---|---|---|
| `--nocachedb` | flag | Disable caching current search database. Use it if your database is too big for system directory, default : True |
| `--usecache` | flag | Use cached search database, turn off if your cached database makes error, default : True |

## save

| Option | Type | Description |
|---|---|---|
| `--tableformat` | str | Default format for table files, [csv, xlsx], default : csv |
| `--nosearchresult` | flag | Do not save blast/mmseqs search matrix, use when dataset gets too big and generates IO bottleneck, default: False |

## setting

| Option | Type | Description |
|---|---|---|
| `--preset` | str | [fast, accurate], or json formatted option config file. Check documentation for each preset, default : fast |

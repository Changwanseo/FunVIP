# FunVIP documentation

### usage: 
```
FunVIP [-h] [--query [QUERY ...]] [--db DB [DB ...]] [--gene [GENE ...]] [--email EMAIL | --api API] 
[--test TEST] [--thread THREAD] [--memory MEMORY] [--outdir OUTDIR] [--runname RUNNAME]
[--mode MODE] [--continue] [--step STEP] [--level LEVEL] [--all] [--confident CONFIDENT] 
[--search SEARCH] [--alignment ALIGNMENT] [--trim TRIM] [--modeltest MODELTEST] [--tree TREE] [--bscutoff BSCUTOFF]
[--highlight HIGHLIGHT] [--heightmultiplier HEIGHTMULTIPLIER] [--maxwordlength MAXWORDLENGTH] 
[--backgroundcolor [BACKGROUNDCOLOR ...]] [--outgroupcolor OUTGROUPCOLOR] [--ftype FTYPE]
[--fsize FSIZE] [--fsize_bootstrap FSIZE_BOOTSTRAP] [--verbose VERBOSE] 
[--maxoutgroup MAXOUTGROUP] [--collapsedistcutoff COLLAPSEDISTCUTOFF]
[--collapsebscutoff COLLAPSEBSCUTOFF] [--bootstrap BOOTSTRAP] [--solveflat] [--regex [REGEX ...]] 
[--cluster-cutoff CLUSTER_CUTOFF] [--cluster-evalue CLUSTER_EVALUE] 
[--cluster-wordsize CLUSTER_WORDSIZE] [--mafft-algorithm MAFFT_ALGORITHM] [--mafft-op MAFFT_OP] 
[--mafft-ep MAFFT_EP] [--trimal-algorithm TRIMAL_ALGORITHM] [--trimal-gt TRIMAL_GT] 
[--allow-innertrimming] [--criterion CRITERION] [--noavx] [--outgroupoffset OUTGROUPOFFSET] 
[--cachedb] [--usecache] [--tableformat TABLEFORMAT] [--nosearchresult] [--preset PRESET] [--version]
```

#### \* If more than one flags are available, either of them can be applied (i.e. both --email and -e are available to use)

## **Meta options**
**\-h, --help**: show this help message

**\--version**: Show version of FunVIP

## **Main options**

**\--query \[QUERY ...\], -q \[QUERY ...\]**:  
Query FASTA or table (csv, tsv, or xlsx) files. Data that you want to analyze. You can either write file names (if they are in current directory), or file path (either full path or relative path is okay). See figure (LINK) to how to generate query files.

**\--db DB \[DB ...\], -d DB \[DB ...\]** **&lt;<REQUIRED!!!!&gt;>**  
Database table (csv, tsv, or xlsx) files. You can either write file names (if they are in current directory), or file path (either full path or relative path is okay). See figure (LINK) to how to generate database files.

**\--gene \[GENE ...\], -g \[GENE ...\]** **&lt;<REQUIRED!!!!&gt;>**  
Gene names to be analyzed. Should be same to corresponding column name (Differences in capitalization are okay)

**\--preset PRESET \[fast\] | \[accurate\] | \[YAML_FORMAT_FILE\]**

**fast**: fast mode, use --auto for MAFFT and Fasttree for tree construction.  
**accurate**: accurate mode, use --localpair (l-ins-i) for MAFFT and RAxML for tree construction.  
**YAML_FORMAT_FILE**: your custom set of options saved in yaml format. See <https://github.com/Changwanseo/FunVIP/tree/main/funvip/preset> to find out how to make it. Default: fast.  

**\--email EMAIL, -e EMAIL**: E-mail notation to download data from GenBank (Without email, NCBI may ban your IP). Mandatory if your data includes GenBank accessions.

## **Test run**

**\--test \[Penicillium\]** | **\[Terrei\]**  
Use test datase.Penicillium dataset compromises data from GenMine paper. It will take 10+ minutes in 12core system.  
Terrei dataset compromises data from FunVIP paper. It will take 3+ minutes in 12c system. Terrei option requires additional email to also test downloading sequence from GenBank

## **Running options**

**\--thread \[THREAD\], -t \[THREAD\]**  
Threads to be used for pipeline. Default: system maximum. If FunVIP uses too much memory, decreasing thread may help reducing memory requirement

**\--memory \[MEMORY\], -m \[MEMORY\]**  
Max memory limit in gigbytes. Write in 'nG' form, ex: '16G'. Should be more than 4G. Default: system maximum. This limit only applies for MMseqs and iqtree. Therefore, if FunVIP exceeds memory in other steps, please decrease thread number instead.

**\--outdir \[OUTDIR\]**  
Out file location, Default: current directory. If you run FunVIP with **\--continue**, this option is mandatory to designate which previous run to use

**\--runname \[RUNNAME\]**  
Name prefix to current run. Default: current timestamp. If you run FunVIP with **\--continue**, this option is mandatory to designate which previous run to use

**\--mode \[MODE\]**  
Mode setup in species identification, **\[validation\] | \[identification\].  
**The mode decides species assignment of zero-length polytomy branch, if species name of the queries wree given. For example, two database sample, Penicillium A, and Penicillium B shows zero-length polytomy branch. In validation mode, if species name of the query is given as Penicillium A, FunVIP decides query as Penicillium A (Believes query information) In identification mode, if query is given as Penicillium A, FunVIP decides query as Penicillium sp. 1, because it is not fully confident that the query belongs to Penicillium A or Penicillium B with query sequence alone.

**\--continue**  
Continue from previous run. Should be used with **\--outdir**, **\--runname**, **\--step** to designate previous run directory and which steps to be continued

**\--step** **\[setup\] | \[search\] | \[cluster\] | \[align\] | \[trim\] | \[concatenate\] | \[modeltest\] | \[tree\] | \[visualize\] | \[report\]  
**Step to continue from previous run. Will be ignored without **\--continue** option.  
**\[setup\]**: Start from initial step, from input validation  
**\[search\]**: Start from BLAST/MMseqs search among sequences  
**\[cluster\]**: Start from group assignment and outgroup selection  
**\[align\]**: Start from MAFFT alignment  
**\[trim\]**: Start from alignment trimming  
**\[concatenate\]**: Start from concatenating alignments 
**\[modeltest\]**: Start from model selection  
**\[tree\]**: Start from phylogenetic tree construction  
**\[visualize\]**: Start from tree interpretation and visualization
**\[report\]**: Start from report

**\--level \[subseries\] | \[series\] | \[subsection\] | \[section\] | \[subtribe\] | \[subfamily\] | \[family\] | \[suborder\] | \[order\] | \[subclass\] | \[class\] | \[subphylum\] | \[phylum\] | \[subdivision\] | \[division\] | \[subkingdom\] | \[kingdom\]**  
Taxonomic level for each phylogenetic tree. Should be higher than species level

**\--all**  
Run FunVIP for all database sequences, regardless of corrresponding sequences exists in query. Default: False.  
If you don't have any query sequences, **\--all** flag will be automatically on.

**\--maxoutgroup MAXOUTGROUP**  
Maximum outgroup numbers to include in phylogenetic analysis, default: 3. We recommend at least 3 to maxoutgroup to confirm quartet hypothesis.

## **Method selection for each step of pipeline**

**\--search \[blast\] | \[mmseqs\]**  
Search methods to be used in selecting genes, groups and outgroups, Default: blast. BLAST is more sensitive, while mmseqs is faster in large number of sequences. MMseqs uses large space of memory, so more than 16GB or memory (RAM) is recommended

**\--trim \[trimal\] \[gblocks\] \[none\]**  
Trimming methods. Default: trimal. If **\[none\]** selected, trimming step will be skipped.  
Default trimming methods are modified to not to remove internal columns. If you want to use original behavior of trimal or gblocks, use with **\--allow-innertrimming**

**\--modeltest** **\[iqtree\] \[modeltestng\] \[none\]**  
Model test methods. Defaul : none. If **\[none\]** selected, modeltest step will be skipped  
Modeltest-ng is currently only available with Linux platform.

**\--tree \[fasttree\] \[iqtree\] \[raxml\]**  
Phylogenetic tree construction method. Default: fasttree

## **Visualization options for drawing phylogenetic tree**

**\--bscutoff \[BSCUTOFF\]**  
Bootstrap cutoff for visualization. Default: 70. If 70, bootstrap values which are equals or above 70 are represented in tree

**\--highlight \[HIGHLIGHT\]**  
Color to highlight query sequences in tree visualization. Default: #AA0000 (Crimson like color). Either in html svg recognizable string or hex code, i.e. “red”, “blue”, “magenta” or “#03030A”.

**\--heightmultiplier \[HEIGHTMULTIPLIER\]**
Height multiplier in drawing collapsed clades (triangles in the tree). Default: 6. Larger heightmultiplier expands phylogenetic tree, while smaller heightmultiplier reduces phylogenetic tree in height.

**\--maxwordlength \[MAXWORDLENGTH\]**  
Maximum letters to be shown in single line of tree annotation. Default: 48.

**\--backgroundcolor \[BACKGROUNDCOLOR ...\]**  
List of background colors to be shown in tree, Default: “#f4f4f4” “#c6c6c6” (Very light pink and pink). Input should be used with quotes, separated by space, and recommended to be used with hex codes. To remove background, use --backgroundcolor "#FFFFFF" "#FFFFFF" (Iterating white and white)

**\--outgroupcolor \[OUTGROUPCOLOR\]**  
Background color to indicate outgroup, default: #999999 (Dark grey). Use --outgroupcolor #FFFFFF to remove background color of outgroup.

**\--ftype \[FTYPE\]**  
Font for phylogenetic tree, default: Arial. The font should be installed in your computer, and some font may not work depending on your operating system.

**\--fsize \[FSIZE\]**  
Letter size for phylogenetic tree in pt. Default: 10

**\--fsize_bootstrap \[FSIZE_BOOTSTRAP\]**  
Letter size for bootstrap support in phylogenetic tree in pt. Default: 9

## **Advanced options**

**\--verbose \[0\] | \[1\] | \[2\] | \[3\], -v \[0\] | \[1\] | \[2\] | \[3\]**  
Verbosity level to inform user and write in log file. Default: 2  
**0**: error (Only inform when error occurs)  
**1**: warning (Inform warnings and errors)  
**2**: info (Default, Inform info, warnings, and errors)  
**3**: debug mode (Inform debug, info, warnings, and errors). Also changes multithreading behaviors to check problems in FunVIP. If you find bugs and seems to be deficiency of FunVIP, please run it again with **\--verbose 3** and send all result file to developer (use email or issue tabs in github repository).  

**\--collapsedistcutoff \[COLLAPSEDISTCUTOFF**\]  
Species delimitation criteria according to phylogenetic distance. Samples within tree distance of this value from least common ancestor will be considered as same species. Default: 0.01 (1% differences from least common ancestor are regarded as same species). This option does not deny topological evidence.

**\--collapsebscutoff \[COLLAPSEBSCUTOFF\]**  
Species delimitation criteria according to tree support percentage (bootstrap for IQTREE and RAxML, local support value for FastTree). Samples higher than this support of this value from least common ancestor will be considered as same species. Default: 100 (Do not use bootstrap as species delimitation criteria). This option does not deny topological evidence.

**\--bootstrap \[BOOTSTRAP\]**
Rapid boostrap number (RAxML) / Ultra-fast bootstrap number (IQTREE) for tree analysis. Will be ignored if fasttree is selected for tree method. Default: 1000. If you use IQTREE, bootstrap should be equal or higher than 1000 to prevent error.

**\--regex \[REGEX ...\]**  
Python regex groups to parse strain numbers from your input (see <https://docs.python.org/3/howto/regex.html>). Maybe useful if your sequence descriptions are dirty.  
For example, if you put --regex "SFC\[0-9\]{8}-\[0-9\]{2}" for option, and your query looks like  
“SFC20180902-01_Hymenochaetales_sp.\_ITS1.ab1”, only “SFC20180902-01” part will be parsed and used for final report.

**\--cluster-cutoff \[CLUSTER_CUTOFF\]**  
Minimum percent identity to be considered as same group in clustering analysis. Should be between 0 and 1.  
Default: 0.97 (Over than 97% identity to database sequence from search result will be regarded as the group). If your final tree includes too many unexpected queries, higher this value. If your final tree does not include expected query, lower this value (This rarely happens in genus level).

**\--cluster-evalue \[CLUSTER_EVALUE\]**  
E-value cutoffs for blast/mmseqs search. Over this value will be calculated and applied to analysis. Default: 0.0000001

**\--cluster-wordsize \[CLUSTER_WORDSIZE\]**  
Word size for blast/mmseqs search. Default: 7

**\--mafft-algorithm \[localpair\] | \[globalpair\] | \[auto\]**  
MAFFT algorithm for alignment, see mafft documents.

**\--mafft-op \[MAFFT_OP\]**  
MAFFT op (gap opening penalty) value. Default: 1.3. It is highly recommended to increase this value to 2 or 3 when your alignment does not seem to be correct. You may have to use 5 in extreme cases.

**\--mafft-ep \[MAFFT_EP**\]  
MAFFT ep (gap extension penalty) value. Default: 0.1. If there are big gaps in your alignments, lowering this value may help you.

**\--trimal-algorithm \[gt\] \[gappyout\] \[strict\] \[strictplus\]** 
TrimAl algorithm for trimming, see TrimAl documents for details (<https://trimal.readthedocs.io/en/latest/usage.html>).  
Will be ignored if Trimal is not selected for trimming option

**\--trimal-gt \[TRIMAL_GT\]**  
gt (Gap threshold) value for TrimAl (<https://trimal.readthedocs.io/en/latest/usage.html>). Default: 0.2. Will be ignored if TrimAl is not selected for trimming option

**\--allow-innertrimming**  
Turn off FunVIP adjustment to not to trim inner alignment columns. Default: False

**\--criterion \[AIC\] \[AICc\] \[BIC\]**  
Modeltest criterion to use, either AIC (Akaike’s Information Criterion), AICc (Akaike’s Information Criterion with a correction for a small design), or BIC (Bayesian information Criterion). Best matching model for designated criterion will be used for tree construction.

**\--noavx**  
Do not use AVX acceleration for RAxML. Default: False. Will be automatically turned on if your system does not have AVX instructions set.

**\--outgroupoffset \[OUTGROUPOFFSET\]**  
Outgroupoffset value in bitscore cutoff. Increasing this value will select more distant outgroup, and decreasing this value will select closer outgroup. If your outgroup is mixed with ingroup in the phylogenetic tree, try increasing this value. Default: 20

## **Saving options**

**\--cachedb**  
Cache current search database to run faster when you are working with the same database multiple times. Turn off if your database is too big for system directory, default: True

**\--usecache**  
Use cached search database, turn off if your previously cached database is making error, default: True

**\--matrixformat** MATRIXFORMAT  
Default format for search matrix files, \[tsv, csv, xlsx\], default: csv

**\--nosearchresult**  
Do not save blast/mmseqs search matrix, use when dataset gets too big and generates IO bottleneck

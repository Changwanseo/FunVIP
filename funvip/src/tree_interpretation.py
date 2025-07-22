# tree interpretation pipeline - collapse and visualize tree
from ete3 import (
    Tree,
    TreeStyle,
    NodeStyle,
    TextFace,
    CircleFace,
    RectFace,
    faces,
)
from Bio import SeqIO
from copy import deepcopy
from time import sleep
import lxml.etree as ET
import pandas as pd
from functools import lru_cache
from funvip.src.tool import sizeof_fmt
from funvip.src.tool import get_id, get_genus_species
from itertools import combinations
import tracemalloc
import dendropy
import collections
import psutil
import os
import re
import sys
import json

# Default zero length branch for concatenation
CONCAT_ZERO = 0  # for better binding

# For colored logging
bold_red = "\x1b[31;1m"
yellow = "\x1b[33;20m"
green = "\x1b[92m"
reset = "\x1b[0m"


## Get maximum tree distance among all leaf pairs in given tree
def get_max_distance(tree):
    max_distance = 0

    # To prevent affecting tree
    tree = deepcopy(tree)
    (farthest_node, max_distance) = tree.detach().get_farthest_node()
    return max_distance


## divide string by new line character to prevent long string from being cut
## by the maximum length of the string
def divide_by_max_len(string, max_len, sep=" "):
    final_string = ""
    tmp_string = ""

    for char in string:
        if char == sep:
            if len(tmp_string) < max_len:
                tmp_string += char
            else:
                tmp_string += char
                final_string += tmp_string + "\n"
                tmp_string = ""
        else:
            tmp_string += char
    final_string += tmp_string
    return final_string


# Per clade information
class Collapse_information:
    def __init__(self):
        self.query_list = []
        self.db_list = []
        self.outgroup = []
        self.clade = None  # partial tree clade
        self.leaf_list = []
        self.clade_cnt = 0  # if clade with same name exists, use this as counter
        self.collapse_type = (
            ""  # line - for single clade / triangle - for multiple clade
        )
        self.color = ""  # color after collapsed
        self.height = ""
        self.width = ""
        self.taxon = ""  # taxon name to be shown
        self.n_db = 0
        self.n_query = 0
        self.n_others = 0
        self.flat = False

    def __str__(self):
        return f"clade {self.taxon} with {len(self.leaf_list)} leaves"

    def __repr__(self):
        return f"clade {self.taxon} with {len(self.leaf_list)} leaves"


## concat two given clade object and return
def concat_clade(
    clade1,
    clade2,
    dist1=CONCAT_ZERO,
    dist2=CONCAT_ZERO,
    support1=1,
    support2=1,
    root_dist=CONCAT_ZERO,
    root_support=0,
):
    tmp = Tree()
    tmp.dist = root_dist
    tmp.support = root_support
    tmp.add_child(clade1, dist=dist1, support=support1)
    tmp.add_child(clade2, dist=dist2, support=support2)
    return tmp


## concat all given branches for concatenation
# clades were given in tuble, and root_dist is given
def concat_all(clade_tuple, root_dist, root_support=0):
    if len(clade_tuple) == 0:
        print("No clade input found, abort")
        raise Exception
    # If one clade were input, return self
    elif len(clade_tuple) == 1:
        return clade_tuple[0].copy("newick")
    # If two clades were input, concat it and return
    elif len(clade_tuple) == 2:
        return_clade = clade_tuple[0].copy("newick")
        return_clade = concat_clade(
            clade1=return_clade,
            clade2=clade_tuple[1].copy("newick"),
            dist1=return_clade.dist,
            dist2=clade_tuple[1].dist,
            support1=return_clade.support,
            support2=clade_tuple[1].support,
            root_dist=root_dist,
            root_support=root_support,
        )
    # If more than 3 clades were input, iteratively concat
    # If more than 2 species exists, and sp included, which taxon sp should be included cannot be decided
    # In that case, move sp clade to last
    elif len(clade_tuple) >= 3:
        return_clade = clade_tuple[0].copy("newick")
        for c in clade_tuple[1:-1]:
            return_clade = concat_clade(
                clade1=return_clade,
                clade2=c.copy("newick"),
                dist1=return_clade.dist,
                dist2=c.dist,
                support1=return_clade.support,
                support2=c.support,
            )
        return_clade = concat_clade(
            clade1=return_clade,
            clade2=clade_tuple[-1].copy("newick"),
            dist1=return_clade.dist,
            dist2=clade_tuple[-1].dist,
            support1=return_clade.support,
            support2=clade_tuple[-1].support,
            root_dist=root_dist,
            root_support=root_support,
        )
    else:
        raise Exception

    return return_clade


# Default tree style
class Tree_style:
    def __init__(self):
        self.ts = TreeStyle()
        self.ts.scale = 1000
        # self.ts.show_branch_length = True
        # self.ts.show_branch_support = True
        self.ts.branch_vertical_margin = 10
        self.ts.allow_face_overlap = True
        self.ts.children_faces_on_top = True
        self.ts.complete_branch_lines_when_necessary = False
        self.ts.extra_branch_line_color = "black"
        self.ts.margin_left = 200
        self.ts.margin_right = 200
        self.ts.margin_top = 200
        self.ts.margin_bottom = 200
        self.ts.show_leaf_name = False


@lru_cache(maxsize=10000)
def decide_type(query_list, db_list, outgroup, string, by="hash", priority="query"):
    query = False
    db = False

    query_hash_list = [FI.hash for FI in query_list]
    db_hash_list = [FI.hash for FI in db_list]
    outgroup_hash_list = [FI.hash for FI in outgroup]

    if by == "hash":
        if string in query_hash_list:
            return "query"
        elif string in db_hash_list:
            return "db"
        elif string in outgroup_hash_list:
            return "outgroup"
        else:
            return "none"

    else:
        print(
            f"{bold_red}[ERROR] DEVELOPMENTAL ERROR, UNEXPECTED by for decide_type{reset}"
        )
        raise Exception


# count number of taxons in the clade
def taxon_count(
    funinfo_dict, query_list, db_list, outgroup, clade, gene, count_query=False
):
    taxon_dict = {}

    for leaf in clade:
        FI = funinfo_dict[leaf.name]
        taxon = None
        if count_query == True:
            taxon = (FI.genus, FI.bygene_species[gene])
        elif (
            decide_type(
                query_list=query_list,
                db_list=db_list,
                outgroup=outgroup,
                string=leaf.name,
            )
            == "db"
            or decide_type(
                query_list=query_list,
                db_list=db_list,
                outgroup=outgroup,
                string=leaf.name,
            )
            == "outgroup"
        ):
            taxon = (FI.genus, FI.bygene_species[gene])

        if not (taxon is None):
            if not (taxon in taxon_dict):
                taxon_dict[taxon] = 1
            else:
                taxon_dict[taxon] += 1

    return taxon_dict


def genus_count(funinfo_dict, gene, clade):
    taxon_dict = {}

    for leaf in clade.iter_leaves():
        FI = funinfo_dict[leaf.name]
        if (
            decide_type(
                query_list=self.query_list,
                db_list=self.db_list,
                outgroup=self.outgroup,
                string=leaf.name,
            )
            == "db"
            or decide_type(
                query_list=self.query_list,
                db_list=self.db_list,
                outgroup=self.outgroup,
                string=leaf.name,
            )
            == "outgroup"
        ):
            if not ((FI.genus, FI.bygene_species[gene]) in taxon_dict):
                taxon_dict[(FI.genus, FI.bygene_species[gene])[0]] = 1
            else:
                taxon_dict[(FI.genus, FI.bygene_species[gene])[0]] += 1

    return taxon_dict


def designate_genus(funinfo_dict, query_list, db_list, outgroup, gene, clade):
    genus_dict = {}

    # Get genus_count
    for leaf in clade.iter_leaves():
        FI = funinfo_dict[leaf.name]
        if (
            decide_type(
                query_list=query_list,
                db_list=db_list,
                outgroup=outgroup,
                string=leaf.name,
            )
            == "db"
            or decide_type(
                query_list=query_list,
                db_list=db_list,
                outgroup=outgroup,
                string=leaf.name,
            )
            == "outgroup"
        ):
            if not ((FI.genus, FI.bygene_species[gene]) in genus_dict):
                genus_dict[(FI.genus, FI.bygene_species[gene])[0]] = 1
            else:
                genus_dict[(FI.genus, FI.bygene_species[gene])[0]] += 1

    if len(genus_dict) >= 2:  # if genus is not clear
        return "AMBIGUOUSGENUS"
    elif len(genus_dict) == 1:
        return list(genus_dict.keys())[0]
    else:
        return designate_genus(
            funinfo_dict=funinfo_dict,
            query_list=query_list,
            db_list=db_list,
            outgroup=outgroup,
            gene=gene,
            clade=clade.up,
        )


def decide_clade(
    funinfo_dict, query_list, db_list, outgroup, clade, gene, count_query=True
):
    taxon_dict = taxon_count(
        funinfo_dict=funinfo_dict,
        query_list=query_list,
        db_list=db_list,
        outgroup=outgroup,
        clade=clade,
        gene=gene,
        count_query=count_query,
    )
    if len(taxon_dict.keys()) == 0:
        return "query"
    else:
        return "db"


# this function finds major species of the clade
def find_majortaxon(
    funinfo_dict, query_list, db_list, outgroup, sp_cnt, clade, gene, opt=None
):
    taxon_dict = taxon_count(
        funinfo_dict=funinfo_dict,
        query_list=query_list,
        db_list=db_list,
        outgroup=outgroup,
        clade=clade,
        gene=gene,
    )
    max_value = 0
    major_taxon = ""

    for taxon in taxon_dict:
        if taxon_dict[taxon] > max_value:
            max_value = taxon_dict[taxon]
            major_taxon = taxon

    if major_taxon == "":
        if opt is None:
            # if major species not selected, try to match genus at least
            major_taxon = (
                designate_genus(
                    funinfo_dict,
                    query_list,
                    db_list,
                    outgroup,
                    gene,
                    clade,
                ),
                f"sp. {sp_cnt}",
            )

        elif opt.mode == "validation":  # in validation mode, try to follow query sp
            taxon_dict = taxon_count(
                funinfo_dict=funinfo_dict,
                query_list=query_list,
                db_list=db_list,
                outgroup=outgroup,
                clade=clade,
                gene=gene,
                count_query=True,
            )
            max_value = 0
            for taxon in taxon_dict:
                if taxon_dict[taxon] > max_value:
                    max_value = taxon_dict[taxon]
                    major_taxon = taxon

            if not (major_taxon[1].startswith("sp")):
                major_taxon = (
                    designate_genus(
                        funinfo_dict,
                        query_list,
                        db_list,
                        outgroup,
                        gene,
                        clade,
                    ),
                    f"sp. {sp_cnt}",
                )

        else:
            major_taxon = (
                designate_genus(
                    funinfo_dict,
                    query_list,
                    db_list,
                    outgroup,
                    gene,
                    clade,
                ),
                f"sp. {sp_cnt}",
            )

    return major_taxon


# decides if the clade is monophyletic
def is_monophyletic(
    funinfo_dict, query_list, db_list, outgroup, opt, sp_cnt, clade, gene, taxon
):
    taxon_dict = taxon_count(
        funinfo_dict=funinfo_dict,
        query_list=query_list,
        db_list=db_list,
        outgroup=outgroup,
        clade=clade,
        gene=gene,
    )
    # if taxon dict.keys() have 0 species: all query
    if len(taxon_dict.keys()) == 0:
        for children in clade.children:
            # if any of the branch length was too long for single clade
            if children.dist > opt.collapsedistcutoff:
                return False
            # or bootstrap is to distinctive
            elif children.support > opt.collapsebscutoff:
                return False
        return True
    # if taxon dict.keys() have 1 species: 1 kinds of species
    elif len(taxon_dict.keys()) == 1:
        # if taxon dict.keys() have only 1 species: group assigned
        for children in clade.children:
            other_children = list(set(clade.children) - set([children]))[0]

            # check query branch
            if find_majortaxon(
                funinfo_dict=funinfo_dict,
                query_list=query_list,
                db_list=db_list,
                outgroup=outgroup,
                sp_cnt=sp_cnt,
                clade=children,
                gene=gene,
            )[1].startswith("sp."):
                if children.dist > opt.collapsedistcutoff:
                    return False
                elif children.support > opt.collapsebscutoff:
                    return False
                elif other_children.dist > opt.collapsebscutoff:
                    return False
                elif other_children.dist > opt.collapsedistcutoff:
                    return False
        return True
    else:
        # more than 2 species : not monophyletic
        return False


def check_monophyletic(
    funinfo_dict, query_list, db_list, outgroup, opt, sp_cnt, clade, gene
):
    ## Start of local_check_monophyletic
    # if clade only has query species or not
    datatype = decide_clade(
        funinfo_dict=funinfo_dict,
        query_list=query_list,
        db_list=db_list,
        outgroup=outgroup,
        clade=clade,
        gene=gene,
        count_query=False,
    )

    # if only one clade, it is confirmly monophyletic
    if len(clade.children) == 1:
        return datatype, True

    # if additional_clustering option is on, check if basal group includes query seqs
    taxon = find_majortaxon(
        funinfo_dict=funinfo_dict,
        query_list=query_list,
        db_list=db_list,
        outgroup=outgroup,
        sp_cnt=sp_cnt,
        clade=clade,
        gene=gene,
    )

    if is_monophyletic(
        funinfo_dict,
        query_list,
        db_list,
        outgroup,
        opt,
        sp_cnt,
        clade,
        gene,
        taxon,
    ):
        return datatype, True
    else:
        return datatype, False


# Main tree information class
class Tree_information:
    def __init__(self, tree, Tree_style, group, gene, opt):
        self.tree_name = tree  # for debugging
        self.t = Tree(tree)
        self.t_publish = (
            None  # for publish tree - will substitute tree_original in long_term
        )
        self.dendro_t = dendropy.Tree.get(
            path=self.tree_name, schema="newick"
        )  # dendropy format for distance calculation

        # if support ranges from 0 to 1, change it from 0 to 100
        # b for branch
        support_set = set()
        for b in self.t.traverse():
            support_set.add(b.support)

        if max(support_set) <= 1:
            for b in self.t.traverse():
                b.support = int(100 * b.support)

        self.query_list = []
        self.db_list = []
        self.outgroup = []
        self.outgroup_leaf_name_list = []  # hash list of outgroup
        self.outgroup_group = []  # list of groups that outgroup sequences designated
        self.funinfo_dict = {}  # leaf.name (hash) : Funinfo

        self.sp_cnt = 1
        self.reserved_sp = set()

        self.Tree_style = Tree_style
        self.group = group
        self.gene = gene
        self.opt = opt

        self.collapse_dict = {}  # { taxon name : collapse_info }

        self.outgroup_clade = None
        self.bgstate = 1
        self.additional_clustering = True
        self.zero = 0.00000100000050002909

        self.flat_clades = []

    # to find out already existing new species number to avoid overlapping
    # e.g. avoid sp 5 if P. sp 5 already exsits in database
    def reserve_sp(self):  # does not seems to be working currently
        for leaf in self.t.iter_leaves():
            FI = self.funinfo_dict[leaf.name]
            taxon = (FI.genus, FI.ori_species)
            sys.stdout.flush()
            if taxon[1].split(" ")[0] in ("sp", "sp."):
                self.reserved_sp.add(" ".join(taxon[1].split(" ")[1:]))

    # Calculate zero length branch length cutoff with given tree and alignment
    def calculate_zero(self, alignment_file, gene, partition_dict):
        # tracemalloc.start()
        # Parse alignment
        seq_list = list(SeqIO.parse(alignment_file, "fasta"))

        # Check if tree leaves and alignments are consensus
        hash_list_tree = [leaf.name for leaf in self.t]
        hash_list_alignment = [seq.id for seq in seq_list]

        if collections.Counter(hash_list_tree) != collections.Counter(
            hash_list_alignment
        ):
            print(
                f"{bold_red}[ERROR] content of tree and alignment is not identical for {self.tree_name}{reset}"
            )
            raise Exception

        # Find identical or including pairs in alignment
        identical_pairs = []
        # different_pairs = []

        # make phylogenetic distance matrix
        pdc = self.dendro_t.phylogenetic_distance_matrix().as_data_table()._data

        diff_min = 999999
        for seq1, seq2 in combinations(seq_list, 2):
            if not (
                str(seq1.id).strip() == str(seq2.id).strip()
                or (seq1.id, seq2.id) in identical_pairs
                or (seq2.id, seq1.id) in identical_pairs
            ):
                """
                # Chenge unusable chars into gap
                seq1_str = str(seq1.seq).lower()
                seq2_str = str(seq2.seq).lower()

                for char in set(seq1_str) - {"a", "t", "g", "c", "-"}:
                    seq1_str = seq1_str.replace(char, "-")

                for char in set(seq2_str) - {"a", "t", "g", "c", "-"}:
                    seq2_str = seq2_str.replace(char, "-")
                """

                seq1_str, seq2_str = str(seq1.seq).lower(), str(seq2.seq).lower()
                seq1_str = "".join(
                    "-" if char not in "atgc-" else char for char in seq1_str
                )
                seq2_str = "".join(
                    "-" if char not in "atgc-" else char for char in seq2_str
                )

                identical_flag = True
                # To prevent distance among different region detected as zero in concatenated analysis
                overlapping_cnt = 0

                # For concatenated sequence alignment, identical sequnece should be checked by each partitions
                if gene == "concatenated":
                    len_dict = partition_dict["len"]
                    gene_order = partition_dict["order"]

                    valid_index = []

                    # calculate valid index to check
                    previous_index = 0

                    for gene in gene_order:
                        start = previous_index
                        end = previous_index + len_dict[gene] - 1

                        # Find the starting point
                        for n in range(previous_index, len_dict[gene] + previous_index):
                            if seq1_str[n] != "-" and seq2_str[n] != "-":
                                start = n
                                break

                        # Compare from the start
                        for n in range(
                            len_dict[gene] + previous_index - 1,
                            previous_index - 1,
                            -1,
                        ):
                            if seq1_str[n] != "-" and seq2_str[n] != "-":
                                end = n
                                break

                        valid_index.extend(range(start, end + 1))

                        previous_index += len_dict[gene]

                    # for valid part
                    for n in valid_index:
                        # connected with or to evaluate insertions or deletions
                        if seq1_str[n] != seq2_str[n]:
                            identical_flag = False
                        else:
                            overlapping_cnt += 1

                else:
                    start = 0
                    end = len(seq1_str) - 1
                    # calculate start and end
                    for n in range(len(seq1_str)):
                        if seq1_str[n] != "-" and seq2_str[n] != "-":
                            start = n
                            break

                    for n in range(len(seq1_str)):
                        if (
                            seq1_str[len(seq1_str) - n - 1] != "-"
                            and seq2_str[len(seq1_str) - n - 1] != "-"
                        ):
                            end = len(seq1_str) - n
                            break

                    # for valid part
                    for n in range(start, end):
                        # connected with or to evaluate insertions or deletions
                        if seq1_str[n] != seq2_str[n]:
                            identical_flag = False
                        else:
                            overlapping_cnt += 1

                # if identical pairs
                id1 = str(seq1.id).strip()
                id2 = str(seq2.id).strip()

                if identical_flag is True and overlapping_cnt > 0:
                    if pdc[id1][id2] > self.zero:
                        self.zero = pdc[id1][id2]

                # if different pairs
                elif identical_flag is False:
                    if pdc[id1][id2] < diff_min:
                        diff_min = pdc[id1][id2]

        if diff_min < self.zero:
            self.zero = diff_min - 0.00000001

        if self.opt.verbose >= 3:
            print(f"[DEBUG] End of calculate zero")
            process = psutil.Process(os.getpid())
            memory_info = process.memory_info()
            print(f"[DEBUG] RAM usage: {memory_info.rss / 1000 / 1000} MB")

        # I think also finding minimal distance between non-identical sequences are also needed
        return self.zero

    def reroot_outgroup(self, out):
        # Rerooting
        # Reroot should be done first because unrooted tree may have 3 children clades
        outgroup_leaves = []

        # Resolve polytomy before rerooting
        self.t.resolve_polytomy()

        # Check if outgroup sequences exists
        print(f"[INFO] Rerooting {self.outgroup} in {self.tree_name}")
        for leaf in self.t:
            if any(outgroup.hash in leaf.name for outgroup in self.outgroup):
                outgroup_leaves.append(leaf)

        self.outgroup_leaf_name_list = [leaf.name for leaf in outgroup_leaves]
        self.outgroup_group = list(
            set(self.funinfo_dict[leaf.name].adjusted_group for leaf in outgroup_leaves)
        )

        # find smallest monophyletic clade that contains all leaves in outgroup_leaves
        # reroot with outgroup_clade
        try:
            # For more than one outgroups, after rerooting, get_common_ancestor of outgroup again
            if len(outgroup_leaves) >= 2:
                self.outgroup_clade = self.t.get_common_ancestor(outgroup_leaves)
                self.t.set_outgroup(self.outgroup_clade)
                self.t.ladderize(direction=1)
                self.outgroup_clade = self.t.get_common_ancestor(outgroup_leaves)
            elif len(outgroup_leaves) == 1:
                self.outgroup_clade = outgroup_leaves[0]
                self.t.set_outgroup(self.outgroup_clade)
                self.t.ladderize(direction=1)
                self.outgroup_clade = outgroup_leaves[0]
            else:
                print(
                    f"{bold_red}[ERROR] no outgroup selected in {self.tree_name}{reset}"
                )
                raise Exception

            # If number of outgroup leaves and outgroup clade does not matches, paraphyletic
            if len(outgroup_leaves) != len(self.outgroup_clade):
                print(
                    f"{yellow}[WARNING] outgroup seems to be paraphyletic in {self.tree_name}{reset}"
                )

        except:
            print(f"{yellow}[WARNING] no outgroup selected in {self.tree_name}{reset}")

            outgroup_flag = False
            # if outgroup_clade is on the root side, reroot with other leaf temporarily and reroot again
            for leaf in self.t:
                if not (leaf in outgroup_leaves):
                    self.t.set_outgroup(leaf)
                    # Rerooting again while outgrouping gets possible
                    try:
                        self.outgroup_clade = self.t.get_common_ancestor(
                            outgroup_leaves
                        )
                        # print(f"Ancestor: {self.outgroup_clade}")
                        self.t.set_outgroup(self.outgroup_clade)
                        outgroup_flag = True
                        break
                    except:
                        pass

            if outgroup_flag is False:
                # never erase this for debugging
                print(
                    f"{bold_red}[ERROR] Outgroup not selected in {self.tree_name}{reset}"
                )
                print(
                    f"{bold_red}[ERROR] local variable outgroup_leaves : {outgroup_leaves}{reset}"
                )
                print(f"{bold_red}[ERROR] tree_info.outgroup : {self.outgroup}{reset}")
                print(
                    f"{bold_red}[ERROR] tree_info.outgroup_clade : {self.outgroup_clade}{reset}"
                )
                raise Exception

        self.Tree_style.ts.show_leaf_name = True

        # Add bootstrap
        # Before doing this, copy the tree to not interfere final tree
        copied_tree = deepcopy(self.t)

        for node in copied_tree.traverse():
            node.img_style["size"] = 0  # removing circles whien size is 0
            if len(node) > 1:  # Prevent bootstrap on single branch
                node.add_face(
                    TextFace(
                        f"{int(node.support)}",
                        fsize=self.opt.visualize.fsize_bootstrap,
                        fstyle="Arial",
                    ),
                    column=0,
                    position="float",
                )
        copied_tree.render(f"{out}", tree_style=self.Tree_style.ts)

        self.Tree_style.ts.show_leaf_name = False

        if self.opt.verbose >= 3:
            print(f"[DEBUG] End of reroot outgroup")
            process = psutil.Process(os.getpid())
            memory_info = process.memory_info()
            print(f"[DEBUG] RAM usage: {memory_info.rss / 1000 / 1000} MB")

    def collapse(self, collapse_info, clade, taxon):
        collapse_info.clade = clade
        collapse_info.taxon = taxon

        if len(clade) == 1:
            collapse_info.collapse_type = "line"
        elif len(clade) >= 2:
            collapse_info.collapse_type = "triangle"
        else:
            raise Exception

        if (
            any(
                decide_type(
                    query_list=self.query_list,
                    db_list=self.db_list,
                    outgroup=self.outgroup,
                    string=leaf.name,
                )
                == "query"
                for leaf in clade.iter_leaves()
            )
            == True
        ):
            collapse_info.color = self.opt.visualize.highlight
        else:
            collapse_info.color = "#000000"

        # all these things were ignored when type is line
        collapse_info.width = (
            get_max_distance(clade) * 1000
        )  # scale problem when visualizing
        collapse_info.height = len(clade) * self.opt.visualize.heightmultiplier

        # count query, db, others
        for leaf in clade.iter_leaves():
            if (
                decide_type(
                    query_list=self.query_list,
                    db_list=self.db_list,
                    outgroup=self.outgroup,
                    string=leaf.name,
                )
                == "db"
                or decide_type(
                    query_list=self.query_list,
                    db_list=self.db_list,
                    outgroup=self.outgroup,
                    string=leaf.name,
                )
                == "outgroup"
            ):
                collapse_info.leaf_list.append((leaf.name, "#000000", leaf.name))
                collapse_info.n_db += 1
            elif (
                decide_type(
                    query_list=self.query_list,
                    db_list=self.db_list,
                    outgroup=self.outgroup,
                    string=leaf.name,
                )
                == "query"
            ):
                collapse_info.leaf_list.append(
                    (leaf.name, self.opt.visualize.highlight, leaf.name)
                )
                collapse_info.n_query += 1
            else:
                print(
                    f"{bold_red}[ERROR] DEVELOPMENTAL ERROR : UNEXPECTED LEAF TYPE FOR {leaf.name}{reset}"
                )
                print(self.tree_name)
                print(f"Query: {sorted([FI.hash for FI in self.query_list])}")
                print(f"DB: {sorted([FI.hash for FI in self.db_list])}")
                print(f"Outgroup: {sorted([FI.hash for FI in self.outgroup])}")
                raise Exception

    # Species level delimitaion on tree
    def tree_search(self, clade, gene, opt=None):
        def local_generate_collapse_information(self, clade, opt=None):
            collapse_info = Collapse_information()
            collapse_info.query_list = self.query_list
            collapse_info.db_list = self.db_list
            collapse_info.outgroup = self.outgroup
            taxon = find_majortaxon(
                funinfo_dict=self.funinfo_dict,
                query_list=self.query_list,
                db_list=self.db_list,
                outgroup=self.outgroup,
                sp_cnt=self.sp_cnt,
                clade=clade,
                gene=gene,
                opt=opt,
            )
            self.collapse(collapse_info, clade, taxon)

            # counting new species
            if taxon[1].startswith("sp."):
                while 1:
                    self.sp_cnt += 1
                    if str(self.sp_cnt) in self.reserved_sp:
                        print(f"Skipping {self.sp_cnt} to avoid overlap in database")
                        continue
                    else:
                        break
            """
            print(
                f"[INFO] Generating collapse information on {self.group} {self.gene} for taxon {taxon}          ",
                end="\r",
            )
            """

            if not (taxon in self.collapse_dict):
                self.collapse_dict[taxon] = [collapse_info]
            else:
                self.collapse_dict[taxon].append(collapse_info)

        ## start of tree_search

        if len(clade.children) == 1:
            local_generate_collapse_information(clade, opt=opt)

            if self.opt.verbose >= 3:
                print(f"[DEBUG] End of Tree search with monophyletic branches")
                process = psutil.Process(os.getpid())
                memory_info = process.memory_info()
                print(f"[DEBUG] RAM usage: {memory_info.rss / 1000 / 1000} MB")

            return

        # In bifurcated clades
        elif len(clade.children) == 2:
            for child_clade in clade.children:
                # Calculate root distance between two childs to check flat
                self.flat = (
                    True if child_clade.dist <= self.opt.collapsedistcutoff else False
                )

                # Check if child clades are monophyletic
                datatype, monophyletic = check_monophyletic(
                    self.funinfo_dict,
                    self.query_list,
                    self.db_list,
                    self.outgroup,
                    self.opt,
                    self.sp_cnt,
                    child_clade,
                    gene,
                )

                # If monophyletic clade, generate collapse_info and finish
                if monophyletic is True:
                    local_generate_collapse_information(self, child_clade, opt=opt)
                # Else, do recursive tree search to divide clades
                else:
                    self.tree_search(child_clade, gene, opt=opt)

            if self.opt.verbose >= 3:
                print(f"[DEBUG] End of Tree search with bifurcated branches")
                process = psutil.Process(os.getpid())
                memory_info = process.memory_info()
                print(f"[DEBUG] RAM usage: {memory_info.rss / 1000 / 1000} MB")

            return

        # if error (more than two branches or no branches)
        else:
            print(
                f"{bold_red}[ERROR] DEVELOPMENTAL ERROR : FAILED TREE SEARCH ON LEAF {clade.children}{reset}"
            )
            raise Exception
        # end of tree_search

    # Reconstruct tree tree to solve flat branches
    def reconstruct(self, clade, gene, opt):
        sys.stdout.flush()  # for logging

        @lru_cache(maxsize=10000)
        def solve_flat(clade):
            # Check if the clade is consists of query db or both
            def consist(c):
                db, query = 0, 0
                for leaf in c:
                    if decide_type(
                        query_list=self.query_list,
                        db_list=self.db_list,
                        outgroup=self.outgroup,
                        string=leaf.name,
                    ) in ("db", "outgroup"):
                        db += 1
                    else:
                        query += 1

                if db == 0 and query == 0:
                    print(
                        f"{bold_red}[ERROR] DEVELOPMENTAL ON CONSIST, {c} {db} {query}{reset}"
                    )
                    raise Exception
                elif db == 0 and query != 0:
                    return "query"
                elif db != 0 and query == 0:
                    return "db"
                else:
                    return "both"

            # Get taxon of the given clade
            # c for clade (to remove redundancy to other variable: clade)
            def get_taxon(c, gene, mode="db"):
                # t for taxon : get taxon of the leaf
                def t(leaf):
                    try:
                        FI = self.funinfo_dict[leaf.name]
                        return (FI.genus, FI.bygene_species[gene])
                    except:
                        print(
                            f"{bold_red}[DEVELOPMENTAL ERROR] in leaf.name tree_interpretation.py line 869 {resety}"
                        )
                        raise Exception

                taxon_dict = {}

                # If only db in the clade
                if mode == "db":
                    for leaf in c:
                        if t(leaf) in taxon_dict:
                            taxon_dict[t(leaf)] += 1
                        else:
                            taxon_dict[t(leaf)] = 1

                    if len(taxon_dict) == 0:
                        print(
                            f"{bold_red}[DEVELOPMENTAL ERROR] in tree_interpretation.py line 884 {taxon_dict}\n {c}{reset}"
                        )
                        raise Exception
                    # If only one species in the clade
                    elif len(taxon_dict) == 1:
                        # If only one taxon here, return the taxon
                        return list(taxon_dict.keys())[0]
                    else:
                        # Else, return False
                        return False

                # If only query in the clade
                elif mode == "query":
                    for leaf in c:
                        if ("", "") in taxon_dict:
                            taxon_dict[("", "")] += 1
                        else:
                            taxon_dict[("", "")] = 1

                    if len(taxon_dict) == 0:
                        print(
                            f"{bold_red}[DEVELOPMENTAL ERROR] Error in tree_interpretation.py line 912 {taxon_dict}\n {c}{reset}"
                        )
                        raise Exception
                    elif len(taxon_dict) == 1:
                        return list(taxon_dict.keys())[0]
                    else:
                        max_taxon = ""
                        maximum = 0
                        for taxon in taxon_dict:
                            if taxon_dict[taxon] > maximum:
                                maximum, max_taxon = taxon_dict[taxon], taxon
                        return max_taxon

                # If db and query mixed in the clade
                elif mode == "both":
                    for leaf in c:
                        condition = False
                        # parse condition
                        if (
                            decide_type(
                                query_list=self.query_list,
                                db_list=self.db_list,
                                outgroup=self.outgroup,
                                string=leaf.name,
                                priority="query",
                            )
                            == "db"
                        ):
                            condition = True
                        elif (
                            decide_type(
                                query_list=self.query_list,
                                db_list=self.db_list,
                                outgroup=self.outgroup,
                                string=leaf.name,
                                priority="query",
                            )
                            == "outgroup"
                        ):
                            condition = True
                        elif (
                            opt.mode == "identification"
                            and decide_type(
                                query_list=self.query_list,
                                db_list=self.db_list,
                                outgroup=self.outgroup,
                                string=leaf.name,
                                priority="query",
                            )
                            == "query"
                        ):
                            condition = True

                        if condition is True:
                            if not (t(leaf) in taxon_dict):
                                taxon_dict[t(leaf)] = 1
                            else:
                                taxon_dict[t(leaf)] += 1

                    if len(taxon_dict) == 0:
                        print(f"{taxon_dict}\n {c}")
                        raise Exception
                    elif len(taxon_dict) == 1:
                        return list(taxon_dict.keys())[0]
                    else:
                        return False

            def seperate_clade(clade, gene, clade_list):
                for c in clade.children:
                    c_tmp = c.copy()
                    # zero clades
                    if c_tmp.dist <= self.zero:
                        # Original version was == instead of >= . Revert if error occurs
                        # What does the "len" means here? -> len means number of tips
                        # If only one tip
                        if len(c_tmp) <= 1:
                            # In the zero branch tip, the query with zero length should move to sp., because they cannot be fully determined
                            clade_list.append(
                                (
                                    get_taxon(c=c_tmp, gene=gene, mode=consist(c_tmp)),
                                    c_tmp,
                                    c_tmp.dist,
                                    c_tmp.support,
                                )
                            )
                        # If more than one tip
                        # I'm not sure if any of the recursion enters here, but just in case
                        else:
                            clade_list = seperate_clade(
                                clade=c_tmp, gene=gene, clade_list=clade_list
                            )

                    # non-zero clades
                    else:
                        c2 = self.reconstruct(c_tmp, gene, opt)
                        clade_list.append(
                            (
                                get_taxon(c2, gene, mode=consist(c2)),
                                c2,
                                c2.dist,
                                c2.support,
                            )
                        )

                return clade_list

            ## Start of function: solve_flat
            root_dist = clade.dist
            root_support = clade.support

            # Divide clade, all in the same level (flat branch)

            clade_list = seperate_clade(clade, gene, [])

            cnt = 0

            # count option.zero clades
            # result is from seperate_clade function
            # each of the result has list of (taxon, clade, dist, support)

            # count number of zero clades
            for result in clade_list:
                if result[2] <= self.zero:
                    cnt += 1

            # when entered to final leaf
            if len(clade_list) == 0:
                return clade

            else:
                # seperating db taxon and query taxon needed
                # clade_dict format: taxon : clade
                clade_dict = {}
                final_clade = []

                for result in clade_list:
                    # if clade does not have taxonomical information
                    # if final clade is mixed?
                    if result[0] is False:
                        final_clade.append(result[1])

                    else:
                        # if new taxa
                        if not (result[0] in clade_dict):
                            clade_dict[result[0]] = [result]
                        # if already checked taxa
                        else:
                            clade_dict[result[0]].append(result)

                candidate_zero_len_taxa = set(clade_dict.keys())

                if ("", "") in candidate_zero_len_taxa:
                    if len(candidate_zero_len_taxa - set({("", "")})) == 1:
                        # Move ("", "") (unknown species) to front to be combined to known species
                        taxon_to_merge = sorted(list(clade_dict.keys()), reverse=False)
                    else:
                        # If ("", "") (unknown species) matches to multiple species, move to back to not be combined with any of the species
                        taxon_to_merge = sorted(list(clade_dict.keys()), reverse=True)
                else:
                    # In other case the order is not important, sort them in any way
                    taxon_to_merge = sorted(list(clade_dict.keys()), reverse=False)

                # seperate this that order should not be affected by non-zero branch
                tmp_final_clade = []

                flat_issue_cnt = 0
                for taxon in taxon_to_merge:
                    l = clade_dict[taxon]  # l for list of results
                    r_list = [r[1] for r in l]  # result clade list
                    r_list.sort(key=lambda r: r.dist, reverse=True)
                    r_tuple = tuple(r_list)

                    # concatenate within taxon clades
                    concatenated_clade = concat_all(
                        clade_tuple=r_tuple, root_dist=0, root_support=0
                    )
                    tmp_final_clade.append(concatenated_clade)

                    if taxon != ("", "") and concatenated_clade.dist <= self.zero:
                        flat_issue_cnt += 1

                final_clade = tmp_final_clade + final_clade

                final = concat_all(
                    clade_tuple=tuple(final_clade),
                    root_dist=root_dist,
                    root_support=root_support,
                )

                # If 2+ species are connected to flat, add flat issue
                if flat_issue_cnt >= 2:
                    # print(taxon_to_merge)
                    for leaf in final:
                        self.flat_clades.append(leaf.name)

                # print(f"Status flat: {self.flat_clades}")

                return final.copy("newick")
            ## end of solve flat

        ## Start of function reconstruct
        if len(clade.children) in (0, 1):
            return clade.copy("newick")

        elif len(clade.children) == 2:
            clade1 = clade.children[0]
            clade2 = clade.children[1]

            # Solve flat
            if clade.dist <= self.zero:
                return solve_flat(clade).copy("newick")
            elif clade1.dist <= self.zero or clade2.dist <= self.zero:
                return solve_flat(clade).copy("newick")
            else:
                r_clade1 = self.reconstruct(clade1, gene, opt)
                r_clade2 = self.reconstruct(clade2, gene, opt)

            concatanated_clade = concat_clade(
                clade1=r_clade1,
                clade2=r_clade2,
                dist1=clade1.dist,
                dist2=clade2.dist,
                support1=clade1.support,
                support2=clade2.support,
                root_dist=clade.dist,
                root_support=clade.support,
            ).copy("newick")

            if self.opt.verbose >= 3:
                print(f"[DEBUG] End of reconstruct")
                process = psutil.Process(os.getpid())
                memory_info = process.memory_info()
                print(f"[DEBUG] RAM usage: {memory_info.rss / 1000 / 1000} MB")

            return concatanated_clade

        else:
            print(f"[ERROR] {clade} {clade.children} {len(clade.children)}")
            raise Exception
        ## end of reconstruct

    def get_bgcolor(self):
        return self.opt.visualize.backgroundcolor[
            self.bgstate % len(self.opt.visualize.backgroundcolor)
        ]

    ### Collapse tree
    def collapse_tree(self):
        # collect taxon string to decide in polishing
        # taxon_string_list = []
        taxon_string_dict = {}

        for collapse_taxon in self.collapse_dict:
            # collapse_info.taxon should be taxon
            for collapse_info in self.collapse_dict[collapse_taxon]:
                clade = collapse_info.clade

                # if only one clade with same name exists
                if len(self.collapse_dict[collapse_taxon]) == 1:
                    taxon_string = " ".join(collapse_info.taxon)
                    taxon_string_dict[taxon_string] = (
                        collapse_info.taxon[0],
                        collapse_info.taxon[1],
                        "",
                    )
                else:
                    collapse_info.clade_cnt = (
                        self.collapse_dict[collapse_taxon].index(collapse_info) + 1
                    )
                    taxon_string = (
                        f'{" ".join(collapse_info.taxon)}-{collapse_info.clade_cnt}'
                    )
                    taxon_string_dict[taxon_string] = (
                        collapse_info.taxon[0],
                        collapse_info.taxon[1],
                        f"-{collapse_info.clade_cnt}",
                    )

                # taxon_string_list.append(taxon_string)

                taxon_text = TextFace(
                    taxon_string,
                    fsize=int(self.opt.visualize.fsize),
                    ftype=self.opt.visualize.ftype,
                    fgcolor=collapse_info.color,
                )

                space_text = TextFace(
                    "  ",
                    fsize=int(self.opt.visualize.fsize),
                    ftype=self.opt.visualize.ftype,
                    fgcolor=collapse_info.color,
                )

                # hash list for further analysis
                string_hash_list = [x[0] for x in collapse_info.leaf_list]

                # order by translated
                string_hash_list.sort(key=lambda x: self.funinfo_dict[x].original_id)

                id_string = divide_by_max_len(
                    ",tmpseperator, ".join(string_hash_list),
                    self.opt.visualize.maxwordlength,
                )

                id_text = TextFace(
                    id_string,
                    fsize=int(self.opt.visualize.fsize),
                    ftype=self.opt.visualize.ftype,
                )

                if collapse_info.collapse_type == "triangle":
                    rectangle = RectFace(
                        width=collapse_info.width,
                        height=collapse_info.height,
                        fgcolor=collapse_info.color,
                        bgcolor=collapse_info.color,
                    )
                    clade.add_face(rectangle, 1, position="branch-right")

                clade.add_face(space_text, 2, position="branch-right")
                clade.add_face(taxon_text, 3, position="branch-right")
                clade.add_face(space_text, 4, position="branch-right")
                clade.add_face(id_text, 5, position="branch-right")

                # Get all tip names of the current working clade
                collapse_leaf_name_list = [x[0] for x in collapse_info.leaf_list]

                # Check if current working clade includes only outgroup sequences

                # Development, color unintended outgroup in outgroup color
                if any(
                    self.funinfo_dict[x].adjusted_group in self.outgroup_group
                    for x in collapse_leaf_name_list
                ):
                    # Change background color to outgroup color
                    clade.img_style["bgcolor"] = self.opt.visualize.outgroupcolor
                    # Do not draw collapsed clades
                    clade.img_style["draw_descendants"] = False

                # If not outgroup sequences
                else:
                    # If any of the sequence in clade considered to be in ingroup
                    if any(
                        self.funinfo_dict[x].adjusted_group == self.group
                        for x in collapse_leaf_name_list
                    ):
                        # color the background
                        clade.img_style["bgcolor"] = self.get_bgcolor()
                        # change background color for next clade to be discrimminated
                        # Currently disabled because it does not looks good
                        self.bgstate += 1
                    # Do not draw collapsed clades
                    clade.img_style["draw_descendants"] = False

        # show branch support above 70%
        for node in self.t.traverse():
            # change this part when debugging flat trees
            node.img_style["size"] = 0  # removing circles whien size is 0

            if node.support >= self.opt.visualize.bscutoff:
                # node.add_face without generating extra line
                # add_face_to_node
                node.add_face(
                    TextFace(
                        f"{int(node.support)}",
                        fsize=self.opt.visualize.fsize_bootstrap,
                        fstyle="Arial",
                    ),
                    column=0,
                    position="float",
                )

        return taxon_string_dict

    ### end of collapse tree

    ### edit svg image from initial output from ete3
    def polish_image(self, out, taxon_string_dict):
        # runname_group_gene.svg file enters here
        # the tree has rectangle collapsed group, tmpseperator, and hash

        # Render it to temporary svg file and re-parse with xml
        self.t.render(f"{out}", tree_style=self.Tree_style.ts)

        tree_xml = ET.parse(f"{out}")

        # in tree_xml, find all group
        _group = list(tree_xml.iter("{http://www.w3.org/2000/svg}g"))
        group_list = list(_group[0].findall("{http://www.w3.org/2000/svg}g"))

        # in tree_xml change all rectangles to polygon (trigangle)
        # This action may change when ete version 4 is available
        for group in group_list:
            if len(list(group.findall("{http://www.w3.org/2000/svg}rect"))) == 1:
                if group.get("fill") in (
                    "#000000",
                    self.opt.visualize.highlight,
                ):
                    rect = list(group.findall("{http://www.w3.org/2000/svg}rect"))[0]
                    rect.tag = "{http://www.w3.org/2000/svg}polygon"
                    rect.set(
                        "points",
                        f'{rect.get("width")},0 0,{int(rect.get("height"))/2} {rect.get("width")},{rect.get("height")}',
                    )

        # for taxons, gather all texts
        text_list = list(tree_xml.iter("{http://www.w3.org/2000/svg}text"))

        # Change this module to be worked with FI hash
        for text in text_list:
            # Decide if string of the tree is bootstrap, scale, taxon or id
            # taxon_list = [" ".join(x) for x in self.collapse_dict.keys()]
            try:
                int(text.text)
                text_type = "bootstrap"
            except:
                if text.text == "0.05":
                    text_type = "scale"
                elif any(
                    taxon.strip() == text.text.strip()
                    for taxon in taxon_string_dict.keys()
                ):
                    text_type = "taxon"
                else:
                    text_type = "hash"

            # relocate text position little bit for better visualization
            text.set("y", f'{int(float(text.get("y")))-2}')

            # Problem here: taxon name with cff., aff. etc includes underscore
            if text_type == "taxon":
                genus = taxon_string_dict[text.text][0]
                species = taxon_string_dict[text.text][1]
                rest = taxon_string_dict[text.text][2]

                # split genus, species, rest of parent into tspan
                text.text = ""
                tspan_list = []
                if genus != "":
                    tspan = ET.SubElement(text, "{http://www.w3.org/2000/svg}tspan")
                    tspan.text = genus + " "
                    tspan.set("font-style", "italic")

                if species != "":
                    try:
                        int(species)
                        tspan = ET.SubElement(text, "{http://www.w3.org/2000/svg}tspan")
                        tspan.text = species
                    except:
                        if "sp." in species:
                            tspan = ET.SubElement(
                                text, "{http://www.w3.org/2000/svg}tspan"
                            )
                            tspan.text = species
                        elif any(
                            term in species
                            for term in ["aff.", "f.", "var.", "cf.", "nom."]
                        ):
                            words = species.split(" ")
                            for n, word in enumerate(words):
                                tspan = ET.SubElement(
                                    text, "{http://www.w3.org/2000/svg}tspan"
                                )
                                tspan.text = word
                                if not (
                                    any(
                                        term in word
                                        for term in [
                                            "aff.",
                                            "f.",
                                            "var.",
                                            "cf.",
                                            "nom.",
                                        ]
                                    )
                                ):
                                    tspan.set("font-style", "italic")

                                if n < len(words) - 1:
                                    tspan.tail = " "

                        else:
                            tspan = ET.SubElement(
                                text, "{http://www.w3.org/2000/svg}tspan"
                            )
                            tspan.text = species
                            tspan.set("font-style", "italic")

                if rest != "":
                    tspan = ET.SubElement(text, "{http://www.w3.org/2000/svg}tspan")
                    tspan.text = rest + " "

            elif text_type == "bootstrap":
                int(text.text)
                # move text a little bit higher position
                text.set("y", f'{int(text.get("y"))-8}')
                text.set("x", f'{int(text.get("x"))+1}')

            elif text_type == "hash":
                words = text.text.split(",tmpseperator, ")
                text.text = ""
                for word in words:
                    tspan = ET.SubElement(text, "{http://www.w3.org/2000/svg}tspan")
                    try:
                        tspan.text = self.funinfo_dict[word.strip()].original_id + "  "
                        if self.funinfo_dict[word.strip()].color is not None:
                            try:
                                tspan.set("fill", self.funinfo_dict[word.strip()].color)
                            except:
                                print("DEVELOPMENTAL ERROR: Failed coloring tree")
                                raise Exception

                        elif (
                            decide_type(
                                query_list=self.query_list,
                                db_list=self.db_list,
                                outgroup=self.outgroup,
                                string=word,
                                by="hash",
                            )
                            == "query"
                        ):
                            tspan.set("fill", self.opt.visualize.highlight)
                    except:
                        pass
        # fit size of tree_xml to svg
        # find svg from tree_xml
        svg = list(tree_xml.iter("{http://www.w3.org/2000/svg}svg"))[0]

        # write to svg file
        tree_xml.write(
            out,
            encoding="utf-8",
            xml_declaration=True,
        )

        if self.opt.verbose >= 3:
            print(f"[DEBUG] End of Tree visualization")
            process = psutil.Process(os.getpid())
            memory_info = process.memory_info()
            print(f"[DEBUG] RAM usage: {memory_info.rss / 1000 / 1000} MB")

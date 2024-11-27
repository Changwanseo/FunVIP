# Dependencies and structures in tree_interpretation.py

def get_max_distance

def divide_by_max_length

class Collapse_information
	def __init__()
	def __str()__()
	def __repr()__()


def concat_clade

def concat_all
 - depends concat_clade

class Tree_style

class tree_information
	def __init__()
	def reserve_sp()
	def decide_type
	def calculate_zero
	def reroot_outgroup
	def taxon_count
	def genus_count
	def designate_genus
	def find_major_taxon
	def collapse
	def decide_clade
	def is_monophyletic
	def check_monophyletic
	def tree_search
		def local_check_monophyletic
			def decide_clade
			def is_monophyletic
		def local_generate_collapse_information
	def reconstruct
		def solve_flat
			def consist
			def get_taxon
			def seperate clade
	def get_bg_color
	def collapse_tree
	def polish_image

========================================
# How it works
main: run tree_interpretation_pipe.pipe_tree_interpretation
tree_interpretation_pipe: pipe_tree_interpretation
- make tree_interpretation_opt
- run pipe module_tree_interpretation

pipe_module_tree_interpretation: 
- move required variables of V to temporary variables
- run initialize_path (for get_genus_species)
- generate tree_information
- run calculate_zero
- run reroot_outgroup
- make hash
- run reserve_sp
- run reconstruct (solve flat)
- run ladderize
- run tree_search
	- run local_check_monophyletic
		- run decide_clade
		- run find majortaxon
		- run is_monophyletic
	- run local_generate_collapse_information if ended
	- run tree_search with subclades if not ended

## is_monophyletic code is same
	dependency: taxon_count, 


## the result differs because of "self"

# Global is_monophyletic code
def is_monophyletic(self, clade, gene, taxon):
    taxon_dict = self.taxon_count(clade, gene)
    if len(taxon_dict.keys()) == 0:
        for children in clade.children:
            if children.dist > self.opt.collapsedistcutoff:
                return False
            elif children.support > self.opt.collapsebscutoff:
                return False
        return True
    elif len(taxon_dict.keys()) == 1:
        for children in clade.children:
            other_children = list(set(clade.children) - set([children]))[0]
            if self.find_majortaxon(children, gene)[1].startswith("sp."):
                if children.dist > self.opt.collapsedistcutoff:
                    return False
                elif children.support > self.opt.collapsebscutoff:
                    return False
                elif other_children.dist > self.opt.collapsebscutoff:
                    return False
                elif other_children.dist > self.opt.collapsedistcutoff:
                    return False
        return True
    else:
        return False

# tree_search is_monophyletic code
def is_monophyletic(self, clade, gene, taxon):
    taxon_dict = self.taxon_count(clade, gene)
    if len(taxon_dict.keys()) == 0:
        for children in clade.children:
            if children.dist > self.opt.collapsedistcutoff:
                return False
            elif children.support > self.opt.collapsebscutoff:
                return False
        return True
    elif len(taxon_dict.keys()) == 1:
        for children in clade.children:
            other_children = list(set(clade.children) - set([children]))[0]
            if self.find_majortaxon(children, gene)[1].startswith("sp."):
                if children.dist > self.opt.collapsedistcutoff:
                    return False
                elif children.support > self.opt.collapsebscutoff:
                    return False
                elif other_children.dist > self.opt.collapsebscutoff:
                    return False
                elif other_children.dist > self.opt.collapsedistcutoff:
                    return False
        return True
    else:
        return False








genus_count(funinfo_dict, gene, clade)
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
		- run local_generate_collapse_information if ended
		- run tree_search with subclades if not ended




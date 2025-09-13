# For final visualization of the FunVIP tree
from ete3 import (
    Tree,
    TreeStyle,
    NodeStyle,
    TextFace,
    CircleFace,
    RectFace,
    faces,
)
from funvip.src.tool import get_genus_species
import pandas as pd
import lxml.etree as ET
import re

# Genus abbreviation
abbreviate = True
# True - "*", False - "T", "NT"
asterisk_type = False

# Tree style designation
ts = TreeStyle()
ts.scale = 5000
ts.branch_vertical_margin = 3
ts.allow_face_overlap = True
ts.children_faces_on_top = True
ts.complete_branch_lines_when_necessary = False
ts.extra_branch_line_color = "black"
ts.margin_left = 200
ts.margin_right = 200
ts.margin_top = 200
ts.margin_bottom = 200

name_offset = -1000


def visualize(
    outgroup,
    out,
    FunVIP_result,
    type_file,
    ml_tree,
    tmp_tree,
    bayesian_tree=None,
):

# Visualize with Lim's FEP style
def visualize(
    V=V,
    path=path,
    opt=opt
):

    # To make name in excel sheet meets to FunVIP format
    def excel_sheetname_legal(string):
        newick_illegal = ["'", "[", "]", ":", "*", "?", "/", "\\", ".", ",", "(", ")"]
        for i in newick_illegal:
            string = string.replace(i, "")

        string = string.replace("  ", " ")
        string = string.replace(" ", "_")
        # string = string.replace("-", "_")

        return string

    # get index from dataframe of leaf
    def get_index(leaf):
        print(leaf.name)

        dup_check = 0
        leaf_id_candidate = []
        for key in dict_indicated_name:
            if key in leaf.name:
                dup_check += 1
                leaf_id_candidate.append(key)

        if dup_check == 0:
            print(dup_check)
            # print(dict_indicated_name)
            print(leaf.name)
            raise Exception
        elif dup_check > 1:
            remove_set = set()
            for _id1 in leaf_id_candidate:
                for _id2 in leaf_id_candidate:
                    if _id1 != _id2 and _id1 in _id2:
                        remove_set.add(_id1)

            for r in remove_set:
                leaf_id_candidate.remove(r)

            if len(leaf_id_candidate) != 1:
                print(leaf_id_candidate)
                raise Exception

        leaf_id = leaf_id_candidate[0]
        for n, i in enumerate(df["ID"]):
            if df["ID"][n] == dict_indicated_name[leaf_id]:
                return n

        print(leaf.name)
        print(leaf_id)
        raise Exception

    ## Read result
    if FunVIP_result.endswith(".csv"):
        encodings = ["UTF-8", "latin-1", "cp1252"]
        for enc in encodings:
            try:
                df = pd.read_csv(FunVIP_result, encoding=enc, quoting=1)
                break
            except UnicodeDecodeError:
                continue
    elif FunVIP_result.endswith(".xlsx"):
        df = pd.read_excel(FunVIP_result)
    else:
        raise Exception

    # Make excel pair
    # {"legalized_ID" : "ID"}
    dict_indicated_name = {}
    for n, i in enumerate(df["ID"]):
        dict_indicated_name[excel_sheetname_legal(i)] = i


    def is_outgroup(leaf):
        if outgroup in leaf.name:
            print(leaf.name)
            return True
        else:
            return False

    # get id of from leaf
    def get_id(leaf):
        if is_outgroup(leaf):
            return leaf.name
        else:
            n = get_index(leaf)
            return df["ID"][n]

    # monophyletic check function
    def monophyletic(clade):
        set_taxon = set()
        for leaf in clade.iter_leaves():
            set_taxon.add(get_taxon(leaf))

        if len(set_taxon) == 1:
            return True
        elif len(set_taxon) == 0:
            raise Exception
        else:
            return False

    # check if clade includes query
    def query_found(clade):
        for leaf in clade.iter_leaves():
            if not (is_outgroup(leaf)):
                n = get_index(leaf)
                if df["DATATYPE"][n] == "query":
                    return True
        return False

    # check if clade is new_speces
    def new_species(clade):
        set_taxon = set()
        for leaf in clade.iter_leaves():
            set_taxon.add(get_taxon(leaf))

        if len(set_taxon) == 1:
            if "sp." in list(set_taxon)[0][1]:
                return True
            else:
                return False
        else:
            raise Exception

    # get taxon from leaf
    def get_taxon(leaf):
        if is_outgroup(leaf):
            return get_genus_species(leaf.name)
        else:
            n = get_index(leaf)
            print(leaf.name, df["SPECIES_ASSIGNED"][n])
            return get_genus_species(df["SPECIES_ASSIGNED"][n])

    # Update background color for monophyletic clade
    # input tree
    def find_and_color_monophyletic_clade(t):
        for clade in t.children:
            # Check if clade is monophyletic
            if is_monophyletic(funinfo_dict, query_list, db_list, outgroup, opt, sp_cnt, clade, gene, taxon):
                # If query is in it
                if query_found(clade):
                    if new_species(clade):
                        # Color with new species color
                        child.img_style["bgcolor"] = "paleturquoise"
                        # Change abbreviate to False if you want to show full name of new species

                        '''
                        for leaf in child.iter_leaves():
                            print(leaf.name)

                        print(f"DEBUG point 1: {list(child.iter_leaves())[0]}")
                        print(
                            f"DEBUG point 2: {get_taxon(list(child.iter_leaves())[0])}"
                        )
                        '''

                        taxon = get_taxon(list(child.iter_leaves())[0])
                        print(f"Taxon: {taxon}")
                        # raise Exception
                        child.add_face(
                            TextFace(
                                taxon[0] + " " + taxon[1],
                                ftype="Arial",
                                fsize=24,
                                fstyle="italic",
                            ),
                            column=0,
                            position="float-behind",
                        )

                    else:
                        # Color with recorded species color
                        child.img_style["bgcolor"] = "lightgrey"
                        taxon = get_genus_species(list(child.iter_leaves())[0].name)
                        child.add_face(
                            TextFace(
                                taxon[0] + " " + taxon[1],
                                ftype="Arial",
                                fsize=24,
                                fstyle="italic",
                            ),
                            column=0,
                            position="float-behind",
                        )

            else:
                find_and_color_monophyletic_clade(child)


    ## Start of the function
    # Read hashed tree
    t_ml = Tree(f"{path.out_tree}/hash/hash_{opt.runname}_{group}_{gene}.nwk")

    '''
    t_ml_leaves = set(str(l.name).replace("-", "_").replace("/", "") for l in t_ml)
    t_ml_leaves = set(
        name.replace("__", "") if name.endswith("__") else name for name in t_ml_leaves
    )
    t_ml_sub = t_ml.get_descendants()  
    '''

    find_and_color_monophyletic_clade(t_ml)

    # patch.py file
    #from patch import patch
    #patch()

    # t_ml.write(format=2, outfile="reallytemporarytreetocheck.nwk")
    # raise Exception

    # Visualization
    # Edit writings
    for leaf in t_ml.iter_leaves():
        leaf.img_style["draw_descendants"] = False
        space_text = TextFace(
            "  ",
            fsize=8,
            ftype="Arial",
            fgcolor="black",
        )

        # Else each of the text should be inspected
        if not ("sp." in get_taxon(leaf)[1]):
            try:
                n = get_index(leaf)
                ifquery = df["DATATYPE"][n]
            except:
                ifquery = "outgroup"
            # If query
            if ifquery == "query":
                id_text = id_text = TextFace(
                    get_id(leaf),
                    fsize=16,
                    ftype="Arial",
                    bold=True,
                    fgcolor="black",
                )

                leaf.add_face(space_text, 1, position="branch-right")
                leaf.add_face(id_text, 2, position="branch-right")

            else:
                genus_text = TextFace(
                    get_taxon(leaf)[0],
                    fsize=16,
                    ftype="Arial",
                    fstyle="italic",
                    fgcolor="black",
                )
                species_text = TextFace(
                    get_taxon(leaf)[1],
                    fsize=16,
                    ftype="Arial",
                    fstyle="italic",
                    fgcolor="black",
                )

                id_text = id_text = TextFace(
                    get_id(leaf),
                    fsize=16,
                    ftype="Arial",
                    fgcolor="black",
                )
                leaf.add_face(space_text, 1, position="branch-right")
                leaf.add_face(genus_text, 2, position="branch-right")
                leaf.add_face(space_text, 3, position="branch-right")
                leaf.add_face(species_text, 4, position="branch-right")
                leaf.add_face(space_text, 5, position="branch-right")
                leaf.add_face(space_text, 6, position="branch-right")
                leaf.add_face(id_text, 7, position="branch-right")

                # Add stars for type
                """
                if excel_sheetname_legal(get_id(leaf)) in dict_type:
                    if dict_type[excel_sheetname_legal(get_id(leaf))]:
                        if asterisk_type is True:
                            type_text = TextFace(
                                "*",
                                fsize=16,
                                ftype="Arial",
                                fgcolor="black",
                            )
                        else:
                            type_text = TextFace(
                                str(
                                    dict_type[excel_sheetname_legal(get_id(leaf))]
                                ).upper(),
                                fsize=12,
                                ftype="Arial",
                                fgcolor="black",
                            )

                        # leaf.add_face(space_text, 8, position="branch-right")
                        leaf.add_face(type_text, 9, position="branch-right")
                else:
                    print(f"{leaf.name} not found in dict_type")
                """

        # For new species, all of them are our species, so do not have to draw species
        else:
            id_text = id_text = TextFace(
                get_id(leaf),
                fsize=16,
                ftype="Arial",
                bold=True,
                fgcolor="black",
            )

            leaf.add_face(space_text, 1, position="branch-right")
            leaf.add_face(id_text, 2, position="branch-right")

        leaf.name = ""


    # show branch support above 70%
    for node in t_ml.traverse():
        # change this part when debugging flat trees
        node.img_style["size"] = 0  # removing circles whien size is 0
        node.img_style["vt_line_width"] = 2
        node.img_style["hz_line_width"] = 2

        if node.support >= 70 and node.support < 100.1:
            # node.add_face without generating extra line
            # add_face_to_node
            if bayesian_tree is None:
                node.add_face(
                    TextFace(
                        f"{int(node.support)}",
                        fsize=12,
                        fstyle="Arial",
                    ),
                    column=0,
                    position="float",
                )

            else:
                node.add_face(
                    TextFace(
                        f"{float(node.support)}",
                        fsize=12,
                        fstyle="Arial",
                    ),
                    column=0,
                    position="float",
                )

        # If bootstrap and bayesian pp both 100
        if node.support == 100.1:
            # node.img_style["vt_line_width"] = 4
            node.img_style["hz_line_width"] = 8

    # Imaging
    for node in t_ml.traverse():
        node.img_style["size"] = 0

    t_ml.render(out, tree_style=ts)

    # Polishing
    tree_xml = ET.parse(f"{out}")
    # in tree_xml, find all group
    _group = list(tree_xml.iter("{http://www.w3.org/2000/svg}g"))
    group_list = list(_group[0].findall("{http://www.w3.org/2000/svg}g"))

    # Find the width attribute in the root <svg> element
    # Hard-coding for Penicillium project, should be fixed later
    svg_width = 5000

    rect_y_coords = []

    # shorten height of background rectangle
    for group in group_list:
        if len(list(group.findall("{http://www.w3.org/2000/svg}rect"))) == 1:
            if group.get("fill") in (
                "#d3d3d3",
                "#afeeee",
            ):
                rect = list(group.findall("{http://www.w3.org/2000/svg}rect"))[0]
                rect.set("width", f'{int(float(rect.get("width"))+2000)}')
                rect.set("height", f'{int(float(rect.get("height"))-5)}')
                rect.set("y", f"{int(rect.get('y'))+2}")
                rect_y_coords.append(int(rect.get("y")))

    # for taxons, gather all texts
    text_list = list(tree_xml.iter("{http://www.w3.org/2000/svg}text"))

    # Change this module to be worked with FI hash
    for text in text_list:
        # Decide if string of the tree is bootstrap, scale, taxon or id
        # taxon_list = [" ".join(x) for x in self.collapse_dict.keys()]
        try:
            float(text.text)
            if float(text.text) == 0.05:
                text_type = "scale"
            else:
                text_type = "bootstrap"
        except:
            try:
                int(text.text)
                text_type = "bootstrap"
            except:
                if text.text == "0.05":
                    text_type = "scale"
                elif str(text.text) == "*" or str(text.text) in ("T", "NT", "ET", "HT"):
                    text_type = "type"
                elif text.get("font-size") == "24pt":
                    text_type = "species_group"
                    print(f"DEBUG species_group {text.text}")

                else:
                    text_type = ""

        # relocate text position little bit for better visualization
        text.set("y", f'{int(float(text.get("y")))-2}')

        if text_type == "bootstrap":
            if bayesian_tree is None:
                # float(text.text)
                # print(f"bootstrap: {text.text}")
                original_text = str(int(float(text.text)))
                # move text a little bit higher position
                text.set("y", f'{int(text.get("y"))+8}')
                text.set("x", f'{int(text.get("x"))-5}')
                # When pp value is under 1

                """
                text.text = text.text.replace(".0", "/0.")
                text.text = text.text.replace(".1", "/*")

                if text.text.endswith("."):
                    text.text = text.text + "00"
                    """

            else:
                float(text.text)
                # print(f"bootstrap: {text.text}")
                original_text = text.text
                # move text a little bit higher position
                text.set("y", f'{int(text.get("y"))+8}')
                text.set("x", f'{int(text.get("x"))-5}')
                # When pp value is under 1
                text.text = text.text.replace(".0", "/0.")

                try:
                    if float(text.text.split("/")[1]) < 0.95:
                        text.text = text.text.split("/")[0] + "/-"
                except:
                    # When pp value does not exists
                    if not (".") in text.text:
                        text.text = text.text + "/-"

                    pass

                text.text = text.text.replace(".1", "/*")

                if text.text.endswith("."):
                    text.text = text.text + "00"

                print(f"Original text {original_text} has been changed to {text.text}")

        if text_type == "species_group":
            parent = text.getparent()

            raw_transform = parent.get("transform")
            # print(raw_transform)
            transform_values = [
                float(x.strip())
                for x in raw_transform.split("(")[1].split(")")[0].split(",")
            ]
            # X value
            # print(transform_values[4])
            transform_values[4] = svg_width - len(text.text) * 16 - name_offset
            transform_values = [str(x) for x in transform_values]

            # print(transform_values)
            parent.set("transform", f'matrix({", ".join(transform_values)})')

            # Y value

            # parent.set("x", )
            # print(f'matrix({", ".join(transform_values)})')
            # text.set("x", f"{int(text.get('x')) + 800}")

        if text_type == "type":
            print("TYPE")
            print(text.text)
            print(text.get("y"), text.get("x"))
            text.set("y", f'{int(text.get("y"))-4}')
            text.set("x", f'{int(text.get("x"))+2}')

    print(svg_width)

    # fit size of tree_xml to svg
    # find svg from tree_xml
    svg = list(tree_xml.iter("{http://www.w3.org/2000/svg}svg"))[0]

    # write to svg file
    tree_xml.write(
        "polished_" + out,
        encoding="utf-8",
        xml_declaration=True,
    )


# This should be edited in FunVIP implemented version
# outgroup = "Coleosporium"
outgroup = "Tylospora"

# Output file name
# out = f"Delastria.svg"
out = f"Rhizopogon.svg"

# FunVIP result
FunVIP_result = f"./hypogeous_final.result.csv"

# Type file name - DB file that includes type (holotype etc) information
type_file = f"./hypogeous_final.result.csv"
ml_tree = f"./hypogeous_final_Rhizopogon_concatenated.nwk"
tmp_tree = "tmp.nwk"

FunVIP_visualizer(
    outgroup=outgroup,
    out=out,
    FunVIP_result=FunVIP_result,
    type_file=type_file,
    ml_tree=ml_tree,
    tmp_tree=tmp_tree,
)

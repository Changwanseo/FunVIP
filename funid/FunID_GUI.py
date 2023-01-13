import streamlit as st
import streamlit.components.v1 as components
from streamlit_tags import st_tags, st_tags_sidebar
from streamlit.script_run_context import get_script_run_ctx

import zipfile
import os
import json
import subprocess
import copy
import numpy as np
import shutil

from Bio import SeqIO
import pandas as pd

from src.initialize import isnan, gettype
from src.validation import validate_seq

INF = 999999

# setups for server hosting
ALLOW_DB_UPLOAD = False
ALLOW_DB_SELECT = False
ALLOW_OPTION_SELECT = False
ALLOW_GENE_SELECT = False
ALLOW_NO_QUERY = False
DESIGNATED_DB = ["MSP_Peni_BenA_reduced.xlsx"]
DESIGNATED_GENE = ["BenA"]
TITLE = "**PeniID**"
SUBTITLE = "A tree-based Penicillium identification pipeline for BenA sequences"
PAGETITLE = "PeniID - A tree based Penicillium identification pipeline"

# Runnning
identified = 0
runname = None
location = None
query_cnt = 0


def update_runname(_runname):
    runname = _runname


def update_location(_location):
    location = _location


# Option generator
def set_option(
    expander,
    obj,
    type_,
    criterion,
    value=np.NaN,
    min_=np.NaN,
    max_=np.NaN,
    default=np.NaN,
    solve=False,
    help_="",
):
    if isnan(obj):
        obj = ""

    if isnan(help_):
        help_ = ""

    if type_ is int:
        obj = int(obj)
        default = int(default)

    if type_ is bool:
        return expander.selectbox(
            label=criterion,
            index=[True, False].index(obj),
            options=[True, False],
            help=help_,
        )

    if not (isnan(value)):

        value = value.replace(" ", "").split(",")
        return expander.selectbox(
            label=criterion, index=value.index(obj), options=value, help=help_
        )

    elif not (isnan(min_)) and not (isnan(max_)):
        if type_ is int:
            return expander.number_input(
                label=criterion,
                value=int(obj),
                min_value=int(min_),
                max_value=int(max_),
                help=help_,
            )
        elif type_ is float:
            return expander.number_input(
                label=criterion,
                value=float(obj),
                min_value=float(min_),
                max_value=float(max_),
                help=help_,
            )
        else:
            raise Exception

    else:
        if type_ is list:
            obj = obj.replace(" ", "").split(",")
            return st_tags_sidebar(label=criterion, value=obj, text=help_)
            # return expander.text_input(label=criterion, value=obj, help=help_)

        else:
            return expander.text_input(label=criterion, value=obj, help=help_)


# for designated value
def set_value(type_, value=np.NaN):

    if isnan(value) is True:
        value = ""

    if type_ is int:
        return int(value)
    elif type_ is float:
        return float(value)
    elif type_ is bool:
        return bool(value)
    elif type_ is str:
        return str(value)
    elif type_ is list:
        return value.replace(" ", "").split(",")
    else:
        print(type_)
        raise Exception


st.set_page_config(page_title=PAGETITLE, layout="wide")

# session_state = Sessionstate.get(identified=False)
st.title(TITLE)
st.subheader(SUBTITLE)


# Parse options
input_option = {}
expander_list = ["BASIC", "METHOD", "VISUALIZE", "OPTIONAL", "ADVANCED"]


# Option sidebars
if ALLOW_OPTION_SELECT is True:
    df = pd.read_excel("./Option_manager.xlsx")  # read setup files
    df.set_index("option", inplace=True)
    if not ("option") in st.session_state:
        st.session_state["option"] = False

    with st.sidebar.form(key="Input_option"):

        st.header("Option")

        option_submit = st.form_submit_button(label="Confirm")

        with st.expander("Upload Option.config file"):
            option = st.file_uploader("Upload Option.config file", key="Option")
        st.subheader("OR")

        with open("./Default_options.config", "r") as fp:
            example_option = json.load(fp)
            input_option = copy.deepcopy(example_option)

        dict_expander = {}

        # make setup dataframe to json
        pre_dict = {}

        for ex in expander_list:
            if ex == "BASIC":
                dict_expander[ex] = st.expander(ex, expanded=True)
            else:
                dict_expander[ex] = st.expander(ex, expanded=False)

            input_option[ex] = {}

        for option_ in df.index:
            if bool(df["host"][option_]) is True:
                input_option[df["class"][option_]][option_] = set_option(
                    dict_expander[df["class"][option_]],
                    obj=df["host_default"][option_],
                    type_=gettype(df["type"][option_]),
                    criterion=option_,
                    value=df["select"][option_],
                    min_=df["host_min"][option_],
                    max_=df["host_max"][option_],
                    default=df["default"][option_],
                    solve=df["solve"][option_],
                    help_=df["help"][option_],
                )
            else:
                input_option[df["class"][option_]][option_] = set_value(
                    type_=gettype(df["type"][option_]),
                    value=df["host_default"][option_],
                )
else:
    st.session_state["option"] = True
    df = pd.read_excel("./Option_manager.xlsx")  # read setup files
    df.set_index("option", inplace=True)
    option_submit = True

    for ex in expander_list:
        input_option[ex] = {}

    # use default option when select is not available
    for option_ in df.index:
        input_option[df["class"][option_]][option_] = set_value(
            type_=gettype(df["type"][option_]),
            value=df["host_default"][option_],
        )


if ALLOW_DB_SELECT is True:
    col_db, col_query = st.columns(2)
else:
    col_query = st.container()


# DB input
if ALLOW_DB_SELECT is True:
    with col_db.form("DB Input"):
        st.subheader("Database")
        DB_select = st.multiselect(
            label="Choose",
            options=[
                file
                for file in os.listdir("./DB")
                if any(
                    file.endswith(x)
                    for x in [".xlsx", ".csv", ".ftr", ".feather", ".parquet"]
                )
                and not (file.startswith("~"))
            ],
        )

        if len(DB_select) > 0:
            dict_expander = {}
            for i in DB_select:
                dict_expander[i] = st.expander(f"View {i}")
                dict_expander[i].dataframe(pd.read_excel(f"./DB/{i}"))

        # get uploaded db
        if ALLOW_DB_UPLOAD is True:
            DB_upload = st.file_uploader(
                "Upload", key="DB", accept_multiple_files=True, type=[".xlsx", ".csv"]
            )
            # DB_upload_confirm = col1.button(label="Confirm", key="DB_upload_confirm")

            if DB_upload:
                with st.spinner("Uploading DB"):
                    # somethings
                    pass

        db_submit = st.form_submit_button(label="DB confirm")
        if not ("DB") in st.session_state:
            st.session_state["DB"] = False
else:
    DB_select = DESIGNATED_DB
    db_submit = True
    st.session_state["DB"] = True

# Query input
with col_query.form("Query Input"):
    st.subheader("Query")
    # get query
    query_upload = st.file_uploader(
        "Upload",
        key="Query upload",
        accept_multiple_files=True,
        type=[".fasta", ".fa", ".fsa", ".fna", ".xlsx", ".csv"],
    )
    query_write = st.text_area(label="Sequence in fasta format")
    query_no = False
    if ALLOW_NO_QUERY is True:
        st.write("OR")
        query_no = st.checkbox("Run without Query")
    query_submit = st.form_submit_button(label="Query confirm")
    if not ("Query") in st.session_state:
        st.session_state["Query"] = False

# Gene select
if ALLOW_GENE_SELECT is True:
    with col_db.form("Gene select"):
        gene_set = set()
        for db in DB_select:
            df = pd.read_excel(f"./DB/{db}")
            for c in list(df):
                gene_set.add(c)
        not_seq = [
            "Accession",
            "Genus",
            "Species",
            "Section",
            "Type",
            "Substrate",
            "Country",
            "Source",
        ]
        for x in not_seq:
            gene_set.discard(x)

        gene_select = st.multiselect(label="Choose", options=list(gene_set))
        gene_submit = st.form_submit_button(label="Gene selection")
        if not ("Gene") in st.session_state:
            st.session_state["Gene"] = False
else:
    gene_select = DESIGNATED_GENE
    gene_submit = True
    st.session_state["Gene"] = True


with st.form("Identify"):
    st.header("Run")
    run = st.form_submit_button(label="Identify")

# When option submit
if option_submit:
    st.session_state["option"] = input_option
    print(input_option)

# When DB submit
if db_submit:
    st.session_state["DB"] = DB_select


# When Query submit
if query_submit:

    # remove previous query files in the location
    # after using session id, we might not need this

    # for file in os.listdir("./Query"):
    #    if not (file.startswith(".")):
    #        os.remove(f"./Query/{file}")

    # For run without query option
    if query_no is True:
        st.session_state["Query"] = True

    else:
        query_flag = 0

        # save input files first
        if len(query_upload) >= 1:
            for file in query_upload:
                with open(
                    f"./Query/{str(get_script_run_ctx().session_id)}_{query_cnt}.fasta",
                    "wb",
                ) as fw:
                    fw.write(file.getbuffer())
                query_cnt += 1

            query_flag = 1

        if len(query_write) >= 1:
            # Check if writen form is available fasta
            write_result = validate_seq(query_write)
            # than save if to file
            if not (write_result[0]) == "invalid":
                if write_result[0] == "seqrecord":
                    print(str(get_script_run_ctx().session_id))
                    print(write_result[1])
                    SeqIO.write(
                        write_result[1],
                        f"./Query/{str(get_script_run_ctx().session_id)}_{query_cnt}.fasta",
                        "fasta",
                    )
                    query_cnt += 1

                else:
                    print("[ERROR] Developmental error in sequence parsing")
                    raise Exception

                query_flag = 1
            else:
                print(write_result)
                st.warning("Invalid sequence")

        if query_flag == 0:
            st.warning("No query uploaded!")
        else:
            st.session_state["Query"] = query_flag
            st.success("Query confirmed!")


if gene_submit:
    st.session_state["Gene"] = gene_select

if run:
    if (
        st.session_state["Query"] is not (False)
        and st.session_state["DB"] is not (False)
        and st.session_state["option"] is not (False)
        and st.session_state["Gene"] is not (False)
    ):
        input_option["BASIC"]["DB"] = st.session_state["DB"]
        input_option["BASIC"]["GENE"] = st.session_state["Gene"]
        with open(f"{get_script_run_ctx().session_id}.config", "w") as f:
            json.dump(input_option, f, indent=4)

        st.session_state["run"] = "run"
        # df = pd.DataFrame("./Result/Penicillium_test_2/Section Assignment.xlsx")
    else:
        if st.session_state["Query"] is False:
            st.warning("Query should be confirmed!")
        if st.session_state["DB"] is False:
            st.warning("DB should be confirmed!")
        if st.session_state["option"] is False:
            st.warning("Option should be confirmed!")
        if st.session_state["Gene"] is False:
            st.warning("Gene should be confirmed!")


# DB selection
if "run" in st.session_state:
    if st.session_state.run == "run":
        # sessionstate.identified = True
        print("Updating run status")

        # Running
        with st.spinner("Running FunID"):

            query_string = " ".join(
                [
                    f"./Query/{get_script_run_ctx().session_id}_{i}.fasta"
                    for i in range(st.session_state["Query"])
                ]
            )

            # push db option
            os.system(
                f"python ./Run.py -q {query_string} -c {get_script_run_ctx().session_id}.config"
            )

        st.session_state["run"] = "finished_run"

if "run" in st.session_state:
    if st.session_state.run == "finished_run":

        # Figure out Runname
        runname = input_option["BASIC"]["RUN_NAME"]
        path_root = f"./Result/{runname}"

        # Check if run with same name exists
        if input_option["METHOD"]["NEW_RUN"] is True:
            if os.path.exists(path_root):
                # if same name exists, try to add numbers at the end to discriminate
                i = 1
                while 1:
                    if os.path.exists(path_root + "_%s" % str(i)):
                        i += 1
                    else:
                        if not (i == 1):
                            path_root = path_root + "_%s" % str(i - 1)
                            runname = runname + "_%s" % str(i - 1)
                            break
                            # Generate working directory

        location = f"./Result/{runname}"  # should be changed by run

        st.header("Identification Result")
        shutil.make_archive(f"./Result/{runname}", "zip", f"./Result/{runname}")
        with open(f"./Result/{runname}.zip", "rb") as f:
            Download = st.download_button(
                "Download all data", f, file_name=f"{runname}.zip"
            )
        Identification_result_table = st.expander("Final Identification result")
        BLAST_result_table = st.expander("BLAST result")

        Tree_viewer = st.expander("View tree")

        gene_list = st.session_state["Gene"]
        if len(st.session_state["Gene"]) > 1 and not (
            "concatenated" in st.session_state["Gene"]
        ):
            gene_list.append("concatenated")

        section_set = set()
        for gene in gene_list:
            svg_list = [
                file
                for file in os.listdir(f"{location}/Tree/")
                if file.startswith(f"{runname}_") and file.endswith(f"_{gene}.svg")
            ]
            for svg in svg_list:
                section_set.add(
                    svg.replace(f"{runname}_", "").replace(f"_{gene}.svg", "")
                )

        tree_gene = Tree_viewer.selectbox(
            label="Gene",
            options=gene_list,
        )

        tree_section = Tree_viewer.selectbox(
            label="Section",
            options=list(section_set),
        )

        Identification_result_table.dataframe(
            pd.read_excel(f"{location}/{runname}_Identification_result.xlsx")
        )

        with st.spinner("Loading tree"):
            Tree_viewer.image(
                f"{location}/Tree/{runname}_{tree_section}_{tree_gene}.svg",
                width=None,
                use_column_width="always",
            )

    # DB to file
    # Query to file
    # Option to file

    # Run command

    # Result
    # Result table
    # Tree data
    # Report
    # Log
    # Zip total result

    # st.download_button("Download", "data")

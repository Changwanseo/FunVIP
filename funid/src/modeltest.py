from funid.src import ext
import os, shutil, copy
import logging

cmd_translator_modeltestng_fasttree = {
    "JC": "",
    "GTR": "-gtr",
    "JC+G4": "-gamma",
    "GTR+G4": "-gtr -gamma",
}
cmd_translator_modelfinder_fasttree = {
    "JC": "",
    "GTR": "-gtr",
    "JC+G4": "-gamma",
    "GTR+G4": "-gtr -gamma",
}
cmd_translator_modeltestng_raxml = {
    "JC": "--JC",
    "JC+G4": "--JC",
    "JC+I": "--JC",
    "JC+G4+I": "--JC",
    "HKY": "--HKY85",
    "HKY+G4": "--HKY85",
    "HKY+I": "--HKY85",
    "HKY+G4+I": "--HKY85",
    "K80": "--K80",
    "K80+G4": "--K80",
    "TrNef": "--GTRCAT",
    "TrNef+G4": "--GTRGAMMA",
    "TrNef+I": "--GTRCATI",
    "TrNef+G4+I": "--GTRGAMMAI",
    "TPM1": "--GTRCAT",
    "TPM1+G4": "--GTRGAMMA",
    "TPM1+I": "--GTRCATI",
    "TPM1+G4+I": "--GTRGAMMAI",
    "TPM2": "--GTRCAT",
    "TPM2+G4": "--GTRGAMMA",
    "TPM2+I": "--GTRCATI",
    "TPM2+G4+I": "--GTRGAMMAI",
    "TPM3": "--GTRCAT",
    "TPM3+G4": "--GTRGAMMA",
    "TPM3+I": "--GTRCATI",
    "TPM3+G4+I": "--GTRGAMMAI",
    "TIM1": "--GTRCAT",
    "TIM1+G4": "--GTRGAMMA",
    "TIM1+I": "--GTRCATI",
    "TIM1+G4+I": "--GTRGAMMAI",
    "TIM2": "--GTRCAT",
    "TIM2+G4": "--GTRGAMMA",
    "TIM2+I": "--GTRCATI",
    "TIM2+G4+I": "--GTRGAMMAI",
    "TIM3": "--GTRCAT",
    "TIM3+G4": "--GTRGAMMA",
    "TIM3+I": "--GTRCATI",
    "TIM3+G4+I": "--GTRGAMMAI",
    "TVMef": "--GTRCAT",
    "TVMef+G4": "--GTRGAMMA",
    "TVMef+I": "--GTRCATI",
    "TVMef+G4+I": "--GTRGAMMAI",
    "SYM": "--GTRCAT",
    "SYM+G4": "--GTRGAMMA",
    "SYM+I": "--GTRCATI",
    "SYM+G4+I": "--GTRGAMMAI",
    "F81": "--GTRCATX",
    "F81+G4": "--GTRGAMMAX",
    "F81+I": "--GTRCATIX",
    "F81+G4+I": "--GTRGAMMAIX",
    "TrN": "--GTRCATX",
    "TrN+G4": "--GTRGAMMAX",
    "TrN+I": "--GTRCATIX",
    "TrN+G4+I": "--GTRGAMMAIX",
    "TPM1uf": "--GTRCATX",
    "TPM1uf+G4": "--GTRGAMMAX",
    "TPM1uf+I": "--GTRCATIX",
    "TPM1uf+G4+I": "--GTRGAMMAIX",
    "TPM2uf": "--GTRCATX",
    "TPM2uf+G4": "--GTRGAMMAX",
    "TPM2uf+I": "--GTRCATIX",
    "TPM2uf+G4+I": "--GTRGAMMAIX",
    "TPM3uf": "--GTRCATX",
    "TPM3uf+G4": "--GTRGAMMAX",
    "TPM3uf+I": "--GTRCATIX",
    "TPM3uf+G4+I": "--GTRGAMMAIX",
    "TIM1uf": "--GTRCATX",
    "TIM1uf+G4": "--GTRGAMMAX",
    "TIM1uf+I": "--GTRCATIX",
    "TIM1uf+G4+I": "--GTRGAMMAIX",
    "TIM2uf": "--GTRCATX",
    "TIM2uf+G4": "--GTRGAMMAX",
    "TIM2uf+I": "--GTRCATIX",
    "TIM2uf+G4+I": "--GTRGAMMAIX",
    "TIM3uf": "--GTRCATX",
    "TIM3uf+G4": "--GTRGAMMAX",
    "TIM3uf+I": "--GTRCATIX",
    "TIM3uf+G4+I": "--GTRGAMMAIX",
    "TVM": "--GTRCATX",
    "TVM+G4": "--GTRGAMMAX",
    "TVM+I": "--GTRCATIX",
    "TVM+G4+I": "--GTRGAMMAIX",
    "GTR": "--GTRCATX",
    "GTR+G4": "--GTRGAMMAX",
    "GTR+I": "--GTRCATIX",
    "GTR+G4+I": "--GTRGAMMAIX",
}


def parse_model(modeltest_file, opt, path):

    model_list = []
    model_cmd = "no model found"

    if opt.method.modeltest == "iqtree":

        # Check iqtree - iqtree pair has misentered here
        if opt.method.tree == "iqtree":
            logging.error(
                f"DEVELOPMENTAL ERROR. MODELTEST {opt.method.modeltest} entered while {opt.method.tree} selected"
            )
            raise Exception

        # read modeltest result file
        try:
            with open(modeltest_file, "r", encoding="UTF-8") as f:
                lines = f.readlines()
            logging.info(f"Successfully parsed {modeltest_file}")
        except:
            logging.error(f"Cannot parse modeltest file for {modeltest_file}")
            raise Exception

        # Parse modeltest file
        flag = 0
        for line in lines:
            # finish parsing from new line break
            if line == "\n" and flag == 1:
                break
            # parse when modeltest result existing region
            if flag == 1:
                model_list.append(line.split(" ")[0])

            # indicing modeltest result start
            if (
                "Model                  LogL         AIC      w-AIC        AICc     w-AICc         BIC      w-BIC"
                in line
            ):
                # parsing flag
                flag = 1

        # Designate modeltest method by tree method
        if opt.method.tree == "fasttree":
            for model in model_list:
                # if model is available model for tree method take it
                if model in cmd_translator_modelfinder_fasttree:
                    # convert to commands
                    model_cmd = cmd_translator_modelfinder_fasttree[model]
                    break

            if model_cmd == "no model found":
                logging.error(f"DEVELOPMENTAL ERROR. FAILED PARSING APPROPRIATE MODEL")
                logging.error(f"{model_list}")
                raise Exception

        elif opt.method.tree == "raxml":
            for model in model_list:
                # if model is available model for tree method take it
                if model in cmd_translator_modelfinder_raxml:
                    # convert to coommands
                    model_cmd = cmd_translator_modelfinder_raxml[model]
                    break

            if model_cmd == "no model found":
                logging.error(f"DEVELOPMENTAL ERROR. FAILED PARSING APPROPRIATE MODEL")
                logging.error(f"{model_list}")
                raise Exception

        else:
            logging.error(
                f"DEVELOPMENTAL ERROR. INAPPROPRIATE TREE METHOD {opt.method.tree} selected"
            )
            raise Exception

    elif opt.method.modeltest == "modeltest-ng":

        ### Testing block ###
        # Parse files

        try:
            with open(modeltest_file, "r", encoding="UTF-8") as f:
                lines = f.readlines()
                status = None

            # Fasttree needs converter
            if opt.method.tree == "fasttree":
                flag = 0  # 1 : Criterion found, 2: --- line found 1st time(while parsing), 0 : parsing ended
                cnt = 1  # line count
                for line in lines:
                    print(f"flag {flag} / [line]: {line}")

                    if (
                        "model              K            lnL          score          delta    weight"
                        in line
                    ):
                        # Check if the right chunk selected
                        if "AICc" in line and opt.criterion == "AICc":
                            flag = 1
                        elif "AIC" in line and opt.criterion == "AIC":
                            flag = 1
                        elif "BIC" in line and opt.criterion == "BIC":
                            flag = 1

                    elif flag == 1 and "--------" in line:
                        # Skipping one line
                        flag = 2
                    elif flag == 2 and not ("------" in line):
                        logging.info(f"[line]")
                        logging.info(line)
                        # Parsing objects
                        model_list.append(line.split(f"{cnt}  ")[1].split(" ")[0])
                        if cnt < 10:
                            cnt += 1
                        else:
                            cnt = 1
                            # Finish all 10 parsed
                            break
                    elif flag == 2 and "------" in line:
                        # reset block
                        flag = 0

                # get commands with available model
                for model in model_list:
                    if model in cmd_translator_modeltestng_fasttree:
                        model_cmd = cmd_translator_modeltestng_fasttree[model]
                        break

            # For IQTREE and RAxML, just easily parse commands
            else:
                for line in lines:
                    # Check current criterion
                    if "AICc" in line:
                        status = "AICc"
                    elif "AIC" in line and not "AICc" in line:
                        status = "AIC"
                    elif "BIC" in line:
                        status = "BIC"

                    if status == opt.criterion:
                        if line.strip().startswith(">"):
                            if opt.method.tree == "fasttree":
                                logging.error(
                                    f"DEVELOPMENTAL ERROR. MODELTEST-NG COMMAND PARSING STEP entered while {opt.method.tree} selected"
                                )
                            elif opt.method.tree == "iqtree":
                                if "iqtree" in line:
                                    model_cmd = f'-m {line.split("-m ")[1]}'
                                    break
                            elif opt.method.tree == "raxml":
                                if "raxmlHPC-SSE3" in line:
                                    model_cmd = (
                                        f'-m {line.split("-m ")[1].split(" -n ")[0]}'
                                    )
                                    break

                # Currently using BIC model
                model_cmd = model_dict["BIC"]

        except:
            logging.warning(f"Cannot parse modeltest file for {modeltest_file}")

    else:
        logging.error(
            f"DEVELOPMENTAL ERROR. INAPPROPRIATE MODELTEST METHOD {opt.method.modeltest} selected"
        )
        raise Exception

    return model_cmd


def cleanup(modeltest_file):
    pass


# main function
def modeltest(V, path, opt) -> dict:

    group_dict = V.dict_dataset
    model_dict = copy.deepcopy(group_dict)

    for group in group_dict:
        for gene in group_dict[group]:
            if opt.method.modeltest.lower() == "modeltest-ng":
                # As modeltest-ng shows only top 10 results,
                # we should reduce number of models to prevent none of the model fits
                if opt.method.tree == "fasttree":
                    models = "-m JC,GTR -h ug"
                else:
                    models = ""

                # run modeltest-ng
                ext.Modeltest_ng(
                    fasta=f"{path.out_alignment}/{opt.runname}_MAFFT_{group}_{gene}.fasta",
                    models=models,
                    out=f"{path.out_modeltest}/{opt.runname}_{group}_{gene}",
                    thread=opt.thread,
                )

                # Change result name
                try:
                    shutil.move(
                        f"{path.out_modeltest}/{opt.runname}_{group}_{gene}.out",
                        f"{path.out_modeltest}/{opt.runname}_{group}_{gene}.modeltest",
                    )
                except:
                    logging.warning(
                        f"Cannot parse modeltest file for {path.out_modeltest}/{opt.runname}_{group}_{gene}.out. Running with default option"
                    )

                # parse model from modeltest-ng result
                model_dict[group][gene] = parse_model(
                    modeltest_file=f"{path.out_modeltest}/{opt.runname}_{group}_{gene}.modeltest",
                    path=path,
                    opt=opt,
                )
            elif opt.method.modeltest.lower() == "iqtree":
                if not (opt.method.tree == "iqtree"):

                    # Run IQTREE ModelFinder
                    ext.ModelFinder(
                        fasta=f"{path.out_alignment}/{opt.runname}_MAFFT_{group}_{gene}.fasta",
                        opt=opt,
                        path=path,
                        thread=opt.thread,
                    )

                    # Move modeltest result to appropriate location
                    try:
                        shutil.move(
                            f"{path.out_alignment}/{opt.runname}_MAFFT_{group}_{gene}.fasta.iqtree",
                            f"{path.out_modeltest}/{opt.runname}_{group}_{gene}.modeltest",
                        )
                    except:
                        logging.warning(
                            f"Cannot parse modeltest file for {path.out_alignment}/{opt.runname}_MAFFT_{group}_{gene}.fasta.iqtree. Running with default option"
                        )

                    # parse model from IQTREE ModelFinder result
                    model_dict[group][gene] = parse_model(
                        modeltest_file=f"{path.out_modeltest}/{opt.runname}_{group}_{gene}.modeltest",
                        path=path,
                        opt=opt,
                    )
                else:
                    logging.info(
                        "IQTREE will perform ModelFinder internally in tree construction step, skipping in modeltest step"
                    )
                    model_dict[group][gene] = "skip"

            else:  # including opt.model_method.lower() == "none":
                if opt.method.tree == "raxml":
                    logging.info("Skipping modeltest. Using GTRGAMMA as default")
                    model_dict[group][gene] = "-m GTRGAMMA"
                else:
                    logging.info("Skipping modeltest.")
                    model_dict[group][gene] = "skip"

    logging.info(model_dict)

    return model_dict

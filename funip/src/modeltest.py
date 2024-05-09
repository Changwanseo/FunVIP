from funip.src import ext
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
cmd_translator_modelfinder_raxml = {
    "JC": "-m JC",
    "JC+G4": "-m JC",
    "JC+I": "-m JC",
    "JC+G4+I": "-m JC",
    "HKY": "-m HKY85",
    "HKY+G4": "-m HKY85",
    "HKY+I": "-m HKY85",
    "HKY+G4+I": "-m HKY85",
    "K80": "-m K80",
    "K80+G4": "-m K80",
    "TrNef": "-m GTRCAT",
    "TrNef+G4": "-m GTRGAMMA",
    "TrNef+I": "-m GTRCATI",
    "TrNef+G4+I": "-m GTRGAMMAI",
    "TPM1": "-m GTRCAT",
    "TPM1+G4": "-m GTRGAMMA",
    "TPM1+I": "-m GTRCATI",
    "TPM1+G4+I": "-m GTRGAMMAI",
    "TPM2": "-m GTRCAT",
    "TPM2+G4": "-m GTRGAMMA",
    "TPM2+I": "-m GTRCATI",
    "TPM2+G4+I": "-m GTRGAMMAI",
    "TPM3": "-m GTRCAT",
    "TPM3+G4": "-m GTRGAMMA",
    "TPM3+I": "-m GTRCATI",
    "TPM3+G4+I": "-m GTRGAMMAI",
    "TIM1": "-m GTRCAT",
    "TIM1+G4": "-m GTRGAMMA",
    "TIM1+I": "-m GTRCATI",
    "TIM1+G4+I": "-m GTRGAMMAI",
    "TIM2": "-m GTRCAT",
    "TIM2+G4": "-m GTRGAMMA",
    "TIM2+I": "-m GTRCATI",
    "TIM2+G4+I": "-m GTRGAMMAI",
    "TIM3": "-m GTRCAT",
    "TIM3+G4": "-m GTRGAMMA",
    "TIM3+I": "-m GTRCATI",
    "TIM3+G4+I": "-m GTRGAMMAI",
    "TVMef": "-m GTRCAT",
    "TVMef+G4": "-m GTRGAMMA",
    "TVMef+I": "-m GTRCATI",
    "TVMef+G4+I": "-m GTRGAMMAI",
    "SYM": "-m GTRCAT",
    "SYM+G4": "-m GTRGAMMA",
    "SYM+I": "-m GTRCATI",
    "SYM+G4+I": "-m GTRGAMMAI",
    "F81+FO": "-m GTRCATX",
    "F81+FO+G4": "-m GTRGAMMAX",
    "F81+FO+I": "-m GTRCATIX",
    "F81+FO+G4+I": "-m GTRGAMMAIX",
    "TrN+FO": "-m GTRCATX",
    "TrN+FO+G4": "-m GTRGAMMAX",
    "TrN+FO+I": "-m GTRCATIX",
    "TrN+FO+G4+I": "-m GTRGAMMAIX",
    "TPM1uf+FO": "-m GTRCATX",
    "TPM1uf+FO+G4": "-m GTRGAMMAX",
    "TPM1uf+FO+I": "-m GTRCATIX",
    "TPM1uf+FO+G4+I": "-m GTRGAMMAIX",
    "TPM2uf+FO": "-m GTRCATX",
    "TPM2uf+FO+G4": "-m GTRGAMMAX",
    "TPM2uf+FO+I": "-m GTRCATIX",
    "TPM2uf+FO+G4+I": "-m GTRGAMMAIX",
    "TPM3uf+FO": "-m GTRCATX",
    "TPM3uf+FO+G4": "-m GTRGAMMAX",
    "TPM3uf+FO+I": "-m GTRCATIX",
    "TPM3uf+FO+G4+I": "-m GTRGAMMAIX",
    "TIM1uf+FO": "-m GTRCATX",
    "TIM1uf+FO+G4": "-m GTRGAMMAX",
    "TIM1uf+FO+I": "-m GTRCATIX",
    "TIM1uf+FO+G4+I": "-m GTRGAMMAIX",
    "TIM2uf+FO": "-m GTRCATX",
    "TIM2uf+FO+G4": "-m GTRGAMMAX",
    "TIM2uf+FO+I": "-m GTRCATIX",
    "TIM2uf+FO+G4+I": "-m GTRGAMMAIX",
    "TIM3uf+FO": "-m GTRCATX",
    "TIM3uf+FO+G4": "-m GTRGAMMAX",
    "TIM3uf+FO+I": "-m GTRCATIX",
    "TIM3uf+FO+G4+I": "-m GTRGAMMAIX",
    "TVM+FO": "-m GTRCATX",
    "TVM+FO+G4": "-m GTRGAMMAX",
    "TVM+FO+I": "-m GTRCATIX",
    "TVM+FO+G4+I": "-m GTRGAMMAIX",
    "GTR+FO": "-m GTRCATX",
    "GTR+FO+G4": "-m GTRGAMMAX",
    "GTR+FO+I": "-m GTRCATIX",
    "GTR+FO+G4+I": "-m GTRGAMMAIX",
    "F81": "-m GTRCATX",
    "F81+G4": "-m GTRGAMMAX",
    "F81+I": "-m GTRCATIX",
    "F81+G4+I": "-m GTRGAMMAIX",
    "TrN": "-m GTRCAT",
    "TrN+G4": "-m GTRGAMMA",
    "TrN+I": "-m GTRCATI",
    "TrN+G4+I": "-m GTRGAMMAI",
    "TPM1uf": "-m GTRCAT",
    "TPM1uf+G4": "-m GTRGAMMA",
    "TPM1uf+I": "-m GTRCATI",
    "TPM1uf+G4+I": "-m GTRGAMMAI",
    "TPM2uf": "-m GTRCAT",
    "TPM2uf+G4": "-m GTRGAMMA",
    "TPM2uf+I": "-m GTRCATI",
    "TPM2uf+G4+I": "-m GTRGAMMAI",
    "TPM3uf": "-m GTRCAT",
    "TPM3uf+G4": "-m GTRGAMMA",
    "TPM3uf+I": "-m GTRCATI",
    "TPM3uf+G4+I": "-m GTRGAMMAI",
    "TIM1uf": "-m GTRCAT",
    "TIM1uf+G4": "-m GTRGAMMA",
    "TIM1uf+I": "-m GTRCATI",
    "TIM1uf+G4+I": "-m GTRGAMMAI",
    "TIM2uf": "-m GTRCAT",
    "TIM2uf+G4": "-m GTRGAMMA",
    "TIM2uf+I": "-m GTRCATI",
    "TIM2uf+G4+I": "-m GTRGAMMAI",
    "TIM3uf": "-m GTRCAT",
    "TIM3uf+G4": "-m GTRGAMMA",
    "TIM3uf+I": "-m GTRCATI",
    "TIM3uf+G4+I": "-m GTRGAMMAI",
    "TVM": "-m GTRCAT",
    "TVM+G4": "-m GTRGAMMA",
    "TVM+I": "-m GTRCATI",
    "TVM+G4+I": "-m GTRGAMMAI",
    "GTR": "-m GTRCAT",
    "GTR+G4": "-m GTRGAMMA",
    "GTR+I": "-m GTRCATI",
    "GTR+G4+I": "-m GTRGAMMAI",
    "F81+F": "-m GTRCAT",
    "F81+F+G4": "-m GTRGAMMA",
    "F81+F+I": "-m GTRCATI",
    "F81+F+G4+I": "-m GTRGAMMAI",
    "TrN+F": "-m GTRCAT",
    "TrN+F+G4": "-m GTRGAMMA",
    "TrN+F+I": "-m GTRCATI",
    "TrN+F+G4+I": "-m GTRGAMMAI",
    "TPM1uf+F": "-m GTRCAT",
    "TPM1uf+F+G4": "-m GTRGAMMA",
    "TPM1uf+F+I": "-m GTRCATI",
    "TPM1uf+F+G4+I": "-m GTRGAMMAI",
    "TPM2uf+F": "-m GTRCAT",
    "TPM2uf+F+G4": "-m GTRGAMMA",
    "TPM2uf+F+I": "-m GTRCATI",
    "TPM2uf+F+G4+I": "-m GTRGAMMAI",
    "TPM3uf+F": "-m GTRCAT",
    "TPM3uf+F+G4": "-m GTRGAMMA",
    "TPM3uf+F+I": "-m GTRCATI",
    "TPM3uf+F+G4+I": "-m GTRGAMMAI",
    "TIM1uf+F": "-m GTRCAT",
    "TIM1uf+F+G4": "-m GTRGAMMA",
    "TIM1uf+F+I": "-m GTRCATI",
    "TIM1uf+F+G4+I": "-m GTRGAMMAI",
    "TIM2uf+F": "-m GTRCAT",
    "TIM2uf+F+G4": "-m GTRGAMMA",
    "TIM2uf+F+I": "-m GTRCATI",
    "TIM2uf+F+G4+I": "-m GTRGAMMAI",
    "TIM3uf+F": "-m GTRCAT",
    "TIM3uf+F+G4": "-m GTRGAMMA",
    "TIM3uf+F+I": "-m GTRCATI",
    "TIM3uf+F+G4+I": "-m GTRGAMMAI",
    "TVM+F": "-m GTRCAT",
    "TVM+F+G4": "-m GTRGAMMA",
    "TVM+F+I": "-m GTRCATI",
    "TVM+F+G4+I": "-m GTRGAMMAI",
    "GTR+F": "-m GTRCAT",
    "GTR+F+G4": "-m GTRGAMMA",
    "GTR+F+I": "-m GTRCATI",
    "GTR+F+G4+I": "-m GTRGAMMAI",
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
        print(modeltest_file)
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
                    # convert to commands
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

    # By group and by gene
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
                    fasta=f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta",
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
                        fasta=f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta",
                        opt=opt,
                        path=path,
                        thread=opt.thread,
                    )

                    # Move modeltest result to appropriate location
                    try:
                        shutil.move(
                            f"{path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta.iqtree",
                            f"{path.out_modeltest}/{opt.runname}_{group}_{gene}.modeltest",
                        )
                    except:
                        logging.warning(
                            f"Cannot parse modeltest file for {path.out_alignment}/{opt.runname}_trimmed_{group}_{gene}.fasta.iqtree. Running with default option"
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

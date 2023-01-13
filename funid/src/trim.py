from funid.src import ext
from funid.src.opt_generator import opt_generator
import multiprocessing as mp
import shutil


def pipe_trimming(V, path, opt):

    trimming_opt = opt_generator(V, opt, path, step="trimming")

    # run multiprocessing start
    if opt.verbose < 3:
        p = mp.Pool(opt.thread)
        if opt.method.trim.lower() == "gblocks":
            trimming_result = p.starmap(ext.Gblocks, trimming_opt)
        elif opt.method.trim.lower() == "trimal":
            trimming_result = p.starmap(ext.Trimal, trimming_opt)
        else:
            trimming_result = p.starmap(shutil.copy, trimming_opt)
        p.close()
        p.join()

    else:
        # non-multithreading mode for debugging
        trimming_result = []
        for option in trimming_opt:
            if opt.method.trim.lower() == "gblocks":
                trimming_result.append(ext.Gblocks(*option))
            elif opt.method.trim.lower() == "trimal":
                trimming_result.append(ext.Trimal(*option))
            else:
                trimming_result.append(shutil.copy(*option))
    return V, path, opt

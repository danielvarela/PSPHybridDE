#!/usr/bin/env python
# coding: utf-8

import glob
import os
import shutil
import sys
from datetime import date

MAIN_OUTPUT = "./results_Time_DefaultEvodock"
HPC = False
MPI = False

pdb_archive = [
    "/easy_dock/1oph_AB.prepack.pdb",
    "/easy_dock/1ml0_AD.prepack.pdb",
    "/easy_dock/1ktz_AB.prepack.pdb",
    "/easy_dock/1qa9_AB.prepack.pdb",
    # "/easy_dock/2hrk_AB_0001.prepack.pdb",
    "/easy_dock/2hrk_AB.prepack.pdb",
    "/easy_dock/1kxp_AD.prepack.pdb",
    "/easy_dock/2hle_AB.prepack.pdb",
    "/easy_dock/1b6c_AB.prepack.pdb",
    "/easy_dock/1ppe_IE.prepack.pdb",
]


def create_combinations(f_options, cr_options):
    combinations = []
    for f in f_options:
        for cr in cr_options:
            combinations.append((f, cr))
    return combinations


def create_output_file(file_name, f, cr, rep=0):
    # fname = file_name.split("/")[-1].split(".")[0].replace("options_", "")
    f = f.replace(".", "")
    cr = cr.replace(".", "")
    evofile = "evolution_f" + f + "_cr" + cr + ".log"
    if rep >= 0:
        os.makedirs(file_name + "/" + str(rep) + "/", exist_ok=True)
        fname = file_name + "/" + str(rep) + "/" + evofile
    else:
        fname = file_name + "/" + evofile
    return fname


def print_file(file_name, input_pdb, output_filename, scheme, np, max_iter, f, cr):
    with open(file_name, "w") as fout:
        fout.write("[inputs]\n")
        # fout.write("reference_input=sample_data/3QDEA.inv\n")
        fout.write("pose_input={}\n".format(input_pdb))
        fout.write("\n")
        fout.write("[outputs]\n")
        fout.write("output_file={}\n".format(output_filename))
        fout.write("\n")
        fout.write("[position]\n")
        fout.write("rot_max_magnitude=180\n")
        fout.write("trans_max_magnitude=70\n")
        fout.write("\n")
        fout.write("[DE]\n")
        fout.write("scheme={}\n".format(scheme))
        fout.write("popsize={}\n".format(np))
        fout.write("mutate={}\n".format(f))
        fout.write("recombination={}\n".format(cr))
        fout.write("maxiter={}\n".format(max_iter))
        fout.write("local_search=mcm_rosetta\n")


def create_filename(prot, f, cr, rep):
    f = f.replace(".", "")
    cr = cr.replace(".", "")
    if rep >= 0:
        name = "options_" + prot + "_f" + f + "_cr" + cr + "_" + str(rep) + ".ini"
    else:
        name = "options_" + prot + "_f" + f + "_cr" + cr + ".ini"
    return name


def select_input_pdb(prot):
    input_pdb = [x for x in pdb_archive if prot in x]
    if len(input_pdb) > 1:
        print("exit, found more than 1 pdb")
        exit()
    elif len(input_pdb) == 0:
        print("exit, pdb not found ", prot)
        exit()
    else:
        input_pdb = input_pdb[0]
    return input_pdb


def print_hpc_header(len_nodes):
    print("#!/bin/sh")
    print("#SBATCH -t 10:00:00")
    print("#SBATCH -o logs/output.log")
    print("#SBATCH -e logs/error.log")
    print("#SBATCH --mem-per-cpu=4000")
    print("#SBATCH -A SNIC2020-5-147")
    if MPI:
        print("#SBATCH -n {}".format(len_nodes * 14))
    else:
        print("#SBATCH -n 1")
    print("cd /home/d/dvarela/pfs/evodock/")
    print(
        "export PYTHONPATH=$PYTHONPATH:/home/d/dvarela/pfs/PyRosetta4.MinSizeRel.python37.linux.release-265/"
    )
    print("export PYTHONPATH=$PYTHONPATH:/home/d/dvarela/pfs/evodock/src/")


def print_for_hpc(filename):
    if MPI:
        line = "srun -N 1 -n 14 python evodock.py {} &> log.txt &".format(filename)
    else:
        line = "srun -N 1 -n 1 python evodock.py {} &> log.txt &".format(filename)
    print(line)


def create_benchmark(benchmark_name, combinations, max_iter, scheme, np, repetitions):
    # config_folder = "./configs/"
    config_folder = benchmark_name + "/configs/"
    os.makedirs(config_folder, exist_ok=True)
    # prots = ["1oph", "1ml0", "1ktz", "1ppe", "1b6c", "2hrk", "1kxp", "2hle", "1qa9"]
    # prots = ["1oph", "1ml0", "1ppe", "1b6c", "1kxp", "2hle", "1qa9"]
    # prots = ["2hle", "1qa9"]
    # prots = ["1oph", "1ml0", "1ktz"]
    # prots = ["1ppe", "1b6c", "2hrk"]
    # prots = ["1ml0"]
    # prots = ["2hrk", "1ml0", "2hle"]
    prots = ["1qa9", "1ktz"]
    if HPC:
        orig_stdout = sys.stdout
        f = open("run_bash.sh", "w")
        sys.stdout = f
        print_hpc_header(len(prots) * len(combinations) * len(repetitions))
    for p in prots:
        for c in combinations:
            for rep in repetitions:
                f = c[0]
                cr = c[1]
                file_name = create_filename(p, f, cr, rep)
                file_name = config_folder + file_name
                output_dir = benchmark_name + "/" + p
                os.makedirs(output_dir, exist_ok=True)
                output_filename = create_output_file(output_dir, f, cr, rep)
                input_pdb = select_input_pdb(p)
                print_file(
                    file_name, input_pdb, output_filename, scheme, np, max_iter, f, cr
                )
                if HPC:
                    print_for_hpc(file_name)
                else:
                    print("python evodock.py {} > log.txt &".format(file_name))
    if HPC:
        print("wait")
        sys.stdout = orig_stdout


def init_output_folder(results_folder):
    today = date.today()
    d1 = today.strftime("%d%m%Y")
    results_folder = d1 + "_" + results_folder
    results_folder = MAIN_OUTPUT + "/" + results_folder
    os.makedirs(results_folder, exist_ok=True)
    return results_folder


class CreateBenchmark:
    def create(self, params):
        benchmark_name = init_output_folder(params["filename"])
        cr_options = params["CR"]
        f_options = params["F"]
        combinations = create_combinations(f_options, cr_options)
        max_iter = params["MAX_ITER"]
        np = params["NP"]
        scheme = params["SCHEME"]
        repetitions = params["REP"]
        create_benchmark(
            benchmark_name, combinations, max_iter, scheme, np, repetitions
        )


def init_console(filename=""):
    params = {}
    params["CR"] = ["0.9"]
    params["F"] = ["0.2"]
    params["MAX_ITER"] = str(25)
    params["NP"] = str(25)
    params["SCHEME"] = "BEST"
    params["REP"] = 0

    if filename == "":
        params["filename"] = filename
    else:
        params["filename"] = ""
        while len(params["filename"]) < 1:
            params["filename"] = input(" insert a benchmark name > ").strip()
            if len(params["filename"]) < 1:
                print("benchmark name needed!")
            similar = [
                x for x in glob.glob(MAIN_OUTPUT + "/*") if params["filename"] in x
            ]
            if len(similar) > 0:
                prev = similar[0]
                print("benchmark name exists {} ".format(prev))
                remove = input(" remove existing? [0] > ").strip()
                remove = 0 if len(remove) < 1 else int(remove)
                if remove:
                    shutil.rmtree(prev)
                else:
                    params["filename"] = ""

        CR = input(" input CR [0.9] > ").strip()
        CR = CR.strip().split(";") if (";" in CR) else [CR.strip()]
        CR = ["0.9"] if len(CR) < 1 else CR
        params["CR"] = CR

        F = input(" input F [0.9] > ").strip()
        F = F.strip().split(";") if (";" in F) else [F.strip()]
        F = ["0.9"] if len(F) < 1 else F
        params["F"] = F

        MAX_ITER = input(" input MAX_ITER [25] > ").strip()
        MAX_ITER = MAX_ITER if len(MAX_ITER) > 1 else 25
        params["MAX_ITER"] = MAX_ITER

        NP = input(" input NP [25] > ").strip()
        NP = NP if len(NP) > 1 else 25
        params["NP"] = NP

        SCHEME = input(" input SCHEME [BEST] > ").strip()
        SCHEME = "BEST" if len(SCHEME) < 1 else SCHEME
        params["SCHEME"] = SCHEME

        REP = input(" repetition? [0] > ").strip()
        REP = list(range(0, int(REP))) if len(REP) > 1 else [-1]
        params["REP"] = REP

    return params


def main():
    if len(sys.argv) == 2:
        params = init_console()
    else:
        results_folder = sys.argv[-1]
        params = init_console(results_folder)

    CreateBenchmark().create(params)


if __name__ == "__main__":
    main()

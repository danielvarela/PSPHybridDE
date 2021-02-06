#!/usr/bin/env python
# coding: utf-8

import glob
import os
import sys

import numpy as np
import pandas as pd
from pyrosetta import Pose, init, pose_from_file
from pyrosetta.rosetta.protocols.moves import PyMOLMover

from proteins_information import prepack_pdb_archive as pdb_archive
from proteins_information import prot_parents
from run_rosetta_local_refinement import build_cmd
from scfxn_fullatom import FAFitnessFunction
from scfxn_zernike import l2_norm
from utils import get_position_info
from visualize_best_ind import visualize_single_individual


def init_options_fa(filename):
    opts = [
        "-mute all",
        # "-mute all -unmute core.pack.rotamers core.pack.rotamer_trials core.pack.task protocols",
        # "-out:level 1000",
        # "-docking:dock_mcm_trans_magnitude 0",
        "-docking:dock_mcm_first_cycles 1",
        "-docking:dock_mcm_second_cycles 1",
        # "-docking:dock_mcm_trans_magnitude 5",
        # "-docking:dock_mcm_rot_magnitude 0",
        "-include_current True",
        "-ex1",
        "-ex2aro",
        "-use_input_sc",
        "-extrachi_cutoff 1",
        "-unboundrot {}".format(filename),
    ]
    return " ".join(opts)


def get_data(files):
    experiments_info = {}
    for best_file in files:
        with open(best_file, "r") as f:
            file_id = int(best_file.split("/")[-2])
            lines = f.readlines()
            experiments_info[file_id] = [
                np.array(l.split(","), dtype=np.float64) for l in lines[1:]
            ]

    df_list = []
    for exp_id, data in experiments_info.items():
        translation_data = [d[3:] for d in data]
        translation_data_X = [d[0] for d in translation_data]
        translation_data_Y = [d[1] for d in translation_data]
        translation_data_Z = [d[2] for d in translation_data]
        df = pd.DataFrame(
            {
                "X": translation_data_X,
                "Y": translation_data_Y,
                "Z": translation_data_Z,
                "exp_id": exp_id,
            }
        )
        df_list.append(df)

    return df_list


def find_last_state(df_list, init_state):
    data = {}
    for df in df_list:
        last_x = df["X"].iloc[-1]
        last_y = df["Y"].iloc[-1]
        last_z = df["Z"].iloc[-1]
        data[df["exp_id"].iloc[-1]] = [last_x, last_y, last_z]

    values = []
    for k, v in data.items():
        l2 = np.linalg.norm(np.array(v) - np.array(init_state))
        values.append((k, l2, v))

    values.sort(key=lambda x: x[1], reverse=False)
    print(values[0])
    return values[0]


def send_best_ind_to_local_refinement(prot, best_pdb):
    print("local_refinement")
    os.makedirs("./rosetta_refinement_evodock/", exist_ok=True)
    os.makedirs("./rosetta_refinement_evodock/" + prot, exist_ok=True)
    os.makedirs("./rosetta_refinement_evodock/" + prot + "/input/", exist_ok=True)
    filename = "./rosetta_refinement_evodock/" + prot + "/input/evodock_best_ind.pdb"
    best_pdb.dump_pdb(filename)
    cmd = build_cmd(prot, "./input/evodock_best_ind.pdb", prot_parents[prot])
    cmd = "cd ./rosetta_refinement_evodock/" + prot + " && " + cmd + " && cd - "
    print(cmd)


def main():
    print("start find closest")
    folder = sys.argv[-1]
    best_files = glob.glob(folder + "*/best*.log")
    df_list = get_data(best_files)

    tokens = folder.split("/")
    keys = pdb_archive.keys()
    prot = ""
    for k in keys:
        if k in tokens:
            prot = k

    if prot != "":
        init(extra_options=init_options_fa(pdb_archive[prot]))

        IP_ADDRESS = "10.8.0.18"
        pymover = PyMOLMover(address=IP_ADDRESS, port=65000, max_packet_size=1400)
        native = Pose()
        pose_from_file(native, pdb_archive[prot])
        scfxn = FAFitnessFunction(native, native, 70)
        native_vec = get_position_info(native)
        positions = native_vec[3:]
        print("score native ", scfxn.scfxn_rosetta.score(native))
        best_ind, l2_value, best_positions = find_last_state(df_list, positions)
        print("value ", l2_value)
        print("best_positions")
        print(best_positions)
        best_vec = best_positions
        best_file = glob.glob(folder + "/" + str(best_ind) + "/best*.log")[-1]
        with open(best_file, "r") as f:
            lines = f.readlines()
            best_vec = np.array(lines[-1].split(","), dtype=np.float64)
        print(best_vec)
        pack_pose, _ = visualize_single_individual(pymover, scfxn, best_vec, best_ind)
        print("goal:")
        print(native_vec)
        pack_pose.dump_pdb("exit_from_GlobalHybridDE_" + prot + ".pdb")
        send_best_ind_to_local_refinement(prot, pack_pose)
    else:
        print(tokens)
        print("protein not found")


if __name__ == "__main__":
    main()

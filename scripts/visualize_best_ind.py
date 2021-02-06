#!/usr/bin/env python
# coding: utf-8

import sys
import time

import matplotlib.pyplot as plt
import numpy as np
from pyrosetta import Pose, init
from pyrosetta.rosetta.protocols.moves import PyMOLMover

from local_search import LocalSearchPopulation
from reader import StructureReader
from scfxn_fullatom import FAFitnessFunction
from scfxn_zernike import FitnessFunction

dict_references = {
    "3QDEA": "./sample_data/3QDEA.inv",
    "6U7IC": "./sample_data/input.inv",
}

dict_natives = {
    "3QDEA": "./sample_data/COMBINED_7.pdb",
    "6U7IC": "./sample_data/COMBINED_0.pdb",
}

pdb_archive = [
    "./easy_dock/1oph_AB.prepack.pdb",
    "./easy_dock/1ml0_AD.prepack.pdb",
    "./easy_dock/1ktz_AB.prepack.pdb",
    "./easy_dock/1qa9_AB.prepack.pdb",
    "./easy_dock/2hrk_AB.prepack.pdb",
    "./easy_dock/1kxp_AD.prepack.pdb",
    "./easy_dock/2hle_AB.prepack.pdb",
    "./easy_dock/1b6c_AB.prepack.pdb",
    "./easy_dock/1ppe_IE.prepack.pdb",
]

IP_ADDRESS = "10.8.0.10"


def get_individual_terms(scfxn, pose):
    valor = scfxn.scfxn_rosetta(pose)
    dict_scores = scfxn.get_dict_scores(pose)
    dict_scores = {
        k: float(v)
        for k, v in sorted(dict_scores.items(), key=lambda item: item[1], reverse=True)
    }
    return dict_scores


def plot_bar_difference(native_terms, best_terms, ind_id=0):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    x_axis = list(range(0, len(native_terms.keys())))
    width = 0.4
    p1 = ax.bar(x_axis, native_terms.values(), width, color="g")
    p2 = ax.bar([x + width for x in x_axis], best_terms.values(), width, color="r")
    ax.set_xticks([x + width for x in x_axis])
    ax.set_xticklabels(native_terms.keys(), rotation=45, ha="right")
    ax.legend((p1[0], p2[0]), ("native", "best"))
    plt.tight_layout()
    plt.savefig("test_barplot_" + str(ind_id) + ".png")


def visualize_single_individual(pymover, scfxn, vec, i):
    gen = scfxn.convert_positions_to_genotype(vec)
    tmp_pose = scfxn.apply_sixD_to_pose(gen)
    ls = LocalSearchPopulation(scfxn, "mcm_rosetta")
    # ls = LocalSearchPopulation(scfxn, "custom_packer")
    # ls = LocalSearchPopulation(scfxn, "custom_rotamer")

    # print(f"{bcolors.WARNING} START LOCAL SEARCH BENCHMARK {bcolors.ENDC}")
    start = time.time()
    energies = []
    for rep in list(range(0, 1)):
        pack_pose = Pose()
        pack_pose.assign(tmp_pose)
        ls.docking.apply(pack_pose)
        energies.append(scfxn.scfxn_rosetta(pack_pose))

    # pack_pose.dump_pdb("problematic_pdb_" + str(i) + ".pdb")
    end = time.time()
    e = np.array(energies)
    ls_time = end - start
    id_ = i
    print("|pdb_id|min|max|avg|std|time|")
    print(
        "|{}|{:.3f}|{:.3f}|{:.3f}|{:.3f}|{:.3f}|".format(
            id_, np.min(e), np.max(e), np.mean(e), np.std(e), ls_time
        )
    )

    ls.docking.apply(pack_pose)
    pack_pose.pdb_info().name("ls_best_pose_" + str(i))
    best_terms = get_individual_terms(scfxn, pack_pose)
    pymover.apply(pack_pose)
    return pack_pose, best_terms


def visualize_best_at_pymol(prot, best_evo):
    pymover = PyMOLMover(address=IP_ADDRESS, port=65000, max_packet_size=1400)
    file_native = [x for x in pdb_archive if prot in x.split("/")[-1].split("_")[0]][0]
    pose = StructureReader().get_from_file(file_native)
    native_pose = Pose()
    native_pose.assign(pose)

    # scfxn = FitnessFunction(file_native, native_pose, pose)
    scfxn = FAFitnessFunction(native_pose, pose, 50)
    native_pose.pdb_info().name("native_pose")
    pymover.apply(native_pose)
    native_terms = get_individual_terms(scfxn, native_pose)
    print("NATIVE SCORE ", str(scfxn.scfxn_rosetta.score(native_pose)))

    with open(best_evo, "r") as data_file:
        data = [x.strip().split(",") for x in data_file.readlines()[1:]]

    data = np.asarray(data, dtype="float")
    ls = LocalSearchPopulation(scfxn, "mcm_rosetta")

    # for i, vec in enumerate(data):
    #     # print(vec)
    #     gen = scfxn.convert_positions_to_genotype(vec)
    #     tmp_pose = scfxn.apply_sixD_to_pose(gen)
    #     tmp_pose.pdb_info().name("best_pose_" + str(i))
    #     ls.docking.apply(tmp_pose)
    #     best_terms = get_individual_terms(scfxn, tmp_pose)
    #     pymover.apply(tmp_pose)
    #     tmp_score = scfxn.scfxn_rosetta.score(tmp_pose)
    #     print("best gen {} - score {} ".format(i, tmp_score))

    # single
    vec = data[-1]
    ind_id = best_evo.split("/")[-2]
    pack_pose, best_terms = visualize_single_individual(pymover, scfxn, vec, ind_id)
    # plot_bar_difference(native_terms, best_terms, ind_id)


def main():
    best_evo = sys.argv[-1]
    prot = sys.argv[-2]
    ind_id = best_evo.split("/")[-2]
    # orig_stdout = sys.stdout
    # f = open("log_" + str(ind_id) + ".org", "w")
    # sys.stdout = f

    dict_options = {
        "-mute all",
        "-out:level 5000",
        "-unmute protocols.rigid.RigidBodyMover",
        "-docking:dock_mcm_first_cycles 1",
        "-docking:dock_mcm_second_cycles 1",
        # "-docking:dock_mcm_rot_magnitude 5",
        "-ex1",
        "-ex2aro",
    }

    str_options = " ".join(dict_options)

    init(extra_options=str_options)
    # init(
    #     extra_options="-mute all -zernike_descriptor:zernike_descriptor_file {}".format(
    #         dict_references[prot]
    #     )
    # )
    visualize_best_at_pymol(prot, best_evo)


if __name__ == "__main__":
    main()

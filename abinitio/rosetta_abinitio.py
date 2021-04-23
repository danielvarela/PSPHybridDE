#!/usr/bin/env python
# coding: utf-8


import configparser
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from pyrosetta import Pose, create_score_function, init
from pyrosetta.rosetta.core.scoring import CA_rmsd

from src.init_random_positions import start_input_poses
from src.local_search import AbinitioBuilder


def get_score_function(stage="None"):
    if stage == "stage1":
        scfxn = create_score_function("score0")
    if stage == "stage2":
        scfxn = create_score_function("score1")
    if stage == "stage3":
        scfxn = create_score_function("score2")
    if stage == "stage4":
        scfxn = create_score_function("score3")
    if stage == "None":
        scfxn = create_score_function("score3")
    return scfxn


def init_options():
    opts = [
        "-mute all",
        # "-out:level 500",
        # "-unmute protocols.abinitio",
        "-ignore_unrecognized_res",
        "-score::weights score3",
        "-in:file:centroid_input",
        "-nonideal true",
        "-corrections:restore_talaris_behavior",
        "-abinitio::rg_reweight 0.5",
        "-abinitio::rsd_wt_helix 0.5",
        "-abinitio::rsd_wt_loop 0.5",
        "-score:weights score3",
        "-abinitio::relax false",
        "-output_secondary_structure true",
        "-do_not_autoassign_SS true",
    ]
    return " ".join(opts)


def main():
    init(extra_options=init_options())

    config = configparser.ConfigParser()
    config.read(sys.argv[-1], encoding="utf-8-sig")

    stages = ["stage1", "stage2", "stage3", "stage4"]

    jobid = "./" + config["outputs"].get("output_file")
    native_input = config["inputs"].get("pose_input")
    native_pose, _ = start_input_poses(native_input, native_input)

    scfxn_rosetta = get_score_function()

    print("native score {} ".format(scfxn_rosetta.score(native_pose)))

    nstruct = config["DE"].getint("maxiter")
    results = []
    for n in range(nstruct):
        pose = Pose()
        pose.assign(native_pose)
        for stage in stages:
            scfxn_rosetta = get_score_function(stage)
            abinitio_builder = AbinitioBuilder(config, stage)
            abinitio = abinitio_builder.set_stage(stage)
            abinitio.apply(pose)
            pose = abinitio_builder.set_ss(pose)
            energy = scfxn_rosetta.score(pose)
            rmsd = CA_rmsd(native_pose, pose)
        results.append({"energy": energy, "rmsd": rmsd})

    df = pd.DataFrame(results)
    ax = sns.scatterplot(x="rmsd", y="energy", data=df, palette=["b"],)
    ax.set_title("Rosetta results")
    fig = ax.get_figure()
    fig.tight_layout()
    output_file = jobid.replace(".log", ".png")
    fig.savefig(output_file)
    plt.close()
    df.to_csv(jobid.replace(".log", ".csv"))


if __name__ == "__main__":
    main()

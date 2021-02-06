#!/usr/bin/env python
# coding: utf-8

import os
import sys

import imageio
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROSETTA_CSV = "~/projects/evodock/rosetta_results_database/rosetta_results_local_prepack_refinement_results.csv"


def adjust_lightness(color, amount=0.5):
    import matplotlib.colors as mc
    import colorsys

    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])


def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys

    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


def animation_evo(popul_file, output_file="evolution.gif"):
    rosetta_df = pd.read_csv(ROSETTA_CSV)
    prot_name = "2hle"
    rosetta_df = rosetta_df[rosetta_df["prot"] == prot_name]
    rosetta_df["gen"] = ["goal" for i in range(len(rosetta_df))]

    rosetta_df = rosetta_df.sample(1)
    rosetta_df = rosetta_df[["energy", "rmsd", "gen"]]
    rosetta_df.to_csv("goal_point.csv", index=False)

    with open(popul_file, "r") as data_file:
        data = [x.strip().split(",") for x in data_file.readlines()[1:]]
    data = np.asarray(data, dtype="float")
    frame = 0
    max_y = max([max(data[i]) for i in range(0, len(data), 2)])
    max_x = max([max(data[i]) for i in range(1, len(data), 2)])
    min_y = min([min(data[i]) for i in range(0, len(data), 2)])
    min_x = min([min(data[i]) for i in range(1, len(data), 2)])

    filenames = []
    # sns.color_palette("rocket_r", cmap="Blues")
    all_df = []
    for i in range(0, len(data), 2):
        df = pd.DataFrame()
        df["energy"] = pd.Series(data[i])
        df["rmsd"] = pd.Series(data[i + 1])
        df["gen"] = [len(all_df)] * len(df)
        df.to_csv(
            "./scripts/evolution_data/energy_rmsd_GEN" + str(len(all_df)) + ".csv",
            index=False,
        )
        all_df.append(df)

    df = pd.concat(all_df)
    df.to_csv("./scripts/evolution_data/energy_rmsd_ALL.csv", index=False)

    norm = mpl.colors.Normalize(vmin=0, vmax=18)
    cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Blues)
    cmap.set_array([])
    list_colors = [cmap.to_rgba(k) for k in range(18)]
    # print(list(list_colors))
    for i in range(0, max(df["gen"]), 1):
        # for i in [1, 4, 9, 17]:
        ax = sns.scatterplot(
            x="rmsd",
            y="energy",
            data=rosetta_df,
            hue="gen",
            palette=["r"],
            alpha=0.2,
            style="gen",
            markers=["s"],
        )

        ax = sns.scatterplot(
            x="rmsd",
            y="energy",
            data=df[df["gen"] == i],
            hue="gen",
            alpha=min(0.5 + 0.05 * i, 1),
            palette=["b"],
            # palette=[list_colors[i]],
            ax=ax,
        )
        ax.set_xlim([max(-1, -1), max_x + 2])
        ax.set_ylim([min_y - 10, max_y + 10])
        ax.set_title("GEN " + str(frame))
        ax.get_legend().remove()
        fig = ax.get_figure()
        fig.tight_layout()
        fig.savefig("frame_" + str(frame) + ".png")
        plt.close()
        filenames.append("frame_" + str(frame) + ".png")
        frame += 1

    images = []
    output_file = popul_file.replace("/", "_").replace(".log", "_ANIMATION.gif")
    for filename in filenames:
        images.append(imageio.imread(filename))
    imageio.mimsave(output_file, images, duration=2.0)
    for filename in filenames:
        os.remove(filename)


def main():
    animation_evo(sys.argv[-1])


if __name__ == "__main__":
    main()

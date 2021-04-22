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


def scatter_evo(popul_file, output_file="evolution.gif"):
    os.makedirs("./scripts/evolution_data/", exist_ok=True)
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
    ax = sns.scatterplot(
        x="rmsd",
        y="energy",
        data=df[df["gen"] == max(df["gen"])],
        hue="gen",
        alpha=min(0.5 + 0.05 * i, 1),
        palette=["b"],
    )
    ax.set_xlim([max(-1, -1), max_x + 2])
    ax.set_ylim([min_y - 10, max_y + 10])
    ax.set_title("GEN " + str(frame))
    ax.get_legend().remove()
    fig = ax.get_figure()
    fig.tight_layout()
    output_file = popul_file.replace("/", "_").replace(".log", "_SCATTER.png")
    fig.savefig(output_file)
    plt.close()


def main():
    scatter_evo(sys.argv[-1])


if __name__ == "__main__":
    main()

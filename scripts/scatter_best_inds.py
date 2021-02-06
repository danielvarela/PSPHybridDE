#!/usr/bin/env python
# coding: utf-8

import glob
import sys
from random import randint

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sb
from mpl_toolkits.mplot3d import Axes3D


def list_of_colors(n):
    color = []
    for i in range(n + 5):
        color.append("#%06X" % randint(0, 0xFFFFFF))
    return color


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
    colors = list_of_colors(len(experiments_info.keys()))
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
                "exp": [colors[exp_id] for x in range(len(translation_data_X))],
                "label": [exp_id for x in range(len(translation_data_X))],
            }
        )
        df_list.append(df)

    return df_list


def plot_different_experiments(df_list):
    plt.clf()
    fig = plt.figure()
    for count, df in enumerate(df_list):
        if count < 9:
            ax = fig.add_subplot(4, 3, count + 1, projection="3d")
            ax.scatter(df["X"], df["Y"], df["Z"], c=df["exp"])
            ax.scatter([0], [0], [0], c="green", marker="s")
            ax.set_xlim(-70, 70)
            ax.set_ylim(-70, 70)
            ax.set_zlim(-70, 70)
    fig.set_size_inches(12, 8)
    fig.tight_layout()
    fig.savefig("best_inds_3d.png")
    plt.close()


def plot_last_state(df_list):
    plt.clf()
    fig = plt.figure()
    last_x = []
    last_y = []
    last_z = []
    last_c = []
    for count, df in enumerate(df_list):
        last_x += [df["X"].iloc[-1]]
        last_y += [df["Y"].iloc[-1]]
        last_z += [df["Z"].iloc[-1]]
        last_c += [df["exp"].iloc[-1]]

    df = {"X": last_x, "Y": last_y, "Z": last_z, "exp": last_c}
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    ax.scatter(df["X"], df["Y"], df["Z"], c=df["exp"])
    ax.scatter([0], [0], [0], c="green", marker="s")
    ax.set_zlim(-70, 70)
    ax.set_xlim(-70, 70)
    ax.set_ylim(-70, 70)
    fig.tight_layout()
    fig.savefig("last_state_best_inds_3d.png")
    plt.close()


def main():
    print("start scatter plot 3D")

    folder = sys.argv[-1]
    best_files = glob.glob(folder + "*/best*.log")

    df_list = get_data(best_files)
    # plot
    plot_last_state(df_list)
    plot_different_experiments(df_list)


if __name__ == "__main__":
    main()

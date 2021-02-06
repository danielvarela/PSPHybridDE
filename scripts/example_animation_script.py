#!/usr/bin/env python
# coding: utf-8

import os

import imageio
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

datasets = {
    "goal": "./evolution_data/goal_point.csv",
    "individuals": "./evolution_data/energy_rmsd_ALL.csv",
}


def find_axes_limits(goal_dataset, individuals_dataset):
    offset = 5
    data = pd.concat([goal_dataset, individuals_dataset])
    min_x = data["rmsd"].min() - offset
    max_x = data["rmsd"].max() + offset
    min_y = data["energy"].min() - offset
    max_y = data["energy"].max() + offset
    return (min_x, max_x), (min_y, max_y)


def set_axes_limits(x_lim, y_lim, ax):
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    return ax


def set_goal_point(goal_dataset, ax):
    ax = sns.scatterplot(
        x="rmsd",
        y="energy",
        data=goal_dataset,
        hue="gen",
        palette=["r"],
        alpha=0.2,
        style="gen",
        markers=["s"],
        ax=ax,
    )
    return ax


def create_animation(framelist, output_file):
    images = []
    for filename in framelist:
        images.append(imageio.imread(filename))
    imageio.mimsave(output_file, images, duration=2.0)
    for filename in framelist:
        os.remove(filename)


def main():
    FRAMES_PATH = "./frames"
    OUTPUT_FILE = "ANIMATION.gif"
    goal_dataset = pd.read_csv(
        datasets["goal"], dtype={"energy": float, "rmsd": float, "gen": str}
    )
    individuals_dataset = pd.read_csv(
        datasets["individuals"], dtype={"energy": float, "rmsd": float, "gen": int}
    )

    fig, ax = plt.subplots(1, 1)

    x_lim, y_lim = find_axes_limits(goal_dataset, individuals_dataset)
    ax = set_axes_limits(x_lim, y_lim, ax)
    ax = set_goal_point(goal_dataset, ax)

    framelist = []
    os.makedirs(FRAMES_PATH, exist_ok=True)
    for gen in range(0, len(individuals_dataset["gen"].unique()), 1):
        ax = sns.scatterplot(
            x="rmsd",
            y="energy",
            data=individuals_dataset[individuals_dataset["gen"] == gen],
            hue="gen",
            palette=["b"],
            ax=ax,
        )
        ax.set_title("GEN " + str(gen))
        ax = set_axes_limits(x_lim, y_lim, ax)
        ax.get_legend().remove()
        fig = ax.get_figure()
        fig.tight_layout()
        frame_name = FRAMES_PATH + "/frame_" + str(gen) + ".png"
        fig.savefig(frame_name)
        ax.clear()
        plt.close()
        framelist.append(frame_name)

    create_animation(framelist, OUTPUT_FILE)
    os.rmdir(FRAMES_PATH)


if __name__ == "__main__":
    main()

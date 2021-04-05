#!/usr/bin/env python
# coding: utf-8

import glob
import shutil
import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

plt.rcParams.update({"font.size": 12})
pd.options.mode.chained_assignment = None  # default='warn'

def plots_during_evolution(df, name):
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(25, 10))
    ax = sns.lineplot(x=df.index, y="best", data=df, ax=axes[0])
    sns.lineplot(x=df.index, y="avg", data=df, ax=ax)
    sns.lineplot(x=df.index, y="rmsd", data=df, ax=axes[1])
    fig.tight_layout()
    fig_name = name.replace(".log", ".png")
    fig.savefig(fig_name)
    plt.close()
    return fig_name.replace(".png", "")


def print_table(df):
    print("** Summary table")
    min_egy = df["energy"].min()
    mean_egy = df["energy"].mean()
    # std_egy = df["energy"].std()
    min_rmsd = df["rmsd"].min()
    mean_rmsd = df["rmsd"].mean()
    # std_rmsd = df["rmsd"].std()
    min_I_sc = df["I_sc"].min()
    mean_I_sc = df["I_sc"].mean()
    min_irmsd = df["irmsd"].min()
    mean_irmsd = df["irmsd"].mean()

    print("|name|e|e|r|r|isc|isc|ir|ir|")
    print("|name|min|avg|min|avg|min|avg|min|avg|")
    print("|-")
    print(
        "|{:.2f}|{:.2f}|{:.2f}|{:.2f}|{:.2f}|{:.2f}|{:.2f}|{:.2f}|".format(
            min_egy,
            mean_egy,
            min_rmsd,
            mean_rmsd,
            min_I_sc,
            mean_I_sc,
            min_irmsd,
            mean_irmsd,
        )
    )


def plot_rmsd_vs_energy(prot, df_results, name):
    df = df_results[~df_results["name"].str.contains("rosetta")]
    df.sort_values(by=["rmsd"], ascending=True, inplace=True)
    sns.set_palette("deep")
    ax = sns.scatterplot(x="rmsd", y="energy", data=df, hue="name")
    ax.legend(bbox_to_anchor=(1.5, 1.0))
    fig = ax.get_figure()
    fig.set_size_inches(8.7, 5.0)
    fig.tight_layout()
    fig.savefig(name)
    print("** plot for all rmsd_energy results")
    print("[[./" + name + "]]")
    plt.close()


def print_info_experiment(evolution_file):
    energy_rmsd_results = []
    print("* Data for prot ")
    with open(evolution_file, "r") as evof:
        exp_data = []
        data = evof.readlines()
        if len(data) > 5:
            process_file = True
            data = data[2:]
            #print(data)
            for l in data:
                gen_info = [x for x in l.strip().split("\t")][1:]
                info = {
                    "gen": int(gen_info[0]),
                    "avg": float(gen_info[1]),
                    "best": float(gen_info[2]),
                    "rmsd": float(gen_info[3]),
                }
                exp_data.append(info)
        else:
            process_file = False
    if process_file:
        df = pd.DataFrame(exp_data)
        plots_during_evolution(df, evolution_file)
        name_legend = "prot" 
        energy_rmsd_results.append(
            {
                "name": name_legend,
                "energy": exp_data[-1]["best"],
                "rmsd": exp_data[-1]["rmsd"],
            }
        )

    df_results = pd.DataFrame(energy_rmsd_results)
    plot_rmsd_vs_energy("prot", df_results, "energy_rmsd.png")
    

def main():
    evolution_file = sys.argv[-1]
    name_png = evolution_file.replace(".log",".png")
    
    print_info_experiment(evolution_file)


if __name__ == "__main__":
    main()

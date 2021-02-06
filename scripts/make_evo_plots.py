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

ROSETTA_CSV = "~/projects/evodock/rosetta_results_database/rosetta_results_local_prepack_refinement_results.csv"
# ROSETTA_CSV = "~/projects/evodock/rosetta_results_database/rosetta_LOWRES_results.csv"

prots = ["1oph", "1ml0", "1ktz", "1ppe", "1b6c", "2hrk", "1kxp", "2hle", "1qa9", "1cgq"]


def print_table(df_results, name):
    df = df_results[~df_results["name"].str.contains("rosetta")]
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
        "|{}|{:.2f}|{:.2f}|{:.2f}|{:.2f}|{:.2f}|{:.2f}|{:.2f}|{:.2f}|".format(
            name,
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


def make_results_dataframe(prot_name, pfolder):
    rosetta_df = pd.read_csv(ROSETTA_CSV)
    rosetta_df = rosetta_df[rosetta_df["prot"] == prot_name]
    rosetta_df["name"] = ["rosetta" for i in range(len(rosetta_df))]

    seed_files = glob.glob(pfolder + "/*/seeds_*.log")
    seed_files += glob.glob(pfolder + "*/seeds_*.log")

    if len(seed_files) > 0:
        files = seed_files
    else:
        files = glob.glob(pfolder + "/*/popul_*.log")
        files += glob.glob(pfolder + "*/popul_*.log")
    df_list = []

    for popul_file in files:
        file_id = popul_file.split("/")[-2]
        with open(popul_file, "r") as data_file:
            data = [x.strip().split(",") for x in data_file.readlines()[1:]][-2:]
        if len(data) > 1:
            df = pd.DataFrame()
            df["energy"] = pd.Series([float(n) for n in data[0]])
            df["rmsd"] = pd.Series([float(n) for n in data[1]])
            df["name"] = pd.Series([str(prot_name + "_" + file_id) for x in data[0]])

            if len(seed_files) > 0:
                interface_file = (
                    "/".join(popul_file.split("/")[:-1]) + "/interface_f05_cr09.log"
                )
            else:
                interface_file = popul_file.replace("popul", "interface")
            with open(interface_file, "r") as data_file:
                data = [x.strip().split(",") for x in data_file.readlines()[1:]][-2:]
            df["I_sc"] = pd.Series([float(n) for n in data[0]])
            df["irmsd"] = pd.Series([float(n) for n in data[1]])
            df_list.append(df)

    df = pd.concat(df_list + [rosetta_df])
    return df


def plots_during_evolution(df, name):
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(25, 10))
    ax = sns.lineplot(x="gen", y="best", data=df, ax=axes[0])
    sns.lineplot(x="gen", y="avg", data=df, ax=ax)
    sns.lineplot(x="gen", y="rmsd", data=df, ax=axes[1])
    fig.tight_layout()
    fig_name = name.replace(".log", ".png")
    fig.savefig(fig_name)
    plt.close()
    return fig_name.replace(".png", "")


def print_evolution_figures_info(figure_names):
    for fig_name in figure_names:
        print("** plot for " + fig_name)
        print("[[./" + fig_name + ".png]]")


def print_last_gen_status_interface(prot_name, pfolder, df_results, ax=None):
    # rosetta_df = pd.read_csv(ROSETTA_CSV)
    # rosetta_df = rosetta_df[rosetta_df["prot"] == prot_name]
    # rosetta_df["name"] = ["rosetta" for i in range(len(rosetta_df))]
    df = df_results[~df_results["name"].str.contains("rosetta")]
    rosetta_df = df_results[df_results["name"].str.contains("rosetta")]

    df.sort_values(by=["irmsd"], ascending=True, inplace=True)
    sns.set_palette("deep")
    ax = sns.scatterplot(
        x="irmsd",
        y="I_sc",
        data=rosetta_df,
        hue="name",
        palette=["b"],
        alpha=0.2,
        style="name",
        markers=["s"],
    )

    ax = sns.scatterplot(x="irmsd", y="I_sc", data=df, hue="name", ax=ax)
    # ax = sns.scatterplot(x="irmsd", y="I_sc", data=rosetta_df, hue="name", ax=ax)
    ax.set_title("interface status " + prot_name)
    ax.legend(bbox_to_anchor=(1.2, 1.0))
    # fig = ax.get_figure()
    # fig.set_size_inches(8.7, 5.0)
    # fig.tight_layout()
    # fig_name = pfolder + "/interface_last_gen_" + prot_name + ".png"
    # fig.savefig(fig_name)
    # plt.close()
    # return fig_name


def print_last_gen_status(prot_name, pfolder, df_results, ax=None):
    # rosetta_df = pd.read_csv(ROSETTA_CSV)
    # rosetta_df = rosetta_df[rosetta_df["prot"] == prot_name]
    # rosetta_df["name"] = ["rosetta" for i in range(len(rosetta_df))]
    df = df_results[~df_results["name"].str.contains("rosetta")]
    rosetta_df = df_results[df_results["name"].str.contains("rosetta")]
    df.sort_values(by=["rmsd"], ascending=True, inplace=True)
    sns.set_palette("deep")
    ax = sns.scatterplot(
        x="rmsd",
        y="energy",
        data=rosetta_df,
        hue="name",
        palette=["black"],
        alpha=0.2,
        style="name",
        markers=["s"],
        ax=ax,
    )
    ax = sns.scatterplot(x="rmsd", y="energy", data=df, hue="name", ax=ax)
    # ax = sns.scatterplot(x="rmsd", y="energy", data=rosetta_df, hue="name", ax=ax)
    ax.set_title("last gen status for " + prot_name)
    ax.get_legend().remove()
    # ax.legend(bbox_to_anchor=(1.5, 1.0))


def plot_rmsd_vs_energy(prot, df_results, name):
    df = df_results[~df_results["name"].str.contains("rosetta")]
    rosetta_df = df_results[df_results["name"].str.contains("rosetta")]
    df.sort_values(by=["rmsd"], ascending=True, inplace=True)
    sns.set_palette("deep")
    ax = sns.scatterplot(
        x="rmsd",
        y="energy",
        data=rosetta_df,
        hue="name",
        palette=["black"],
        alpha=0.2,
        style="name",
        markers=["s"],
    )

    ax = sns.scatterplot(x="rmsd", y="energy", data=df, hue="name", ax=ax)
    ax.legend(bbox_to_anchor=(1.5, 1.0))
    fig = ax.get_figure()
    fig.set_size_inches(8.7, 5.0)
    fig.tight_layout()
    fig.savefig(name)
    print("** plot for all rmsd_energy results")
    print("[[./" + name + "]]")
    plt.close()


def print_info_experiment(results_folder, zip_archive):
    energy_rmsd_results = []
    prot_folders = glob.glob(results_folder + "/*")
    for pfolder in prot_folders:
        energy_rmsd_results = []
        figure_names = []
        prot = pfolder.split("/")[-1]
        all_results = make_results_dataframe(prot, pfolder)
        if prot in prots:
            print("* Data for {} ".format(prot))
            files = glob.glob(pfolder + "/evolution_*.log")
            files += glob.glob(pfolder + "/*/evolution_*.log")
            for count, f in enumerate(files):
                exp_data = []
                with open(f, "r") as evof:
                    data = evof.readlines()
                    if len(data) > 5:
                        process_file = True
                        data = data[2:]
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
                    figure_names.append(
                        (exp_data[-1]["rmsd"], plots_during_evolution(df, f))
                    )
                    name_legend = (
                        f.split("/")[-1].replace("evolution_", "").replace(".log", "")
                    )
                    name_legend = name_legend + " ({:.2f})".format(exp_data[-1]["rmsd"])
                    name_legend = str(count) + "_" + name_legend
                    energy_rmsd_results.append(
                        {
                            "name": name_legend,
                            "energy": exp_data[-1]["best"],
                            "rmsd": exp_data[-1]["rmsd"],
                        }
                    )

            df_results = pd.DataFrame(energy_rmsd_results)
            print_table(all_results, zip_archive)
            plot_rmsd_vs_energy(prot, df_results, pfolder + "/energy_rmsd.png")
            fig_name = join_global_interface_plots(prot, pfolder, all_results)
            print("** last gen status")
            print("[[./" + fig_name + "]]")
            figure_names.sort(key=lambda x: [0])
            figure_names = [f for (rmsd, f) in figure_names]
            print_evolution_figures_info(figure_names)

        shutil.make_archive(zip_archive, "zip", "./" + results_folder)


def join_global_interface_plots(prot, pfolder, all_results):
    plt.clf()
    plt.figure()
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
    print_last_gen_status(prot, pfolder, all_results, ax1)
    print_last_gen_status_interface(prot, pfolder, all_results, ax2)
    # fig.set_size_inches(8.7, 5.0)
    fig.set_size_inches(12, 5)
    fig.tight_layout()
    fig_name = pfolder + "/final_status_" + prot + ".png"
    fig.savefig(fig_name)
    plt.close()
    return fig_name


def main():
    results_folder = sys.argv[-1]
    if results_folder[-1] != "/":
        results_folder = results_folder + "/"
    zip_archive = results_folder.split("/")[-2]

    orig_stdout = sys.stdout
    f = open(zip_archive + ".org", "w")
    sys.stdout = f

    print_info_experiment(results_folder, zip_archive)


if __name__ == "__main__":
    main()

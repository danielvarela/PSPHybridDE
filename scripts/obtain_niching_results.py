#!/usr/bin/env python
# coding: utf-8

import glob
import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from proteins_information import prepack_pdb_archive

plt.rcParams.update({"font.size": 12})

ROSETTA_CSV = "~/projects/evodock/rosetta_results_database/rosetta_results_local_prepack_refinement_results.csv"
rosetta_df = pd.read_csv(ROSETTA_CSV)
rosetta_df["name"] = ["rosetta" for i in range(len(rosetta_df))]


def search_at_seeds_evolution(input_folder):
    print("* Information for {}".format(input_folder))
    for folder in glob.glob(input_folder + "/*"):
        figures = []
        if folder.split("/")[-1] in prepack_pdb_archive.keys():
            prot = folder.split("/")[-1]
            evolution_df = []
            seed_evo_files = glob.glob(folder + "/*/seeds_evolution*.log")
            for seed_evo_file in seed_evo_files:
                with open(seed_evo_file, "r") as data_file:
                    data = [x.strip().split(",") for x in data_file.readlines()[1:]]
                for i in range(0, len(data), 2):
                    df = pd.DataFrame()
                    df["energy"] = pd.Series(data[i], dtype="float")
                    df["rmsd"] = pd.Series(data[i + 1], dtype="float")
                    evolution_df.append(df)

            evolution_df = pd.concat(evolution_df)
            figures.append((prot, print_scatter(evolution_df, input_folder, prot)))

            min_row = evolution_df[
                evolution_df["energy"].eq(evolution_df["energy"].min())
            ].iloc[0]
            energy_min = min_row["energy"]
            rmsd_min_energy = min_row["rmsd"]
            print("best for {} : {} ({}) ".format(prot, energy_min, rmsd_min_energy))

        for f in figures:
            print("** scatter final ", f[0])
            print("[[{}]]".format(f[1]))


def print_scatter(df, pfolder, prot):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    prot_rosetta_df = rosetta_df[rosetta_df["prot"] == prot]
    df.sort_values(by=["rmsd"], ascending=True, inplace=True)
    sns.set_palette("deep")
    ax = sns.scatterplot(
        x="rmsd",
        y="energy",
        data=prot_rosetta_df,
        hue="name",
        palette=["black"],
        alpha=0.2,
        style="name",
        markers=["s"],
        ax=ax,
    )
    ax = sns.scatterplot(x="rmsd", y="energy", data=df, ax=ax)
    # ax = sns.scatterplot(x="rmsd", y="energy", data=rosetta_df, hue="name", ax=ax)
    ax.set_title("last gen status for " + prot)
    ax.get_legend().remove()
    fig.set_size_inches(12, 5)
    fig.tight_layout()
    fig_name = pfolder + "/final_status_" + prot + ".png"
    fig.savefig(fig_name)
    plt.close()
    return fig_name


def get_evolution_information(input_folder):
    print("* Information for {}".format(input_folder))
    return_data = []
    for folder in glob.glob(input_folder + "/*"):
        prot = [x for x in folder.split("/") if x in prepack_pdb_archive.keys()]
        if len(prot) > 0:
            prot = prot[0]
            evolution_df = []
            seed_evo_files = glob.glob(folder + "/*/evolution*.log")
            seed_evo_files += glob.glob(folder + "*/*/evolution*.log")
            for seed_evo_file in seed_evo_files:
                with open(seed_evo_file, "r") as data_file:
                    data = [
                        x.strip().split(" ")[-2:] for x in data_file.readlines()[2:]
                    ]

                df = pd.DataFrame()
                df["energy"] = pd.Series([x[0] for x in data], dtype="float")
                df["rmsd"] = pd.Series([x[1] for x in data], dtype="float")
                evolution_df.append(df)

            evolution_df = pd.concat(evolution_df)

            min_row = evolution_df[
                evolution_df["energy"].eq(evolution_df["energy"].min())
            ].iloc[0]
            energy_min = min_row["energy"]
            rmsd_min_energy = min_row["rmsd"]
            print(
                "best for {} : {:.2f} ({:.2f}) ".format(
                    prot, energy_min, rmsd_min_energy
                )
            )
            return_data.append(
                {
                    "exp": input_folder,
                    "prot": prot,
                    "energy": energy_min,
                    "rmsd": rmsd_min_energy,
                }
            )
    return return_data


def main():
    input_folder = sys.argv[-1]
    search_at_seeds_evolution(input_folder)
    # folders_to_compare = [
    #     "results_SeedsDE/06102020_species_F01CR09",
    #     "results_SeedsDE/06102020_species_F05CR09",
    #     "results_SeedsDE/06102020_species_F09CR03",
    #     "results_Crowding/05102020_crowding05_F09CR03",
    #     "results_Crowding/05102020_crowding100_F09CR03",
    #     "results_Crowding/05102020_firstTry_F09CR03",
    #     "results_FixSeedsDE/08102020_fixrecalculate_F09CR09",
    #     "results_kb_stage2/24092020_MCMrosetta_F09CR03",
    #     "results_mpi/29092020_highsearch_F09CR03",
    #     "results_stage2/24092020_MCMrosetta_F09CR03",
    #     "results_SeedsDE/13102020_superDiversity_F05CR09",
    # ]

    # data = []
    # for input_folder in folders_to_compare:
    #     return_data = get_evolution_information(input_folder)
    #     data += return_data

    # # print(data)
    # for prot in ["2hrk"]:
    #     print("* Data for ", prot)
    #     print("| exp | energy | rmsd |")
    #     for row in data:
    #         if row["prot"] == prot:
    #             print("|{}|{}|{}|".format(row["exp"], row["energy"], row["rmsd"]))


if __name__ == "__main__":
    main()

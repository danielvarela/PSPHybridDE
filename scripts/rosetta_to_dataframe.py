#!/usr/bin/env python
# coding: utf-8

import glob
import sys

import pandas as pd


def read_results(prot_info):
    score_file = prot_info[1] + "/score.sc"
    with open(score_file, "r") as sc:
        lines = sc.readlines()[2:]
        lines = [l.split() for l in lines]

    energies = [float(l[1]) for l in lines]
    rmsds = [float(l[2]) for l in lines]
    iscore = [float(l[5]) for l in lines]
    irmsds = [float(l[6]) for l in lines]
    filenames = [l[-1] for l in lines]
    prot = [prot_info[0] for l in lines]
    df = pd.DataFrame(
        {
            "energy": energies,
            "rmsd": rmsds,
            "filename": filenames,
            "I_sc": iscore,
            "irmsd": irmsds,
            "prot": prot,
        }
    )
    return df


def main():
    input_folder = sys.argv[-1]
    prot_folders = [(f.split("/")[-1], f) for f in glob.glob(input_folder + "/*")]
    df_list = []
    for info in prot_folders:
        df_list.append(read_results(info))
    df = pd.concat(df_list)
    df.to_csv("rosetta_results_" + input_folder.split("/")[-1] + ".csv", index=False)


if __name__ == "__main__":
    main()

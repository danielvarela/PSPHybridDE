#!/usr/bin/env python
# coding: utf-8

import os
import shutil
import sys

from proteins_information import prepack_pdb_archive as pdb_archive
from proteins_information import prot_parents


def build_cmd(prot, filename, partners):
    cmd = " ".join(
        [
            "$ROS/docking_protocol.linuxgccrelease -overwrite",
            "-s {}".format(filename),
            "-nstruct 50",
            "-partners {}".format(partners),
            "-docking_local_refine",
            "-ex1 -ex2aro",
            "-use_input_sc",
            "-unboundrot {}".format(filename),
        ]
    )

    return cmd


def main():
    print("Rosetta local refinement")
    prots = ["1oph", "1ml0", "1ktz", "1ppe", "1b6c", "2hrk", "1kxp", "2hle", "1qa9"]
    list_of_cmd = []

    destiny_folder = "./local_prepack_refinement_results/"
    os.makedirs(destiny_folder, exist_ok=True)
    for p in prots:
        os.makedirs(destiny_folder + p, exist_ok=True)
        os.makedirs(destiny_folder + p + "/easy_dock/", exist_ok=True)
        shutil.copy(pdb_archive[p], destiny_folder + p + "/easy_dock/")
        cmd = build_cmd(p, pdb_archive[p], prot_parents[p])
        cmd = "cd " + destiny_folder + p + " && " + cmd + " && cd - "
        list_of_cmd.append(cmd)

    f = open("run_local_refinement.sh", "w")
    sys.stdout = f
    print("#!/bin/bash")
    print("")
    print("\n".join([c + " & " for c in list_of_cmd]))


if __name__ == "__main__":
    main()

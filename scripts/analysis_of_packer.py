#!/usr/bin/env python
# coding: utf-8

import os
import sys

from proteins_information import pdb_archive, prot_parents
from repacking import run_repacking


def create_options(filename, prot, option="standard"):
    opts = " ".join(
        [
            "-mute all",
            # "-unmute core.pack.pack_rotamers core.pack.dunbrack",
            "-partners {}".format(prot_parents[prot]),
        ]
    )

    extra_opt = " -ex1 -ex2aro"
    extrachi_opt = " -extrachi_cutoff 1"
    inputsc_opt = " -use_input_sc"
    unbound_opt = " -unboundrot {}".format(filename)

    if option == "extra" or option == "all":
        opts += extra_opt
    if option == "chi" or option == "all":
        opts += extrachi_opt
    if option == "inputsc" or option == "all":
        opts += inputsc_opt
    if option == "current" or option == "all":
        opts += " -include_current True"
    else:
        opts += " -include_current False"
    if option == "unbound" or option == "all":
        opts += unbound_opt
    print(opts)
    return opts


def main():
    # prots = ["2hrk"]
    prots = ["1oph", "1ml0", "1ktz", "1ppe", "1b6c", "2hrk", "1kxp", "2hle", "1qa9"]
    packers = ["standard", "extra", "current", "chi", "input_sc", "unbound", "all"]
    # packers = ["all"]

    os.makedirs("./packer_pdbs/", exist_ok=True)
    for prot in prots:
        os.makedirs("./packer_pdbs/" + prot, exist_ok=True)
        energy_results = {}
        for packer_config in packers:
            filename = pdb_archive[prot]
            opts = create_options(filename, prot, packer_config)
            final_energy, packed_pdb = run_repacking(filename, opts)
            packed_pdb.dump_pdb(
                "./packer_pdbs/" + prot + "/" + prot + "." + packer_config + ".pdb"
            )
            energy_results[packer_config] = final_energy

        orig_stdout = sys.stdout
        f = open("packer_analysis.org", "a")
        sys.stdout = f
        print("** prot ", prot)
        print("|prot|" + " | ".join(["{}".format(p) for p in packers]) + "|")
        print("|-")
        print(
            "|"
            + prot
            + "|"
            + " | ".join(["{:.2f}".format(v) for k, v in energy_results.items()])
            + "|"
        )
        sys.stdout = orig_stdout


if __name__ == "__main__":
    main()

#!/usr/bin/env python
# coding: utf-8

import glob
import sys

from pyrosetta import Pose, init, pose_from_file
from pyrosetta.rosetta.core.scoring import ScoreFunction, ScoreType
from pyrosetta.rosetta.protocols.motifs import name3_from_oneletter
from pyrosetta.rosetta.protocols.moves import PyMOLMover


def main():
    opts = " ".join(
        ["-mute all", "-partners A_B", "-ex1", "-ex2aro", "-extrachi_cutoff 1"]
    )
    init(extra_options=opts)

    pdb_folder = sys.argv[-1]
    scorefxn_low = ScoreFunction()
    scorefxn_low.set_weight(ScoreType.fa_rep, 0.55)

    pymover = PyMOLMover(address="10.8.0.6", port=65000, max_packet_size=1400)
    f = open("packer_analysis_problematic_residues.org", "a")
    sys.stdout = f
    print("* Complex ", [d for d in pdb_folder.split("/") if len(d) > 1][-1])
    for pdb in glob.glob(pdb_folder + "/*.pdb"):
        config = pdb.split(".")[-2]
        pose = Pose()
        pose_from_file(pose, pdb)
        pymover.apply(pose)
        sequence = pose.sequence()
        scorefxn_low.score(pose)
        energies = pose.energies()
        residue_energies = [
            (i, energies.residue_total_energy(i))
            for i in range(1, pose.total_residue() + 1)
        ]

        print("** problematic for ", config)
        for d in residue_energies:
            if d[1] > 10:
                res_name = name3_from_oneletter(sequence[d[0]])
                print("*** {}_{} {:.2f}".format(res_name, d[0], d[1]))


if __name__ == "__main__":
    main()

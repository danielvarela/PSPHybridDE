#!usr/bin/env python

import sys

from pyrosetta import Pose, Vector1, init, standard_packer_task
from pyrosetta.rosetta.core.pack.task import PackerTask, TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import (
    IncludeCurrent, InitializeFromCommandline, NoRepackDisulfides,
    RestrictToRepacking)
from pyrosetta.rosetta.protocols.loops import get_fa_scorefxn
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.protocols.moves import PyMOLMover

from local_search import LocalSearchPopulation
from reader import StructureReader
from scfxn_fullatom import FAFitnessFunction


class PosePacker:
    def __init__(self, filename):
        self.input_pose = StructureReader().get_fa_from_file(filename)

        self.output_name = filename.replace(".pdb", ".prepack.pdb")

    def run(self):
        pack_pdb = self.pack(self.input_pose)
        pack_pdb.dump_pdb(self.output_name)
        return pack_pdb

    def pack(self, pose):
        pack_pose = Pose()
        pack_pose.assign(pose)

        scorefxn = get_fa_scorefxn()
        pose_packer = standard_packer_task(pack_pose)

        pose_packer.restrict_to_repacking()
        pose_packer.or_include_current(False)

        # local_tf = TaskFactory()
        # local_tf.push_back(InitializeFromCommandline())
        # local_tf.push_back(IncludeCurrent())
        # local_tf.push_back(RestrictToRepacking())
        # local_tf.push_back(NoRepackDisulfides())
        # conformer_full_repack = PackRotamersMover(scorefxn)
        # conformer_full_repack.task_factory(local_tf)
        # packmover = PackRotamersMover(scorefxn, pose_packer)
        # # print("Pre packing score : {} ".format(scorefxn(pose)))
        packmover.apply(pack_pose)
        # conformer_full_repack.apply(pack_pose)
        # print("Post packing score : {} ".format(scorefxn(pack_pose)))
        return pack_pose


def main():
    # pymover = PyMOLMover(address="10.8.0.14", port=65000, max_packet_size=1400)
    pose_input = sys.argv[-1]
    init(
        extra_options="-mute all -unmute core.pack.pack_rotamers core.pack.dunbrack -ex1 -ex2aro -partners A_B -unboundrot {}".format(
            pose_input
        )
    )

    pose = StructureReader().get_fa_from_file(pose_input)
    pose.pdb_info().name("INIT_STATE")
    scfxn = FAFitnessFunction(pose, pose, 50)
    print("INIT ENERGY ", scfxn.scfxn_rosetta.score(pose))
    # pymover.apply(pose)
    packer = PosePacker(pose_input)

    ls = LocalSearchPopulation(scfxn, "custom_packer")
    pack_pose = Pose()
    pack_pose.assign(pose)
    # ls.docking.apply(pack_pose)
    # print("FINAL ENERGY ", scfxn.scfxn_rosetta.score(pack_pose))
    packed_pdb = packer.run()
    print("FINAL ENERGY ", scfxn.scfxn_rosetta.score(packed_pdb))
    pack_pose = Pose()
    pack_pose.assign(packed_pdb)
    # ls.docking.apply(pack_pose)

    # packed_pdb.pdb_info().name("PACKED")

    # pymover.apply(packed_pdb)


if __name__ == "__main__":
    main()

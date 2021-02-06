import logging
import time

from pyrosetta import (MoveMap, Pose, SwitchResidueTypeSetMover,
                       create_score_function, pose_from_file)
from pyrosetta.rosetta.core.fragment import ConstantLengthFragSet
from pyrosetta.rosetta.core.pose import read_psipred_ss2_file
from pyrosetta.rosetta.core.scoring import CA_rmsd
from pyrosetta.rosetta.protocols.abinitio import ClassicAbinitio

from differential_evolution import Individual


class AbinitioBuilder:
    def __init__(self, config, stage="None"):
        if config.has_option("DE", "local_search"):
            local_search_option = config["DE"].get("local_search")
        else:
            local_search_option = "None"
        self.cycles = 0.01
        self.local_stage = local_search_option
        short_frag_filename = config["inputs"].get("frag3_input")
        long_frag_filename = config["inputs"].get("frag9_input")
        self.native = Pose()
        pose_from_file(self.native, config["inputs"].get("pose_input"))
        to_centroid = SwitchResidueTypeSetMover("centroid")
        to_centroid.apply(self.native)
        self.fragset_long = ConstantLengthFragSet(9, long_frag_filename)
        self.fragset_short = ConstantLengthFragSet(3, short_frag_filename)
        ss2_file = config["inputs"].get("ss2_input")
        read_psipred_ss2_file(self.native, ss2_file)
        self.movemap = MoveMap()
        self.movemap.set_bb(True)
        self.abinitio = self.set_stage(local_search_option)

    def init_abinitio(self):
        abinitio = ClassicAbinitio(self.fragset_short, self.fragset_long, self.movemap)
        abinitio.set_cycles(self.cycles)
        abinitio.bSkipStage1_ = True
        abinitio.bSkipStage2_ = True
        abinitio.bSkipStage3_ = True
        abinitio.bSkipStage4_ = True
        return abinitio

    def set_stage(self, stage):
        abinitio = self.init_abinitio()
        abinitio.init(self.native)
        abinitio.bSkipStage1_ = True
        abinitio.bSkipStage2_ = True
        abinitio.bSkipStage3_ = True
        abinitio.bSkipStage4_ = True
        if stage == "stage1":
            abinitio.bSkipStage1_ = False
        if stage == "stage2":
            abinitio.bSkipStage2_ = False
        if stage == "stage3":
            abinitio.bSkipStage3_ = False
        if stage == "stage4":
            abinitio.bSkipStage4_ = False
        abinitio.set_cycles(self.cycles)
        return abinitio


class LocalSearchPopulation:
    # Options:
    # None: only score and return the poses
    # only_slide: just slide_into_contact
    # custom_rotamer: RotamerTrialsMover
    # custom_packer: stochastic default packer from rosetta
    # mcm_rosetta: mcm protocol mover (high res) from rosetta (2 cycles)
    def __init__(self, scfxn, config, abinitio_builder):
        if config.has_option("DE", "local_search"):
            local_search_option = config["DE"].get("local_search")
        else:
            local_search_option = "None"

        self.local_stage = local_search_option
        self.abinitio_builder = abinitio_builder
        self.scfxn = scfxn
        self.abinitio = self.abinitio_builder.abinitio

    def energy_score(self, pose):
        score = self.scfxn.scfxn_rosetta(pose)
        return score

    def process_individual(self, ind, local_search=False):
        pose = self.scfxn.convert_genotype_to_ind_pose(ind)
        before = self.energy_score(pose)
        if self.local_stage != "None":
            self.abinitio.apply(pose)
            after = self.energy_score(pose)
        else:
            after = before
        rmsd = self.scfxn.get_rmsd(pose)
        interface = 0
        irms = 0
        genotype = self.scfxn.convert_positions_to_genotype(self.scfxn.get_dofs(pose))
        result_individual = Individual(genotype, after, rmsd, interface, irms)
        return result_individual, before, after

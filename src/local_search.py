import logging
import time

from pyrosetta import (MoveMap, Pose, SwitchResidueTypeSetMover,
                       create_score_function, pose_from_file)
from pyrosetta.rosetta.core.fragment import ConstantLengthFragSet
from pyrosetta.rosetta.core.pose import read_psipred_ss2_file
from pyrosetta.rosetta.core.scoring import CA_rmsd
from pyrosetta.rosetta.protocols.abinitio import ClassicAbinitio
from pyrosetta.rosetta.protocols.membrane import get_secstruct
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.protocols.simple_moves import ShearMover, SmallMover

from differential_evolution import Individual


def transform_ss(code):
    if code == "C":
        return "L"
    if code == "E":
        return "E"
    if code == "H":
        return "H"


def get_ss_from_file(file_ss):
    with open(file_ss, "r") as f:
        lines = [l.strip() for l in f.readlines()[2:] if len(l) > 0]
        codes = [transform_ss(l.split(" ")[2]) for l in lines]
    return codes


def get_ss(pose):
    ss = []
    for i in range(1, pose.total_residue() + 1):
        ss.append(pose.secstruct(i))
    return ss


class AbinitioBuilder:
    def __init__(self, config, stage="None"):
        if config.has_option("DE", "local_search"):
            local_search_option = config["DE"].get("local_search")
        else:
            local_search_option = "None"
        # self.cycles = 0.01
        if config.has_option("DE", "cycles"):
            self.cycles = float(config["DE"].get("cycles"))
        else:
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
        self.ss2_struct = get_ss_from_file(ss2_file)
        self.movemap = MoveMap()
        self.movemap.set_bb(True)
        self.abinitio = self.set_stage(local_search_option)

    def set_ss(self, pose):
        for i in range(1, pose.total_residue() + 1):
            pose.set_secstruct(i, self.ss2_struct[i - 1])
        return pose

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


class MCMminimization:
    def __init__(self, scfxn):
        self.pose = scfxn.native_pose
        kT = 1.0
        n_moves = 50
        movemap = MoveMap()
        movemap.set_bb(True)
        self.small_mover = SmallMover(movemap, kT, n_moves)
        self.small_mover.scorefxn(scfxn.scfxn_rosetta)
        self.shear_mover = ShearMover(movemap, kT, n_moves)
        self.shear_mover.scorefxn(scfxn.scfxn_rosetta)

    def apply(self, pose):
        self.small_mover.apply(pose)
        self.shear_mover.apply(pose)
        return pose


class GradientMinimizer:
    def __init__(self, scfxn):
        kT = 1.0
        n_moves = 50
        movemap = MoveMap()
        movemap.set_bb(True)
        self.small_mover = SmallMover(movemap, kT, n_moves)
        self.small_mover.scorefxn(scfxn.scfxn_rosetta)
        self.shear_mover = ShearMover(movemap, kT, n_moves)
        self.shear_mover.scorefxn(scfxn.scfxn_rosetta)
        self.min_mover = MinMover()
        self.min_mover.movemap(movemap)
        self.min_mover.score_function(scfxn.scfxn_rosetta)

    def apply(self, pose):
        self.small_mover.apply(pose)
        self.shear_mover.apply(pose)
        self.min_mover.apply(pose)
        return pose


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
        pose = self.scfxn.convert_ind_to_pose(ind)
        backup = Pose()
        backup.assign(pose)
        before = self.energy_score(pose)
        if self.local_stage != "None":
            self.abinitio.apply(pose)
            pose = self.abinitio_builder.set_ss(pose)
            after = self.energy_score(pose)
        else:
            after = before

        if after > before:
            pose.assign(backup)
            after = before

        rmsd = self.scfxn.get_rmsd(pose)
        omegas = [pose.omega(i) for i in range(1, len(pose) + 1)]
        ss = [pose.secstruct(i) for i in range(1, len(pose) + 1)]
        genotype = self.scfxn.convert_positions_to_genotype(self.scfxn.get_dofs(pose))
        result_individual = Individual(genotype, after, rmsd, omegas, ss)
        return result_individual, before, after

#!usr/bin/env python

import logging
from math import sqrt

import numpy as np
from pyrosetta import Pose, create_score_function
from pyrosetta.rosetta.core.scoring import CA_rmsd

from differential_evolution import Individual
from utils import convert_range


def l2_norm(zd_reference, zd_current):
    L2_norm = 0
    for i, v in enumerate(zd_current):
        L2_norm += (zd_current[i] - zd_reference[i]) * (zd_current[i] - zd_reference[i])
    sqrt_L2_norm = sqrt(L2_norm)
    return sqrt_L2_norm


class PSPFitnessFunction:
    def __init__(self, stage, reference, native_pose, input_pose):
        self.logger = logging.getLogger("evodock.scfxn")
        self.native_pose = native_pose
        self.logger.setLevel(logging.INFO)
        # self.pymover = PyMOLMover(address=IP_ADDRESS, port=65000, max_packet_size=1400)
        self.scfxn_rosetta = self.get_score_function(stage)
        self.dock_pose = Pose()
        self.dock_pose.assign(input_pose)
        self.dock_pose.pdb_info().name("INIT_STATE")

    def get_score_function(self, stage):
        if stage == "stage1":
            scfxn = create_score_function("score0")
        if stage == "stage2":
            scfxn = create_score_function("score1")
        if stage == "stage3":
            scfxn = create_score_function("score2")
        if stage == "stage4":
            scfxn = create_score_function("score3")
        if stage == "None":
            scfxn = create_score_function("score3")
        return scfxn

    def get_native_ss(self):
        ss = [
            self.native_pose.secstruct(i) for i in range(1, len(self.native_pose) + 1)
        ]
        return ss

    def get_native_omegas(self):
        omegas = [
            self.native_pose.omega(i) for i in range(1, len(self.native_pose) + 1)
        ]
        return omegas

    def get_dofs(self, pose):
        phis = [pose.phi(i) for i in range(1, len(pose) + 1)]
        psis = [pose.psi(i) for i in range(1, len(pose) + 1)]
        return phis + psis

    def get_rmsd(self, pose):
        rmsd = CA_rmsd(self.native_pose, pose)
        return rmsd

    def render_models(self, pdb_id, dofs, is_best=None):
        pose = self.apply_dofs_to_pose(dofs)
        dst = self.scfxn_rosetta.score(pose)
        if np.isnan(dst):
            dst = 10000
        prot_name = "popul" if is_best is None else is_best
        pose.pdb_info().name(prot_name + "_pose_" + str(pdb_id))
        # self.pymover.apply(pose)
        rmsd = self.get_rmsd(pose)
        return dst, rmsd

    def score(self, dofs_vector):
        pose = self.apply_dofs_to_pose(dofs_vector)
        try:
            dst = self.scfxn_rosetta.score(pose)
        except ValueError:
            dst = 10000
        if np.isnan(dst):
            dst = 10000
        return dst

    def size(self):
        return self.native_pose.total_residue() * 2

    def get_sol_string(self, sol):
        return " , ".join([str(e) for e in sol])

    def convert_positions_to_genotype(self, positions):
        # return self.converter.convert_positions_to_genotype(positions)
        pos = []
        for p in positions:
            pos.append(convert_range(p, (-180, 180), (-1, 1)))
        return pos

    def convert_genotype_to_positions(self, genotype):
        gen = []
        for i, g in enumerate(genotype):
            gen.append(convert_range(g, (-1, 1), (-180, 180)))
        return gen
        # return convert_range(genotype)

    def apply_dofs_to_pose(self, genotype):
        ind = Individual(
            genotype, 0, 1000, self.get_native_omegas(), self.get_native_ss()
        )
        return self.convert_ind_to_pose(ind)

    def convert_ind_to_pose(self, ind):
        dofs = self.convert_genotype_to_positions(ind.genotype)
        phis = dofs[: int(len(dofs) / 2)]
        psis = dofs[int(len(dofs) / 2) :]
        ind_pose = Pose()
        ind_pose.assign(self.dock_pose)
        for i in range(1, len(phis) + 1):
            ind_pose.set_phi(i, phis[i - 1])
            ind_pose.set_psi(i, psis[i - 1])
            ind_pose.set_omega(i, ind.omegas[i - 1])
            ind_pose.set_secstruct(i, ind.ss[i - 1])

        # now is time to score the join_pose
        return ind_pose

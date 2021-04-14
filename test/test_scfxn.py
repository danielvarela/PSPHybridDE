#!/usr/bin/env python
# coding: utf-8

from pyrosetta import create_score_function, init

from init_random_positions import start_input_poses
from scfxn_psp import PSPFitnessFunction


def init_options(ss2):
    opts = [
        "-mute all",
        "-in:file:centroid_input",
        "-ignore_unrecognized_res",
        "-nonideal true",
        "-out:pdb",
        "-use_filters true",
        "-corrections:restore_talaris_behavior",
        "-abinitio::rg_reweight 0.5",
        "-abinitio::rsd_wt_helix 0.5",
        "-abinitio::rsd_wt_loop 0.5",
        "-psipred_ss2 {}".format(ss2),
        "-score:weights score3",
        "-abinitio::relax false",
        "-output_secondary_structure true",
        "-do_not_autoassign_SS true",
    ]
    return " ".join(opts)


def main():
    ss2_file = "./input_files/info_1c8c/vf_1c8cA.psipred_ss2"
    init(init_options(ss2_file))
    pose_input = "./input_files/1c8ca_mejor_score.pdb"

    native_pose, init_state_pose = start_input_poses(pose_input, pose_input)

    score3 = create_score_function("score3")

    print("score3 : ", score3.score(native_pose))
    scfxn = PSPFitnessFunction("stage4", native_pose, native_pose, native_pose)
    dofs = scfxn.get_dofs(native_pose)
    native_genotype = scfxn.convert_positions_to_genotype(dofs)

    print("len pose {} ".format(native_pose.total_residue()))
    print("dofs ", dofs[0:10])
    print("native_genotype {} ".format(len(native_genotype)))
    print("native_genotype ", native_genotype[0:10])
    print("my score result : ", scfxn.score(native_genotype))


if __name__ == "__main__":
    main()

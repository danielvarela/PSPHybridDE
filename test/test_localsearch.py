#!/usr/bin/env python
# coding: utf-8

import configparser

from pyrosetta import create_score_function, init

from differential_evolution import Individual
from init_random_positions import start_input_poses
from scfxn_psp import PSPFitnessFunction
from src.local_search import AbinitioBuilder
from src.population import ScorePopulation


def read_config(ini_file):
    config = configparser.ConfigParser()
    config.read(ini_file, encoding="utf-8-sig")
    return config


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

    score3_value = score3.score(native_pose)
    print("score3 : ", score3_value)
    scfxn = PSPFitnessFunction("stage4", native_pose, native_pose, native_pose)
    dofs = scfxn.get_dofs(native_pose)
    native_genotype = scfxn.convert_positions_to_genotype(dofs)

    print("len pose {} ".format(native_pose.total_residue()))
    print("dofs ", dofs[0:10])
    print("native_genotype {} ".format(len(native_genotype)))
    print("native_genotype ", native_genotype[0:10])
    print("my score result : ", scfxn.score(native_genotype))

    config = read_config("./configs/sample_psp_1c8cA.ini")
    stage = "stage4"
    abinitio_builder = AbinitioBuilder(config, stage)
    score_popul = ScorePopulation(
        scfxn, "./evolution_test.log", abinitio_builder, config
    )

    native_ind = Individual(
        native_genotype, 0, 1000, scfxn.get_native_omegas(), scfxn.get_native_ss()
    )

    new_ind, before, after = score_popul.local_search.process_individual(native_ind)

    min_value = after
    min_dofs = scfxn.convert_genotype_to_positions(new_ind.genotype)
    for i in range(1000):
        new_ind, before, after = score_popul.local_search.process_individual(new_ind)
        if abs(before - min_value) > 0.01:
            print("not correct before to min_value")
            print("{} vs {}".format(before, min_value))
            current_dofs = scfxn.convert_genotype_to_positions(new_ind.genotype)
            for d in range(len(current_dofs)):
                if current_dofs[d] != min_dofs[d]:
                    print("different dofs at {}".format(d))
                    print("{} {} ".format(current_dofs[d], min_dofs[d]))
                    exit()
            exit()
        if after < min_value:
            print("new {} -> {} ".format(min_value, after))
            min_dofs = scfxn.convert_genotype_to_positions(new_ind.genotype)
            min_value = after
        if before < after:
            print("before {} after {} ".format(before, after))

        # if (i % 100) == 0:
        #     print("final after: {}".format(after))

    print("final after: {}".format(after))


if __name__ == "__main__":
    main()

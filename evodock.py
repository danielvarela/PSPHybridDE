#!/usr/bin/env python
# coding: utf-8

import configparser
import logging
import os
import sys

import numpy as np
# from mpi4py import MPI
from pyrosetta import init

from differential_evolution import DifferentialEvolutionAlgorithm as DE
from init_random_positions import start_input_poses
from local_search import AbinitioBuilder
from population import ScorePopulation
from scfxn_psp import PSPFitnessFunction
from single_process import SingleMasterProcess as MasterProcess

# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# size = comm.Get_size()


MAIN_PATH = os.getcwd()

logging.basicConfig(level=logging.ERROR)
# logging.disable(logging.INFO)


def init_options(reference_input):
    opts = [
        "-mute all",
        "-ignore_unrecognized_res",
        "-nonideal true",
        "-corrections:restore_talaris_behavior",
        "-abinitio::rg_reweight 0.5",
        "-abinitio::rsd_wt_helix 0.5",
        "-abinitio::rsd_wt_loop 0.5",
        "-score:weights score3",
        "-abinitio::relax false",
        "-output_secondary_structure true",
        "-do_not_autoassign_SS true",
    ]
    return " ".join(opts)


def read_config(ini_file):
    config = configparser.ConfigParser()
    config.read(ini_file, encoding="utf-8-sig")
    return config


def get_translation_max(native_pose, dock_pose):
    jump_num = 1
    flexible_jump = dock_pose.jump(jump_num)
    translation = np.asarray(flexible_jump.get_translation())
    return max(translation) + 70


def main():
    config = read_config(sys.argv[-1])
    if config.has_option("inputs", "reference_input"):
        # reference_input = MAIN_PATH + "/" + config["inputs"].get("reference_input")
        reference_input = config["inputs"].get("reference_input")
    else:
        reference_input = ""

    # pose_input = MAIN_PATH + "/" + config["inputs"].get("pose_input")
    pose_input = config["inputs"].get("pose_input")

    init(extra_options=init_options(reference_input))

    logger = logging.getLogger("evodock")
    logger.setLevel(logging.INFO)

    # --- DE PARAMS -----------------------------------+

    scheme = config["DE"].get("scheme")
    popsize = config["DE"].getint("popsize")
    mutate = config["DE"].getfloat("mutate")
    recombination = config["DE"].getfloat("recombination")
    maxiter = config["DE"].getint("maxiter")

    if config.has_option("DE", "local_search"):
        local_search_option = config["DE"].get("local_search")
    else:
        logger.info("DANGER: local_search is None")
        local_search_option = "None"

    # --- INIT -----------------------------------------+
    jobid = "./" + config["outputs"].get("output_file")
    native_input = config["inputs"].get("pose_input")
    native_pose, init_state_pose = start_input_poses(pose_input, native_input)

    rank = 0
    size = 1
    # -- TESTING PURPOUSES ----
    stages = ["stage1", "stage2", "stage3", "stage4"]

    for stage in stages:
        abinitio_builder = AbinitioBuilder(config, stage)
        scfxn = PSPFitnessFunction(stage, reference_input, native_pose, init_state_pose)

        logger.info(
            " init_state pose {:.2f}".format(scfxn.scfxn_rosetta.score(init_state_pose))
        )

        score_popul = ScorePopulation(scfxn, jobid, abinitio_builder, config)
        master_calculator = MasterProcess(size, score_popul)

        if rank == 0:
            logger.info("==============================")
            logger.info(" init the params ")
            logger.info(" starts stage {}".format(stage))
            if stage == "stage1":
                alg = DE(
                    master_calculator,
                    scheme,
                    popsize,
                    mutate,
                    recombination,
                    int(maxiter / (len(stages) - 1)),
                    jobid,
                )
            else:
                alg.set_popul_calculator(master_calculator)

            # --- RUN -----------------------------------------+
            logger.info("==============================")
            logger.info(" starts high res evolution")
            logger.info(
                " native pose {:.2f}".format(scfxn.scfxn_rosetta.score(native_pose))
            )
            init_population = alg.init_population()
            high_res_population = master_calculator.run(init_population)
            if stage != "stage1":
                init_population = alg.main(high_res_population)
            else:
                init_population = high_res_population
            if stage == "stage4":
                master_calculator.terminate()
        # else:
        #     Worker(rank, size, score_popul).run()


if __name__ == "__main__":
    main()

import copy
import time

import numpy as np
from mpi4py import MPI

from differential_evolution import Individual
from init_random_positions import start_input_poses
from local_search import AbinitioBuilder
from population import ScorePopulation
from scfxn_psp import PSPFitnessFunction

comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# size = comm.Get_size()


class IndividualMPI:
    def __init__(self, idx, ind):
        self.idx = idx
        if type(ind) is list:
            len_ind = int(len(ind) / 2)
            self.genotype = ind
            self.score = 1000
            self.rmsd = 0
            self.omegas = [0 for f in range(len(len_ind))]
            self.ss = ["L" for f in range(len(len_ind))]
        else:
            self.genotype = ind.genotype
            self.score = ind.score
            self.rmsd = ind.rmsd
            self.omegas = ind.omegas
            self.ss = ind.ss

    def convert_to_np(self):
        a = self.genotype
        a.append(-10000)
        a.append(self.score)
        a.append(-20000)
        a.append(self.rmsd)
        a.append(-30000)
        a = np.concatenate((a, np.array(self.omegas)), axis=0)
        a = np.concatenate((a, np.array([-40000])), axis=0)
        a = np.concatenate((a, np.array(convert_ss_to_code(self.ss))), axis=0)
        a = np.concatenate((a, np.array([-50000])), axis=0)
        a = np.concatenate((a, np.array([-self.idx])), axis=0)
        a = np.ravel(a)
        return a

    def convert_to_ind(self):
        return Individual(self.genotype, self.score, self.rmsd, self.omegas, self.ss)


def convert_code_to_ss(code):
    ss = []
    code = [int(round(c)) for c in code]
    for c in code:
        if c == 0:
            ss.append("L")
        if c == 1:
            ss.append("H")
        if c == 2:
            ss.append("E")
    return ss


def convert_ss_to_code(ss):
    code = []
    for s in ss:
        if s == "L":
            code.append(0)
        if s == "H":
            code.append(1)
        if s == "E":
            code.append(2)
    return code


def np_to_ind(a):
    a = list(a)
    genotype_idx = a.index(-10000)
    score_idx = a.index(-20000)
    rmsd_idx = a.index(-30000)
    omegas_idx = a.index(-40000)
    ss_idx = a.index(-50000)

    genotype = a[:genotype_idx]
    score = a[score_idx - 1]
    rmsd = a[rmsd_idx - 1]
    omegas = a[rmsd_idx + 1 : omegas_idx]
    ss = a[omegas_idx + 1 : ss_idx]
    ss = convert_code_to_ss(ss)

    ind = Individual(genotype, score, rmsd, omegas, ss)

    idx = int(round(a[-1]))
    return (idx, ind)


class MasterProcess:
    def __init__(self, size, cost_func, stage):
        self.rank = 0
        self.size = size
        self.cost_func = cost_func
        self.stage = int(stage.replace("stage", ""))

    def __make_chunks(self, popul, size):
        lst = []
        for d in range(0, len(popul)):
            lst.append(IndividualMPI(d, popul[d]).convert_to_np())
        size = size
        chunks = [lst[i : i + size] for i in range(0, len(lst), size)]
        return chunks

    def run(self, popul, with_slide=False):
        list_popul = popul
        req = []
        ind_per_process = int(len(list_popul) / self.size)
        data = self.__make_chunks(list_popul, ind_per_process)
        ind_size = len(data[0][0])
        recv_data = np.zeros(
            shape=(self.size, ind_per_process, ind_size), dtype=np.float64
        )
        run_score_function = self.stage
        for d in range(1, self.size):
            comm.send(run_score_function, dest=d)
            comm.send(data[d], dest=d)
            r = comm.irecv(buf=recv_data[d], source=d, tag=d)
            req.append(r)

        convert_pop = [np_to_ind(x) for x in data[0]]
        mpi_pop = []
        for idx, ind in convert_pop:
            scored_ind, before, after = self.cost_func.local_search.process_individual(
                ind
            )
            mpi_pop.append(IndividualMPI(idx, scored_ind).convert_to_np())

        for idx, ind in enumerate(mpi_pop):
            recv_data[0][idx] = ind

        all_finish = False
        while all_finish is False:
            check_finish = []
            for r in req:
                s = r.Test()
                check_finish.append(s)
            all_finish = all([s is True for s in check_finish])

        # print("recv data ", recv_data)
        result_pop = popul
        for proc in recv_data:
            for arr in proc:
                idx, ind = np_to_ind(arr)
                result_pop[idx] = ind

        return result_pop

    def terminate(self):
        terminate_function = -1
        for d in range(1, comm.Get_size()):
            comm.send(terminate_function, dest=d)


class Worker:
    def __init__(self, rank, size, config):
        self.rank = rank
        self.size = size
        jobid = "./" + config["outputs"].get("output_file")
        native_input = config["inputs"].get("pose_input")
        pose_input = config["inputs"].get("pose_input")
        reference_input = config["inputs"].get("reference_input")
        native_pose, init_state_pose = start_input_poses(pose_input, native_input)
        self.score_stages = {}
        for stage in ["stage1", "stage2", "stage3", "stage4"]:
            abinitio_builder = AbinitioBuilder(config, stage)
            scfxn = PSPFitnessFunction(
                stage, reference_input, native_pose, init_state_pose
            )
            score_popul = ScorePopulation(scfxn, jobid, abinitio_builder, config)
            self.score_stages[stage] = score_popul

        self.cost_func = self.score_stages["stage1"]

    def run(self):
        # functions = {1: "run_score", 2: "run_score_rmsd", 3: "terminate"}
        while True:
            action = comm.recv(source=0)
            if action == -1:
                exit()
            else:
                if action < 5:
                    self.cost_func = self.score_stages["stage{}".format(action)]
            data = comm.recv(source=0)
            # process data here
            convert_pop = [np_to_ind(x) for x in data]
            mpi_pop = []
            for idx, ind in convert_pop:
                (
                    scored_ind,
                    before,
                    after,
                ) = self.cost_func.local_search.process_individual(ind)
                mpi_pop.append(IndividualMPI(idx, scored_ind).convert_to_np())

            send_data = np.array(mpi_pop, dtype=np.float64)
            comm.Send(send_data, dest=0, tag=self.rank)

import numpy as np

from mpi_utils import IndividualMPI, np_to_ind


class SingleMasterProcess:
    def __init__(self, size, cost_func):
        self.rank = 0
        self.size = size
        self.cost_func = cost_func

    def __make_chunks(self, popul, size):
        lst = []
        for d in range(0, len(popul)):
            lst.append(IndividualMPI(d, popul[d]).convert_to_np())
        size = size
        chunks = [lst[i : i + size] for i in range(0, len(lst), size)]
        return chunks

    def run(self, popul, with_slide=False):
        list_popul = popul
        ind_per_process = int(len(list_popul) / self.size)
        data = self.__make_chunks(list_popul, ind_per_process)
        ind_size = len(data[0][0])
        recv_data = np.zeros(
            shape=(self.size, ind_per_process, ind_size), dtype=np.float64
        )
        convert_pop = [np_to_ind(x) for x in data[0]]
        mpi_pop = []
        for idx, ind in convert_pop:
            if ind.score == 1000:
                scored_ind, _, _ = self.cost_func.local_search.process_individual(
                    ind.genotype, with_slide
                )
                mpi_pop.append(IndividualMPI(idx, scored_ind).convert_to_np())
            else:
                mpi_pop.append(IndividualMPI(idx, ind).convert_to_np())

        for idx, ind in enumerate(mpi_pop):
            recv_data[0][idx] = ind

        result_pop = popul
        for proc in recv_data:
            for arr in proc:
                idx, ind = np_to_ind(arr)
                result_pop[idx] = ind

        # result_pop = self.cost_func.local_search.set_new_max_translations(result_pop)
        return result_pop

    def terminate(self):
        # terminate
        pass

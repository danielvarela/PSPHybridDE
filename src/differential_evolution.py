import logging
import random
import time

from selection import (CrowdingSelection, EliteSelection, GreedySelection,
                       nearest_neighbor)


def ensure_bounds(vec, bounds):
    vec_new = []
    # cycle through each variable in vector
    for i, k in enumerate(vec):

        # variable exceedes the minimum boundary
        if vec[i] < bounds[i][0]:
            vec_new.append(bounds[i][0])

        # variable exceedes the maximum boundary
        if vec[i] > bounds[i][1]:
            vec_new.append(bounds[i][1])

        # the variable is fine
        if bounds[i][0] <= vec[i] <= bounds[i][1]:
            vec_new.append(vec[i])

    return vec_new


class Individual:
    def __init__(self, genotype, score, rmsd=1000, omegas=[], ss=[]):
        self.genotype = genotype
        self.score = score
        self.rmsd = rmsd
        self.omegas = omegas
        self.ss = ss


# --- MAIN ---------------------------------------------------------------------+


class DifferentialEvolutionAlgorithm:
    def __init__(
        self,
        popul_calculator,
        scheme,
        popsize,
        mutate,
        recombination,
        maxiter,
        jobid,
        selection_opt,
    ):
        self.scheme = scheme
        self.logger = logging.getLogger("evodock.de")
        self.logger.setLevel(logging.INFO)
        self.popul_calculator = popul_calculator
        self.popsize = popsize
        self.mutate = mutate
        self.recombination = recombination
        self.maxiter = maxiter
        self.job_id = jobid
        self.ind_size = popul_calculator.cost_func.size()
        self.bounds = [(-1, 1)] * self.ind_size
        self.file_time_name = self.job_id.replace("evolution", "time")
        self.init_file()

        if selection_opt == "Greedy":
            self.selection = GreedySelection()
        if selection_opt == "Crowding":
            self.selection = CrowdingSelection(self.popul_calculator.cost_func)

    def set_popul_calculator(self, popul_calculator):
        self.popul_calculator = popul_calculator

    def init_population(self, popsize=None):
        # --- INITIALIZE A POPULATION (step #1) ----------------+
        if popsize is None:
            popsize = self.popsize
        self.logger.info(" init population")
        population = []
        native_omega = self.popul_calculator.cost_func.scfxn.get_native_omegas()
        native_ss = self.popul_calculator.cost_func.scfxn.get_native_ss()
        for i in range(0, popsize):
            indv = []
            for j in range(len(self.bounds)):
                indv.append(random.uniform(self.bounds[j][0], self.bounds[j][1]))
            population.append(Individual(indv, 0, 1000, native_omega, native_ss))

        population = self.popul_calculator.run(population, False)

        # scores = [p.score for p in population]
        # print(scores)
        # rmsds = [p.rmsd for p in population]
        # df = pd.DataFrame({"score": scores, "rmsd": rmsds})

        # df.plot.scatter(x="rmsd", y="score")
        # plt.savefig("init_population.png")

        return population

    def evaluate_population(self, popul):
        with_slide = False
        population = self.popul_calculator.run(popul, with_slide)
        return population

    def init_file(self):
        with open(self.job_id, "w") as file_object:
            file_object.write(
                "CONF: maxiter : {}, np : {}, f {}, cr {} \n".format(
                    self.maxiter, self.popsize, self.mutate, self.recombination
                )
            )
            file_object.write("INIT EVOLUTION\n")
        with open(self.file_time_name, "w") as file_time:
            file_time.write(
                "CONF: maxiter : {}, np : {}, f {}, cr {} \n".format(
                    self.maxiter, self.popsize, self.mutate, self.recombination
                )
            )
            file_time.write("INIT EVOLUTION\n")

    def main(self, population):

        self.logger.info(" DE")
        # --- SOLVE --------------------------------------------+

        self.popul_calculator.cost_func.print_information(population)
        # _ = input("continue > ")
        # cycle through each generation (step #2)
        for i in range(1, self.maxiter + 1):
            self.logger.info(" GENERATION:" + str(i))
            start = time.time()
            file_object = open(self.job_id, "a")
            file_time = open(self.file_time_name, "a")

            # file_object = open(outdir + "evolution_example.txt", "a")
            file_object.write("GENERATION: \t" + str(i) + "\t")

            # before_sc = [ind.score for ind in population]
            if (i % 10) == 0:
                population = self.popul_calculator.run(population)
            gen_scores = [ind.score for ind in population]
            # for idx in range(len(before_sc)):
            #     if (abs(before_sc[idx]) - abs(gen_scores[idx])) < -0.02:
            #         print(idx)
            #         print("problem with frags")
            #         print("bf {} af {} ".format(before_sc[idx], gen_scores[idx]))
            #         print(population[idx].genotype)
            #         exit()
            # cycle through each individual in the population
            trials = []
            for j in range(0, self.popsize):

                # --- MUTATION (step #3.A) ---------------------+

                # select three random vector index positions [0, self.popsize),
                # not including current vector (j)
                candidates = list(range(0, self.popsize))
                candidates.remove(j)
                random_index = random.sample(candidates, 3)

                best_index = gen_scores.index(min(gen_scores))
                if self.scheme == "CURRENT":
                    x_1 = population[j].genotype
                if self.scheme == "RANDOM":
                    x_1 = population[random_index[0]].genotype
                    x_1_idx = random_index[0]
                if self.scheme == "BEST":
                    x_1 = population[best_index].genotype
                    x_1_idx = best_index

                x_2 = population[random_index[1]].genotype
                x_3 = population[random_index[2]].genotype
                x_t = population[j].genotype  # target individual

                # subtract x3 from x2, and create a new vector (x_diff)
                x_diff = [x_2_i - x_3_i for x_2_i, x_3_i in zip(x_2, x_3)]

                # multiply x_diff by the mutation factor (F) and add to x_1

                v_donor = [
                    x_1_i + self.mutate * x_diff_i
                    for x_1_i, x_diff_i in zip(x_1, x_diff)
                ]

                v_donor = ensure_bounds(v_donor, self.bounds)

                # --- RECOMBINATION (step #3.B) ----------------+

                v_trial = []
                for k, obj in enumerate(x_t):
                    crossover = random.random()
                    if crossover <= self.recombination:
                        v_trial.append(v_donor[k])
                    else:
                        v_trial.append(x_t[k])

                trial_ind = Individual(
                    v_trial,
                    0,
                    1000,
                    population[x_1_idx].omegas,
                    population[x_1_idx].ss,
                )
                trials.append(trial_ind)

            # --- SELECTION (step #3.C) -------------+
            trial_inds = self.popul_calculator.run(trials)
            self.popul_calculator.cost_func.print_information(trial_inds, True)

            population, gen_scores, trial_scores = self.selection.apply(
                trial_inds, population
            )

            if any([ind.score == 0 for ind in population]):
                print("hay score == 0")
                exit()
            # --- SCORE KEEPING --------------------------------+
            # self.logger.info("POPUL SCORES : ", gen_scores)
            gen_avg = sum(gen_scores) / self.popsize  # current generation avg. fitness
            gen_best = min(gen_scores)  # fitness of best individual

            trial_avg = sum(trial_scores) / self.popsize
            trial_best = min(trial_scores)

            gen_sol = population[gen_scores.index(min(gen_scores))]

            self.logger.info("   > GENERATION AVERAGE: %f " % gen_avg)
            self.logger.info("   > GENERATION BEST: %f " % gen_best)
            self.logger.info(
                "   > TRIAL INFO: {:.2f} {:.2f} ".format(trial_best, trial_avg)
            )

            # problematic_terms = self.popul_calculator.cost_func.get_problematic_terms(trial_inds)
            # self.logger.info("   > TOP PROBLEM TERMS {}".format(str(problematic_terms)))

            best_SixD_vector, best_rmsd = self.popul_calculator.cost_func.render_best(
                i, gen_sol, population
            )
            best_sol_str = self.popul_calculator.cost_func.scfxn.get_sol_string(
                best_SixD_vector
            )
            # self.logger.info("   > BEST SOL: {} ".format(best_sol_str))
            self.popul_calculator.cost_func.print_information(population)
            # self.popul_calculator.cost_func.pymol_visualization(population)
            # _ = input("continue > ")
            file_object.write("%f \t" % gen_avg)
            file_object.write("%f \t" % gen_best)
            file_object.write("%f \n" % best_rmsd)
            file_object.close()
            end = time.time()
            file_time.write("%f \n" % (end - start))
            file_time.close()

        # self.popul_calculator.cost_func.print_popul(population)
        return population

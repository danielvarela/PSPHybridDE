import os

from pyrosetta.rosetta.protocols.moves import PyMOLMover

from local_search import LocalSearchPopulation
from utils import IP_ADDRESS


class IndividualZD:
    def __init__(self, SixD_vector=[], score=0, ZD=[]):
        self.SixD_vector = SixD_vector
        self.ZD = ZD
        self.score = score


class ScorePopulation:
    def __init__(self, scfxn, jobid, abinitio_builder, config=None):
        self.name = "ScorePopulation"
        self.scfxn = scfxn
        jobid_folder = "/".join([x for x in jobid.split("/")[:-1]])
        os.makedirs("./" + jobid_folder, exist_ok=True)
        self.log_best = jobid.replace("evolution", "best")
        self.log_popul = jobid.replace("evolution", "popul")
        self.log_trials = jobid.replace("evolution", "trials")
        self.local_search = LocalSearchPopulation(scfxn, config, abinitio_builder)
        with open(self.log_best, "w") as file_object:
            file_object.write("#{}\n".format(jobid))
        with open(self.log_trials, "w") as file_object:
            file_object.write("#{}\n".format(jobid))
        with open(self.log_popul, "w") as file_object:
            file_object.write("#{}\n".format(jobid))

    def score(self, popul):
        score_popul = []
        for ind in popul:
            score = self.scfxn.score(ind.SixD_vector)
            score_popul.append(IndividualZD(ind.SixD_vector, score, ind.ZD))
        return score_popul

    def convert_ind_to_pose(self, ind):
        return self.scfxn.convert_ind_to_pose(ind)

    def convert_genotype(self, genotype):
        return self.scfxn.convert_genotype_to_positions(genotype)

    def apply(self, genotypes):
        print("APPLY NO SE USA")
        exit()
        popul = []
        for g in genotypes:
            gen = self.convert_genotype(g)
            popul.append(IndividualZD(gen))
        popul_score = self.score(popul)
        return [p.score for p in popul_score]

    def render_best(self, gen, best_solution, population):
        SixD_vector = self.convert_genotype(best_solution.genotype)
        rmsd = best_solution.rmsd
        with open(self.log_best, "a") as file_object:
            # file_object.write("gen:\t{}\t".format(gen))
            vector_str = ",".join(["{}".format(i) for i in SixD_vector])
            file_object.write("{}\n".format(vector_str))

        return SixD_vector, rmsd

    def size(self):
        return self.scfxn.size()

    def popul_to_pdbs(self, popul):
        os.makedirs("./folder_pdbs/", exist_ok=True)
        for i, ind in enumerate(popul):
            # pose = self.scfxn.apply_sixD_to_pose(ind.SixD_vector)
            pose = ind.pose
            pose.dump_pdb("./folder_pdbs/ind_" + str(i) + ".pdb")

    def apply_local_search(self, popul):
        new_popul = self.local_search.apply(popul)
        return new_popul

    # -- debug --#
    def get_problematic_terms(self, popul):
        terms = {}
        for ind in popul:
            # pose = self.scfxn.apply_sixD_to_pose(ind)
            pose = ind.pose
            self.scfxn.scfxn_rosetta(pose)
            dict_scores = self.scfxn.get_dict_scores(pose)
            dict_scores = {
                k: v
                for k, v in sorted(
                    dict_scores.items(), key=lambda item: item[1], reverse=True
                )
            }

            problematic_terms = list(dict_scores.keys())[:3]
            for t in problematic_terms:
                if t in list(terms.keys()):
                    terms[t] += dict_scores[t]
                else:
                    terms[t] = dict_scores[t]
        avg_terms = {}
        for t, k in terms.items():
            avg_terms[t] = k / len(popul)
        return avg_terms

    # --- only for debug --- #

    def print_popul_info(self, popul, destiny, trial_popul=False):
        popul_dst = [ind.score for ind in popul]
        popul_rmsd = [ind.rmsd for ind in popul]

        with open(destiny, "a") as file_object:
            popul_dst_str = ",".join(["{:.2f}".format(i) for i in popul_dst])
            popul_rmsd_str = ",".join(["{:.2f}".format(i) for i in popul_rmsd])
            file_object.write("{}\n".format(popul_dst_str))
            file_object.write("{}\n".format(popul_rmsd_str))

    def print_information(self, popul, trial_popul=False):
        if trial_popul is False:
            destiny = self.log_popul
        else:
            destiny = self.log_trials
        self.print_popul_info(popul, destiny, trial_popul)

    def pymol_visualization(self, popul):
        pymover = PyMOLMover(address=IP_ADDRESS, port=65000, max_packet_size=1400)
        for ind, p in enumerate(popul):
            # gen = scfxn.convert_positions_to_genotype(p.genotype)
            gen = p.genotype
            tmp_pose = self.scfxn.apply_sixD_to_pose(gen)
            tmp_pose.pdb_info().name("popul_" + str(ind))
            pymover.apply(tmp_pose)

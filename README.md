

# Protein structure prediction using HybridDE

A hybrid version between differential evolution and the fragment replacement technique was defined for protein structure prediction [Varela2020]. The coarse-grained atomic model of the Rosetta system was used for protein representation. The high-dimensional and multimodal nature of protein energy landscapes requires an efficient search for obtaining the native structures with minimum energy. However, the energy model of Rosetta presents an additional difficulty, since the best energy area in the landscape does not necessarily correspond to the closest conformations to the native structure. A strategy is to obtain a diverse set of protein conformations that correspond to different minima in the landscape. The incorporation of the crowding niching method into the hybrid evolutionary algorithm allows addressing the problem of the energy landscape deceptiveness, allowing to obtain a set of optimized and diverse protein folds.

# Dependencies

* imageio==2.9.0
* matplotlib==3.3.3
* numpy==1.19.4
* pandas==1.1.4
* seaborn==0.11.0
* rosetta==0.3
* mpi4py==3.0.3

# Installation

This package is only compatible with Python 3.4 and above. To install this package, please follow the instructions below:

* Install the previous descripted dependencies
* Download and install PyRosetta following the instructions found at http://www.pyrosetta.org/dow
* Install the package itself:

```console
git clone https://github.com/danielvarela/PSPHybridDE.git
cd PSPHybridDE
pip install -r requirements.txt
```

# Basic Usage

1. Create a configuration file following the example found at https://github.com/danielvarela/PSPHybridDE/blob/main/configs/sample_psp.ini
```dosini
[inputs]
# pdb, secondary structure file and fragment files
pose_input=./input_files/info_1wit/vf_1wit.pdb
ss2_input=./input_files/info_1wit/vf_1wit_.psipred_ss2
frag3_input=./input_files/info_1wit/boinc_vf_aa1wit_03_05.200_v1_3
frag9_input=./input_files/info_1wit/boinc_vf_aa1wit_09_05.200_v1_3

[outputs]
# output folder and evolution file log
output_file=/resultados/evolution_sample.log

[DE]
# evolution algorithm selection strategy [Greedy, Crowding] 
selection=Greedy
# evolution algorithm parent strategy [RANDOM, BEST] 
scheme=RANDOM 
# population size
popsize=4
# mutation rate (weight factor F) 
mutate=0.001
# crossover probability (CR) 
recombination=0.1
# maximum number of generations/iterations (stopping criteria)
maxiter=12
# hybrid local search strategy [None, stage1, stage2, stage3, stage4]
local_search=stage4

```
information about the DE parameters can be found at https://en.wikipedia.org/wiki/Differential_evolution

2. Run with the algorithm with the desired configuration

```console
python evodock.py configs/sample_psp.ini
```

run with mpi4py

```console
mpirun -np 2 python evodock.py configs/sample_psp.ini
```



# Protein Structure Prediction

One of the most important problems in molecular biology is to obtain the native structure of a protein from its primary structure, i.e., the amino acids chain. Ab-initio methods adopt different approaches for the protein structure representation. For example, coarse-grained protein representation models considered the phi, psi and omega angles of the backbone structure while the sidechains are representated by a centroid. In the ab initio protein structure prediction problem (PSP) many authors have been working on the use of search methods, specially evolutionary algorithms, employing the coarse-grained representation model provided by the Rosetta software suit.

# Differential Evolution Algorithm

Differential Evolution [Price97] is a population-based search method. DE creates new candidate solutions by combining existing ones according to a simple formula of vector crossover and mutation, and then keeping whichever candidate solution has the best score or fitness on the optimization problem at hand.


# Bibliography

* Varela, D., Santos, J. Protein structure prediction in an atomic model with differential evolution integrated with the crowding niching method. Nat Comput (2020). https://doi.org/10.1007/s11047-020-09801-7

* Storn, R., Price, K. Differential Evolution – A Simple and Efficient Heuristic for global Optimization over Continuous Spaces. Journal of Global Optimization 11, 341–359 (1997). https://doi.org/10.1023/A:1008202821328 

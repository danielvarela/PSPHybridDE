from pyrosetta import *
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.toolbox import cleanATOM
from rosetta.core.pack.task import TaskFactory, operation

init()

pose = pose_from_pdb(sys.argv[-1])
tf = TaskFactory()
tf.push_back(operation.RestrictToRepacking())
extrarot = operation.ExtraRotamersGeneric()
extrarot.ex1(True)
extrarot.ex2aro(True)
tf.push_back(extrarot)
scorefxn = get_fa_scorefxn()
packer = PackRotamersMover(scorefxn)
packer.task_factory(tf)
packer.apply(pose)
pose.dump_scored_pdb("packed.pdb", scorefxn)

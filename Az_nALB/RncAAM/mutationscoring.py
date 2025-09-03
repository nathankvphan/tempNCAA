import optparse
from rosetta import *
from pyrosetta import *
from rosetta.core.pack.task import TaskFactory
import numpy as np
import pandas as pd 
from rosetta.core.kinematics import MoveMap
from rosetta.protocols import minimization_packing as pack_min
from rosetta.protocols import relax as rel
from rosetta.protocols.antibody.residue_selector import CDRResidueSelector
from rosetta.protocols.antibody import *
from rosetta.protocols.loops import *
from rosetta.protocols.relax import FastRelax
from rosetta.core.scoring import ScoreFunctionFactory
from pyrosetta.rosetta.core.kinematics import Jump, Stub
from pyrosetta.rosetta.numeric import xyzVector_double_t
from rosetta.core.select import residue_selector as selections
from rosetta.core.pose import add_variant_type_to_pose_residue
from rosetta.core import select
from rosetta.core.select.movemap import *
from pyrosetta.rosetta.protocols.interface import select_interface_residues
param_path = "/users/PAS2252/nathankvphan4/Projects/ncPPI/numbatch"
path = '../Parameters'
init(f'-extra_res_path {path} -mute all -constant_seed');

def mutationfxn(pose, chain, pos, ncaa, radius, trial, mutation_relax = True, mutation_dump = True):
    rts = rosetta.core.chemical.ChemicalManager.get_instance().residue_type_set("fa_standard")
    
    base = Pose()
    base.assign(pose)
    ncaa_scorefxn = ScoreFunctionFactory.create_score_function("ref2015")
    wt_scorefxn = ScoreFunctionFactory.create_score_function("ref2015")
    wt = Pose()
    wt.assign(base)
    mutant = mutater(base, chain, ncaa, pos)
    if mutation_relax:
        mutant = chillaxer(mutant, chain, ncaa_scorefxn, pos, radius)
        #mutant = minimizer(mutant, chain, ncaa_scorefxn, pos, radius)
    else:
        pass
    if mutation_dump:
        mutant.dump_pdb(f"../outputfiles/peptide_pos{pos}_trial{trial}_{ncaa}.pdb")
    else:
        pass
    wt_score = binding_energy(wt, wt_scorefxn, chain)
    mutant_score = binding_energy(mutant, ncaa_scorefxn, chain)
    return(mutant_score - wt_score), wt, mutant

def Get_backbone_delta(wt,mutant):
    return pyrosetta.rosetta.core.scoring.CA_rmsd(wt, mutant)

def binding_energy(pose, scorefxn, moving_chain):
    test_pose = Pose()
    test_pose.assign(pose)
    before_score = scorefxn(test_pose)
    translation_vec=xyzVector_double_t(1000, 0, 0)
    moving_residues = [i for i in range(1, test_pose.total_residue() + 1)
                       if test_pose.pdb_info().chain(i) == moving_chain]
    if not moving_residues:
        raise ValueError(f"Chain '{moving_chain}' not found in pose.")
    for r in moving_residues:
        res = test_pose.residue(r)
        for a in range(1, res.natoms() + 1):
            current_xyz = res.xyz(a)
            res.set_xyz(a, current_xyz + translation_vec)
    after_score = scorefxn(test_pose)
    return before_score - after_score

def mutater(pose, chain, ncaa, pos):
    mutant = Pose()
    mutant.conformation().detect_disulfides()
    mutant.assign(pose)
    chm = core.chemical.ChemicalManager.get_instance()
    rts = chm.residue_type_set("fa_standard")
    ncaa_res = rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map( ncaa ) )
    pos = mutant.pdb_info().pdb2pose(chain, pos)
    mutant.replace_residue(pos, ncaa_res, True)
    return mutant

def g_chillaxer(pose, chain):
    not_chill = Pose()
    not_chill.assign(pose)
    mover = MoveMap()
    scorefxn = ScoreFunctionFactory.create_score_function("ref2015")
    nb_chain = chain 
    for i in range(1, pose.total_residue() + 1):
        if pose.pdb_info().chain(i) == nb_chain:
            mover.set_bb(i, True)
            mover.set_chi(i, True)
        else:
            mover.set_bb(i, False)
            mover.set_chi(i, False)
    chillaxer = FastRelax()
    chillaxer.set_movemap(mover)
    chillaxer.constrain_relax_to_start_coords(True)
    chillaxer.coord_constrain_sidechains(True)
    chillaxer.ramp_down_constraints(False)
    chillaxer.min_type('lbfgs_armijo_nonmonotone')
    chillaxer.set_scorefxn(scorefxn)
    chillaxer.apply(not_chill)
    return not_chill
def minimizer(pose, chain, scorefxn, pos, radius):
    not_chill = Pose()
    not_chill.assign(pose)
    mover = MoveMap()
    center = pose.residue(pos).nbr_atom_xyz()

    nb_chain = chain 

    for i in range(1, pose.total_residue() + 1):
        if pose.pdb_info().chain(i) == nb_chain and pose.residue(i).nbr_atom_xyz().distance_squared(center) < radius**2:
            mover.set_bb(i, True)
            mover.set_chi(i, True)
        else:
            mover.set_bb(i, False)
            mover.set_chi(i, False)

    minimizer = rosetta.protocols.minimization_packing.MinMover()
    minimizer.movemap(mover)
    minimizer.score_function(scorefxn)
    minimizer.min_type('lbfgs_armijo_nonmonotone')
    minimizer.tolerance(0.001)
    minimizer.apply(not_chill)

    return not_chill

def chillaxer(pose, chain, scorefxn, pos, radius):
    not_chill = Pose()
    not_chill.assign(pose)
    mover = MoveMap()
    center = pose.residue(pos).nbr_atom_xyz()
    nb_chain = chain 
    Interface = 'A_B'
    interface_residues = select_interface_residues(
        not_chill,  Interface , 10)
    outputs = []
    for i in range(1, pose.total_residue() + 1):
        if (pose.pdb_info().chain(i) == nb_chain and pose.residue(i).nbr_atom_xyz().distance_squared(center) < radius**2) or (pose.pdb_info().chain(i) == nb_chain and interface_residues[i]):
            mover.set_bb(i, True)
            mover.set_chi(i, True)
            if interface_residues[i] and pose.pdb_info().chain(i) == nb_chain:
                outputs.append(i)
        else:
            mover.set_bb(i, False)
            mover.set_chi(i, False)
    chillaxer = FastRelax()
    chillaxer.set_movemap(mover)
    chillaxer.constrain_relax_to_start_coords(True)
    chillaxer.coord_constrain_sidechains(True)
    chillaxer.ramp_down_constraints(False)
    chillaxer.min_type('lbfgs_armijo_nonmonotone')
    chillaxer.set_scorefxn(scorefxn)
    chillaxer.apply(not_chill)

    return not_chill

def scoring_chump(wt, chain, pos, radius, trials, positions, ncaa, ncaa_str, mutation_relax=True, mutation_dump=False):
    pose = Pose()
    pose.assign(wt)

    if pose.pdb_info().chain(1) != 'A':
        chains = pose.split_by_chain()
        # print(chains[1].residue(1))
        # print(chains[2].residue(1))
        new_pose = Pose()
        new_pose.assign(chains[2])  # chain B (indexing starts from 0)
        new_pose.append_pose_by_jump(chains[1], new_pose.total_residue())

        wt = new_pose

    nanobody = Pose()
    pdb_info = wt.pdb_info()

    Mutation_residues = [i for i in range(1, wt.size() + 1) if wt.pdb_info().chain(i) == chain]
    Target_residues = [i for i in range(1, wt.size() + 1) if wt.pdb_info().chain(i) != chain]

    if not Mutation_residues:
        raise ValueError("Target mutation chain not found in pose.")
    
    Mutation_first_residue = Mutation_residues[0]
    Mutation_last_residue = Mutation_residues[-1]

    Target_first_residue = Target_residues[0]
    Target_last_residue = Target_residues[-1]

    scorefxn = ScoreFunctionFactory.create_score_function("ref2015")
    
    data_line = []
    for trial in trials:
        for pos in positions:
            ddg_score, wt, mutant = mutationfxn(wt, chain, pos, radius, trial, mutation_relax, mutation_dump)
            data_entry = [pos, 
                            str(trial), 
                            str(wt.sequence())[Mutation_first_residue-1:Mutation_last_residue], 
                            str(mutant.sequence())[Mutation_first_residue-1:Mutation_last_residue], 
                            str(ncaa), 
                            float(ddg_score), 
                            Get_backbone_delta(wt, mutant), 
                            float(scorefxn(mutant) - scorefxn(wt)), 
                            str(ncaa_str)]
            data_line.append(data_entry)
    df = pd.DataFrame(np.array(data_line))
    df.columns = ["Position", "Trial", "WT Sequence", "Mutant Sequence", "NCAA", "ddG", "Backbone Delta", "Score Difference", "NCAA String"]
    return df
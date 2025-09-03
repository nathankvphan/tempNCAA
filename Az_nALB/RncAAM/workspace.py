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
init()


def scanning(wt, chain, positions, ncaa, radius, output, trials):


    nanobody = Pose()
    pdb_info = pose.pdb_info()
    
    selected_residues = [i for i in range(1, pose.size() + 1) if pose.pdb_info().chain(i) == chain]
    
    if not selected_residues:
        raise ValueError(f"Chain {chain_id} not found in pose.")
    
    
    print('='*80)
    selected_residues = np.array(selected_residues)-1
    print(selected_residues)
    print(pose.residue(1))
    print('='*80)
    start = selected_residues[0]
    
    end = selected_residues[-1]
    

    mutation_pos = positions

    wt_sequence = str(wt.sequence())[start:end]
    
    
    print('='*80)
    print(wt_sequence[0])
    print('='*80)
    
    scorefxn = ScoreFunctionFactory.create_score_function("ref2015")

    final_array = []
    for trial in trials:
        for pos in mutation_pos:
            #if 0 in mutation_pos:
                #mutated_sequence = str([ncaa]) + wt_sequence[1:]
            #elif len(wt_sequence) in mutation_pos:
                #mutated_sequence = wt_sequence[:-1] + str([ncaa])
            #elif 2 in mutation_pos:
                #mutated_sequence = wt_sequence[:1] + str([ncaa]) + wt_sequence[2:]
            #else:
                #mutated_sequence = wt_sequence[:pos-1] + str([ncaa]) + wt_sequence[pos:]
            score_num, wt, mutant = score(wt, chain, pos, radius, trial)
            entry_i = [pos, str(trial), str(wt_sequence), str(mutant.sequence())[start:end], str(ncaa),
                       score_num, Get_backbone_delta(wt, mutant), float(scorefxn(mutant)-scorefxn(wt)) ]
            final_array.append(entry_i)
    df = pd.DataFrame(np.array(final_array))
    sorted = df.sort_values(by=[0, 1], ascending=True)
    sorted.columns = ['Position', 'Trial', 'WT Sequence', 'Mutated Sequence', 'Mutant', 'Energy Score', 'Backbone RMSD', 'Folded Energy Delta']
    return sorted
           
def score(pose, chain, pos, radius, trial):
    base = Pose()
    base.assign(pose)
    ncaa_scorefxn = ScoreFunctionFactory.create_score_function("ref2015")
    wt_scorefxn = ScoreFunctionFactory.create_score_function("ref2015")
    
    wt = Pose()
    wt.assign(base)

    mutant = mutater(base, chain, ncaa, pos)
    mutant = chillaxer(mutant, chain, ncaa_scorefxn, pos, radius)
    #mutant = g_chillaxer(mutant, chain)
    mutant.dump_pdb(f"/users/PAS2252/nathankvphan4/ncPPI/OutputFiles/pos{pos}_trial{trial}_{ncaa}.pdb")

    wt_score = binding_energy(wt, wt_scorefxn)
    mutant_score = binding_energy(mutant, ncaa_scorefxn)

    return(mutant_score - wt_score), wt, mutant

def Get_backbone_delta(wt,mutant):
    return pyrosetta.rosetta.core.scoring.CA_rmsd(wt, mutant)

def binding_energy(pose, scorefxn):
    test_pose = Pose()
    test_pose.assign(pose)

    before = scorefxn(test_pose)

    xyz = rosetta.numeric.xyzVector_double_t(1000, 0,0)

    chain2starts = len(pose.chain_sequence(1)) + 1
    for r in range(chain2starts, test_pose.total_residue() + 1):
        for a in range(1, test_pose.residue(r).natoms() + 1):
            test_pose.residue(r).set_xyz(a,
                test_pose.residue(r).xyz(a) + xyz)
    after = scorefxn(test_pose)
    # test_pose.dump_pdb(f"/users/PAS2252/nathankvphan4/ncPPI/OutputFiles/binding_energy_{pose.pdb_info().name()}.pdb")

    return before - after

def mutater(pose, chain, ncaa, pos):
    mutant = Pose()
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
    not_chill.dump_pdb("example.pdb")

    return not_chill
    




#=============================================================================================================================================

parser = optparse.OptionParser()
parser.add_option('-i', dest = 'input')
parser.add_option('-n', dest = 'nanobody_chain')
parser.add_option('-p', dest = 'positions')
parser.add_option('-t', dest = 'trials')
parser.add_option('-r', dest = 'radius')
parser.add_option('-o', dest = 'output')
(options, args) = parser.parse_args()
input = options.input
chain = options.nanobody_chain
positions = list(map(int,options.positions.split(',')))
trials = range(int(options.trials))
radius = int(options.radius)
output = options.output

pose = Pose()
pose_from_file(pose, input)
print(pose.residue(1))

if pose.pdb_info().chain(1) == 'A':
    pose.conformation().detect_disulfides()
    pose = g_chillaxer(pose, chain)
else:
    chains = pose.split_by_chain()
    # print(chains[1].residue(1))
    # print(chains[2].residue(1))
    new_pose = Pose()
    new_pose.assign(chains[2])  # chain B (indexing starts from 0)
    new_pose.append_pose_by_jump(chains[1], new_pose.total_residue())
    print(new_pose.residue(1))

    pose = new_pose
    pose = g_chillaxer(new_pose, chain)


import time


ncaadf = pd.read_csv('residue_names.csv')
# ncaa_array = ncaadf['Residue Name'].to_numpy()
ncaa_array = ['B30', 'A94'] # for C1=CC=C(C=C1)C(=O)C2=CC=C(C=C2)C[C@@H](C(=O)O)N

ncaa_array = ['DEM']
#ncaa_array = ['BPA']
outdf_list = []
errors = []
start_time = time.time()
for ncaa in reversed(ncaa_array):
    if ncaa in ['C00','B95', 'ZOT', 'Z7X', 'Z6K', 'YJM', 'YEU', 'Y1Z', 'XO7', 'XCU', 'WIN', 'WGR', 'V8M', 'V4Y', 'URF', 'UKS', 'TVN', 'TOL', 'TC9', 'SHY', 'SHN', 'SGA', 'SFE', 'RVH', 'RLQ', 'QF1', 'PLS', 'PCK', 'P4H', 'NFC', 'N16', 'MNJ', 'M60', 'M4N', 'LFA', 'L2Y', 'L2N', 'KW8', 'J1S', 'I70', 'H0B', 'GV1', 'GAY', 'FQC', 'FEQ', 'F0Z', 'EZ0', 'EYJ', 'D09', 'CU9', 'C4A', '73Y', '31B', '9GN', '9EF', '7EI', '1SU']:
        continue
    # try:
    finaldf = scanning(pose, chain, positions, ncaa,
    radius, output, trials)
    outdf_list.append(finaldf)
    print('='*90)
    print(outdf_list)
    outdf = pd.concat(outdf_list, ignore_index=True) 
    end_time = time.time()
    elapsed_time = end_time - start_time
    # except Exception as e:
    #     print(f"[ERROR] ncaa: {ncaa} — {e}")
    #     errors.append(str(ncaa))
if outdf_list:
    outdf = pd.concat(outdf_list, ignore_index=True)
else:
    print("[WARNING] No successful outputs to concatenate.")
    outdf = pd.DataFrame()  # or handle accordingly
print('='*90)
print(outdf_list)
print("ncaa length:" + str(len(ncaa_array)))
outdf = pd.concat(outdf_list, ignore_index=True) 
outdf.to_csv('nbonlyextendedb1.csv', index=False)
print('Completed')
print('Elapsed time: {:.2f} seconds'.format(elapsed_time))
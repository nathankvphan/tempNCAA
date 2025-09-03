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
from RncAAM.mutationscoring import mutationfxn, Get_backbone_delta, binding_energy, chillaxer, mutater

import os
path = '../Parameters'
init(f'-extra_res_path {path} -mute all -constant_seed');
class MutationScorer():
    def __init__(self, wt, chain, positions, ncaa, ncaa_str, radius, 
                 trials, mutationrelx = True, mutation_dump = True, mutationfxn=mutationfxn):
        self.wt = wt
        self.chain = chain
        self.positions = positions
        self.ncaa = ncaa
        self.radius = radius
        self.trials = trials
        self.ncaa_str = ncaa_str
        self.mutationfxn = mutationfxn
        self.Get_backbone_delta = Get_backbone_delta
        self.mutationrelx = mutationrelx
        self.mutation_dump = mutation_dump


    def score(self, wt, chain, pos, ncaa, radius):
        pose = Pose()
        pose.assign(wt)
        pose.conformation().detect_disulfides()
        # print(chain)
        # print(pose.pdb_info().chain(1))

        if pose.pdb_info().chain(1) != chain:
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
        for trial in self.trials:
            for pos in self.positions:
                ddg_score, wt, mutant = self.mutationfxn(wt, chain, pos, ncaa, radius, trial, self.mutationrelx, self.mutation_dump)
                data_entry = [pos, 
                              str(trial), 
                              str(wt.sequence())[Mutation_first_residue-1:Mutation_last_residue], 
                              str(mutant.sequence())[Mutation_first_residue-1:Mutation_last_residue], 
                              str(wt.sequence())[Target_first_residue-1:Target_last_residue],
                              str(self.ncaa), 
                              float(ddg_score), 
                              self.Get_backbone_delta(wt, mutant), 
                              float(scorefxn(mutant) - scorefxn(wt)), 
                              str(self.ncaa_str)]
                data_line.append(data_entry)
        df = pd.DataFrame(np.array(data_line))
        df.columns = ["Position", "Trial", "WT Sequence", "Mutant Sequence", "Target Sequence","NCAA", "ddG", "Backbone Delta", "Score Difference", "NCAA String"]
        return df


def get_scanner(wt, chain, positions, ncaa, ncaa_str, radius, 
                 trials, mutationrelx = True, mutation_dump = True, mutationfxn=mutationfxn):
    return MutationScorer(
        wt,
        chain,
        positions,
        ncaa,
        ncaa_str,
        radius,
        trials,
        mutationrelx,
        mutation_dump
    )
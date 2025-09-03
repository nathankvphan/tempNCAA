from pyrosetta import *
from rosetta import *
import pandas as pd
from RncAAM.mutationscoring import g_chillaxer
from RncAAM.scanning import get_scanner

import os

path = '../Parameters'
#Added constant seed. remove for production
init(f'-extra_res_path {path} -mute all -constant_seed');

def filter_by_std(group):
        mean = group["ddG"].mean()
        std = group["ddG"].std()
        return group[(group["ddG"] >= mean - std) & (group["ddG"] <= mean + std)]

def filter_by_std_fold(group):
        mean = group["Score Difference"].mean()
        std = group["Score Difference"].std()
        return group[(group["Score Difference"] >= mean - std) & (group["Score Difference"] <= mean + std)]

def group(df):
    df = df.sort_values(by=["Position"])
    df["ddG"] = pd.to_numeric(df["ddG"])
    filtered_df = df.groupby(['Position', 'Mutant Sequence', 'NCAA String'], group_keys=False).apply(filter_by_std)
    
    # Optional: Return the mean or std of the filtered energy scores per group
    grouped = filtered_df.groupby(['Position', 'Mutant Sequence', 'NCAA String'])["ddG"].mean().reset_index()

    grouped_df = grouped.sort_values(by=["ddG"], ascending=True)
    score = grouped_df['ddG'][0]
    
    return float(score)

def folding_group_score(df):
    df = df.sort_values(by=["Position"])
    df["Score Difference"] = pd.to_numeric(df["Score Difference"])
    # print(df['Score Difference'])
    filtered_df = df.groupby(['Position', 'Mutant Sequence', 'NCAA String'], group_keys=False).apply(filter_by_std_fold)
    
    grouped = filtered_df.groupby(['Position', 'Mutant Sequence', 'NCAA String'])["Score Difference"].mean().reset_index()

    grouped_df = grouped.sort_values(by=["Score Difference"], ascending=True)
    score = df['Score Difference'].mean()
    
    return float(score)

def Mutation_Wrapper(complexed_pose, chain, positions, ncaa, smiles, radius_relax = 7, trials = 3, mutationrelx = False, mutation_dump = False):
    trials = range(trials)

    # complex_p = g_chillaxer(complex_p, chain)
    large_op = []
    for i in range(len(list(ncaa))):
        ncaa = ncaa
        smiles = smiles


        scanner = get_scanner(complexed_pose, chain, positions, ncaa, smiles, radius_relax, trials, mutationrelx = mutationrelx, mutation_dump = mutation_dump)
    
        df = scanner.score(complexed_pose, chain, positions, ncaa, radius_relax)

        large_op.append(df)

        fdf = pd.concat(large_op, ignore_index=True)
    
    print(folding_group_score(fdf))
   
 

    return group(fdf), folding_group_score(fdf)
    
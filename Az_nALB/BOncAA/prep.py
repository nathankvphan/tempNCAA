import torch
import numpy as np
from rosetta import *
from pyrosetta import *
from pyrosetta.rosetta.protocols.interface import select_interface_residues
from pyrosetta.rosetta.core.scoring.sasa import rel_per_res_sc_sasa

def coordinates(target):
    coordinates = []
    for i in range(1, target.size()+1):
        residue = target.residue(i)
        if residue.has("CA"):
            ca_res = residue.xyz("CA")
            coordinates.append([ca_res.x, ca_res.y, ca_res.z])
    return coordinates

def binding_coordinates(target, binding_region):
    b_coordinates = []
    non_binding_coordiantes = []
    for i in range(1,target.size()+1):
        residue = target.residue(i)
        if residue.has("CA") and i in binding_region:
            ca_res = residue.xyz("CA")
            b_coordinates.append([ca_res.x, ca_res.y, ca_res.z])
        else:
            ca_res = residue.xyz("CA")
            non_binding_coordiantes.append([ca_res.x, ca_res.y, ca_res.z])
    return np.array(b_coordinates), np.array(non_binding_coordiantes)

def get_coordiantes(target, region):
    coordinates = []
    if type(region) == type(torch.tensor(region)):
        region = region.numpy()

    for i in range(1,target.size()+1):
        residue = target.residue(i)
        if residue.has("CA") and i in region:
            ca_res = residue.xyz("CA")
            coordinates.append([ca_res.x, ca_res.y, ca_res.z])
    return np.array(coordinates)

def find_binding_region(pose, target_chain_letter = 'B'):
    Interface = 'A_B'
    interface_residues = select_interface_residues(
        pose,  Interface , 7)
    target_interface = []
    non_interface = []
    for i in range(1, pose.size()+1):
        if interface_residues[i] and pose.pdb_info().chain(i) == target_chain_letter:
            target_interface.append(i)
        elif pose.pdb_info().chain(i) == target_chain_letter:
            non_interface.append(i)


    if target_interface[0] != 1:
        for i in range(1,pose.size()+1):
            if pose.pdb_info().chain(i) == target_chain_letter:
                firstval = i-1
                target_interface = np.array(target_interface) - np.ones(len(target_interface))* firstval
                non_interface = np.array(non_interface) - np.ones(len(non_interface)) * firstval
                break
    return np.array(target_interface, dtype=int), np.array(non_interface, dtype=int)

def binding_coordinates(target, binding_region):
    binding_coordinates = []
    non_binding_coordiantes = []
    for i in range(1,target.size()+1):
        residue = target.residue(i)
        if residue.has("CA") and i in binding_region and target.pdb_info().chain(i):
            ca_res = residue.xyz("CA")
            binding_coordinates.append([ca_res.x, ca_res.y, ca_res.z])
        elif residue.has("CA"):
            ca_res = residue.xyz("CA")
            non_binding_coordiantes.append([ca_res.x, ca_res.y, ca_res.z])
    return np.array(binding_coordinates), np.array(non_binding_coordiantes)

def distance_to_region(region_A, region_B):
    #Region A is this distance away from region B
    region_A = torch.from_numpy(region_A)
    region_B = torch.from_numpy(region_B)
    distances = torch.cdist(region_A, region_B, p = 2)
    return distances.mean(dim = 1)

def get_surface_exposures(target, target_range):
    surface_exposure_vector = np.array(rel_per_res_sc_sasa(target))
    cystine_free_vec = []
    for i in range(1, len(surface_exposure_vector)+1):
        if i in target_range:
            cystine_free_vec.append(surface_exposure_vector[i-1])
        else: 
            continue
    return torch.tensor(cystine_free_vec)
            


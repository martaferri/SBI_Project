#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 17 15:46:30 2018

@author: Aida
"""
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import Bio.PDB
import numpy
import os
import collections as col
import string
import utilities

def get_atoms_list(chain):
    atoms = chain.get_atoms()
    atoms_list = []
    for atom in atoms:
        if atom.get_name() == 'CA':
            atoms_list.append(atom)
    return atoms_list


inputs_dir = os.getcwd()+"/inputs/"
inputs_files = os.listdir(inputs_dir)

PDBparser = Bio.PDB.PDBParser()
super_imposer = Bio.PDB.Superimposer()

structures_list = []
chains_dic = {} #key: chain id / value: chain object
allchains = {'A': 'A', 'B': 'B', 'C': 'A', 'D': 'B'} #IMPORT

structureA = PDBparser.get_structure('hemo_chainA', 'chainA.pdb')
structureD = PDBparser.get_structure('hemo_chainD', 'chainD.pdb')
structureE = PDBparser.get_structure('nucleosoma_chainE_Mn', '3kuy_E.pdb')

chainA = structureA[0]['A']
chainD = structureD[0]['D']
chainE = structureE[0]['E']



def three_to_one(three_res_list): # returns seq string (one format). Used in the get_seq_from_pdbchain function
    one_res_list = []

    for res in three_res_list:
        try:
            one = utilities.three_to_one[res]
            one_res_list.append(one)
        except KeyError:
            return False
    return "".join(one_res_list)

def get_seq_from_pdbchain(chain): # returns the pdb seq in one letter format to make the alignment. Used in refine_for_superimpose
    three_res_list = []
    for res in chain:
        residues_atoms = res.get_atoms()
        for atom in residues_atoms:
            if atom.get_name() == 'CA':
                residue = atom.get_parent()
                three_res_list.append(residue.get_resname())
    return three_to_one(three_res_list) #three_to_one function


def one_to_three(chain_refined): # returns seq list of a chain (three format) from a seq in one letter format. Used in refine_for_superimpose
    three_res_list = []

    for res in chain_refined:
            three = utilities.one_to_three[res]
            three_res_list.append(three)
    return three_res_list


def refine_for_superimpose(fixedchain, movingchain): # returns a pattern of 0 and 1 to include or exclude the residue of the original chain
    fixedchain_seq = get_seq_from_pdbchain(fixedchain)
    movingchain_seq = get_seq_from_pdbchain(movingchain)


    alignments = pairwise2.align.globalxx(fixedchain_seq, movingchain_seq)

    fixedchain_seq = alignments[0][0]
    movingchain_seq = alignments[0][1]

    fixedchain_pattern = []
    movingchain_pattern = []

    for res1, res2 in zip(fixedchain_seq, movingchain_seq):
        if res1 == res2:
            fixedchain_pattern.append(1)
            movingchain_pattern.append(1)
        elif res1 == '-':
            movingchain_pattern.append(0)
        elif res2 == '-':
            fixedchain_pattern.append(0)

    chains_pattern = (fixedchain_pattern, movingchain_pattern)
    return chains_pattern


# def refine_for_superimpose(fixedchain, movingchain):  # returns the seq of the fixed and moving chain refined in one letter format
#     chains_refined = tuple()
#     fixedchain_seq = get_seq_from_pdbchain(fixedchain)
#     movingchain_seq = get_seq_from_pdbchain(movingchain)
#
#     alignments = pairwise2.align.globalxx(fixedchain_seq, movingchain_seq)
#
#     fixedchain_seq = alignments[0][0]
#     movingchain_seq = alignments[0][1]
#
#     fixedchain_refined = ""
#     movingchain_refined = ""
#
#     for res1, res2 in zip(fixedchain_seq, movingchain_seq):
#         if res1 == res2:
#             fixedchain_refined += res1
#             movingchain_refined += res2
#
#     chains_refined = (one_to_three(fixedchain_refined), one_to_three(movingchain_refined))
#     return chains_refined

# chains_refined = refine_for_superimpose(chainA, chainD)
#
# fixed_refined = chains_refined[0]
# moving_refined = chains_refined[1]

chains_pattern = refine_for_superimpose(chainA, chainD)

fixed_pattern = chains_pattern[0]
moving_pattern = chains_pattern[1]


# def get_chain_refined(chain_original, chain_refined_seq): # creates new chain objects filtering the residues of the original chain, to get their atoms later, with get_atoms_list() and superimpose
#     new_chain = Bio.PDB.Chain.Chain('X')
#
#     for residue in chain_original.get_residues():
#         for residue2 in chain_refined_seq:
#             if residue.get_resname() == residue2:
#                 new_chain.add(residue.copy())
#                 break
#     return new_chain

def get_chain_refined(chain_original, chain_pattern): # creates new chain objects filtering the residues of the original chain, to get their atoms later, with get_atoms_list() and superimpose
    new_chain = Bio.PDB.Chain.Chain('X')

    for residue, pattern in zip(chain_original.get_residues(), chain_pattern):
            if pattern == 1:
                new_chain.add(residue.copy())
    return new_chain


def get_atoms_list(chain):
    atoms = chain.get_atoms()
    atoms_list = []
    for atom in atoms:
        if atom.get_name() == 'CA':
            atoms_list.append(atom)
    return atoms_list

newchainA_result = get_chain_refined(chainA, fixed_pattern)
newchainD_result = get_chain_refined(chainD, moving_pattern)

newchainA_atoms = get_atoms_list(newchainA_result)
newchainD_atoms = get_atoms_list(newchainD_result)





print("END")

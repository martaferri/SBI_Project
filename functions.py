from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import Bio.PDB
import numpy
import os
import collections as col
import string
import utilities
import copy
import itertools

def get_distance(list1, list2):
    """Given two lists of coordinates, it calculates the distance."""
    import math
    distance_result = (list1[0] - list2[0]) ** 2 + (list1[1] - list2[1]) ** 2 + (list1[2] - list2[2]) ** 2
    return math.sqrt(abs(distance_result))

def check_type(chain):
    """Function to check if the chain is a prot or a nucleic acid"""
    atoms = chain.get_atoms()
    type_chain = ""
    list_c = []
    for element in atoms:
        list_c.append(element.get_name())
    if "CB" in list_c:
        type_chain = "protein"
    else:
        type_chain = "nucleic_acid"
    return type_chain

def adapt_chain(chain):
    """Checks the type of the chain and changes the C1' to CB if it is nucleic acid"""
    type_chain = check_type(chain)
    name = chain.id
    if type_chain == "nucleic_acid":
        new_chain = Bio.PDB.Chain.Chain(name)
        chain = copy.copy(chain)
        for residue in chain:
            new_chain.add(residue.copy())

        for residue in new_chain:
            for atom in residue:
                if atom.id == "C1'":
                    atom.id = "CB"
                    residue.add(atom.copy())
        return new_chain
    else:
        return chain

def get_atoms_list(chain):
    type_chain = check_type(chain)
    if type_chain == "protein":
        atom_id = "CA"
    elif type_chain == "nucleic_acid":
        atom_id = "P"
    atoms = chain.get_atoms()
    atoms_list = []
    for atom in atoms:
        if atom.id == atom_id:
            atoms_list.append(atom)
    return atoms_list

def three_to_one(three_res_list):
    """Returns seq string (one format). Used in the get_seq_from_pdbchain function"""
    one_res_list = []

    for res in three_res_list:
        try:
            one = utilities.three_to_one[res]
            one_res_list.append(one)
        except KeyError:
            return False
    return "".join(one_res_list)

def get_seq_from_pdbchain(chain):
    """Returns the pdb seq (in one letter format if it is protein) to make the alignment. Used in refine_for_superimpose"""
    type_chain = check_type(chain)
    if type_chain == "protein":
        three_res_list = []
        for res in chain:
            residues_atoms = res.get_atoms()
            for atom in residues_atoms:
                if atom.get_name() == 'CA':
                    residue = atom.get_parent()
                    three_res_list.append(residue.get_resname())
        return three_to_one(three_res_list)  # three_to_one function
    else:
        nucleic_acid_res = []
        for res in chain:
            residues_atoms = res.get_atoms()
            for atom in residues_atoms:
                if atom.get_name() == 'P':
                    residue = atom.get_parent()
                    nucleic_acid_res.append(residue.get_resname())
        nucleic_acid_seq = [x[2] for x in nucleic_acid_res]
        return "".join(nucleic_acid_seq)


#
# def get_seq_from_pdbchain(chain):
#     """Returns the pdb seq in one letter format to make the alignment. Used in refine_for_superimpose"""
#     three_res_list = []
#     for res in chain:
#         residues_atoms = res.get_atoms()
#         for atom in residues_atoms:
#             if atom.get_name() == 'CA':
#                 residue = atom.get_parent()
#                 three_res_list.append(residue.get_resname())
#
#     return three_to_one(three_res_list)  # three_to_one function

def one_to_three(chain_refined):
    """Returns seq list of a chain (three format) from a seq in one letter format. Used in refine_for_superimpose"""
    three_res_list = []

    for res in chain_refined:
        three = utilities.one_to_three[res]
        three_res_list.append(three)
    return three_res_list

def refine_for_superimpose(fixedchain, movingchain):
    """Returns a pattern of 0 and 1 to include or exclude the residue of the original chain"""
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

def get_chain_refined(chain_original, chain_pattern):
    """Creates new chain objects filtering the residues of the original chain,
    to get their atoms later, with get_atoms_list() and superimpose"""
    new_chain = Bio.PDB.Chain.Chain('X')

    for residue, pattern in zip(chain_original.get_residues(), chain_pattern):
        if pattern == 1:
            new_chain.add(residue.copy())
    return new_chain

def check_clashes(current_model, newchain):
    neighbors2 = set()
    atom_list = Bio.PDB.Selection.unfold_entities(current_model, 'A')
    ns = Bio.PDB.NeighborSearch(atom_list)
    # chain_list = Bio.PDB.Selection.unfold_entities(current_model, 'C')

    for at in newchain.get_atoms():
        center = at.get_coord()
        neighbors = ns.search(center, 1.5, level='C')
        for element in neighbors:
            if element != newchain:
                neighbors2.add(element)

    if len(neighbors2) != 0:
        return True
    else:
        return False

def unique_chains_fasta (number_to_letter):
    fo = open("unique_chains_fasta.mfa", "w")
    for key, value in number_to_letter.items():
        name_chain = "chain_" + value
        fo.write(">%s\n%s\n" % (name_chain, get_seq_from_pdbchain(key)))
    fo.close()

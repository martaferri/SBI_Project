import Bio.PDB
import numpy
import os
import collections as col
import utilities
import copy
from functions import *

def get_distance(list1, list2):
    """Given two lists of coordinates, it calculates the distance."""
    import math

    distance_result = (list1[0] - list2[0]) ** 2 + (list1[1] - list2[1]) ** 2 + (list1[2] - list2[2]) ** 2
    return math.sqrt(distance_result)


inputs_dir = os.getcwd() + "/interactions_results/"
inputs_files = os.listdir(inputs_dir)

PDBparser = Bio.PDB.PDBParser()

list_of_pairs = []
for pair in inputs_files:
    input_pair_list = []
    pdb_code = os.path.basename(os.path.splitext(pair)[0])
    pdb_filename = inputs_dir + pair
    structure = PDBparser.get_structure(pdb_code, pdb_filename)
    for model in structure:
        for chain in model:
            input_pair_list.append(chain)
    list_of_pairs.append(input_pair_list)

cutoff = 0.9
c = 0
redundancy_list = []
for pair0 in list_of_pairs:
    c += 1
    for pair1 in list_of_pairs[c:]:
        print("COMPARING %s AND %s" % (pair0, pair1))
        for chain in pair1:
            pair00_seq = get_seq_from_pdbchain(pair0[0])
            chain_seq = get_seq_from_pdbchain(chain)
            alignments = pairwise2.align.globalxx(pair00_seq, chain_seq)

            if len(alignments) > 0:
                score_chain = alignments[0][2]
                len_max = max(len(pair00_seq), len(chain_seq))
                cutoff_chain = score_chain / len_max
                if cutoff_chain >= cutoff:
                    fixedchain = pair0[0]
                    movingchain = chain
                    pair1_to_pop = copy.copy(pair1)
                    pair1_to_pop.pop(pair1.index(movingchain))
                    altchain = pair1_to_pop[0]
                    comparingchain = pair0[1]

                    fixed_atoms_list = get_atoms_list(fixedchain)
                    moving_atoms_list = get_atoms_list(movingchain)
                    alt_atoms_list = get_atoms_list(altchain)

                    new_chain = Bio.PDB.Chain.Chain("X")
                    altchain = copy.copy(altchain)
                    for residue in altchain.get_residues():
                        new_chain.add(residue.copy())

                    if len(fixed_atoms_list) != len(moving_atoms_list):
                        chains_pattern = refine_for_superimpose(fixedchain, movingchain)
                        fixed_pattern = chains_pattern[0]
                        moving_pattern = chains_pattern[1]

                        fixedchain = get_chain_refined(fixedchain, fixed_pattern)
                        movingchain = get_chain_refined(movingchain, moving_pattern)

                        fixed_atoms_list = get_atoms_list(fixedchain)
                        moving_atoms_list = get_atoms_list(movingchain)

                    super_imposer = Bio.PDB.Superimposer()
                    super_imposer.set_atoms(fixed_atoms_list, moving_atoms_list)
                    print("\t RMSD : %.40f" % (numpy.abs(super_imposer.rms)))
                    super_imposer.apply(new_chain.get_atoms())

                    new_chain_atoms_list = get_atoms_list(new_chain)
                    comparing_chain_atoms_list = get_atoms_list(comparingchain)

                    if len(new_chain_atoms_list) != len(comparing_chain_atoms_list):
                        chains_pattern = refine_for_superimpose(new_chain, comparingchain)
                        new_pattern = chains_pattern[0]
                        comparing_pattern = chains_pattern[1]

                        new_chain = get_chain_refined(new_chain, new_pattern)
                        comparingchain = get_chain_refined(comparingchain, comparing_pattern)

                        new_chain_atoms_list = new_chain.get_atoms()
                        comparing_chain_atoms_list = comparingchain.get_atoms()

                    distances = []
                    info = []
                    distances_info = []
                    for new_atom, comparing_atom in zip(new_chain_atoms_list, comparing_chain_atoms_list):
                        coords_new = new_atom.get_coord()
                        coords_comparing = comparing_atom.get_coord()

                        # print("TEST COORDS: %s %s" % (coords_new, coords_comparing))
                        # print(list(coords_new.tolist()), list(coords_comparing.tolist()))

                        distance = get_distance(list(coords_new.tolist()), list(coords_comparing.tolist()))
                        distances.append(distance)
                        # print("new_chain coords: %s\ncomparingchain coords: %s\nDistance: %.4f\n" % (coords_new, coords_comparing, distance))


                    max_distance = max(distances)
                    if max_distance < 9:
                        pairs = [pair0, pair1]
                        print("The pairs %s and %s are REDUNDANT")
                        redundancy_list.append(pairs)
            else:
                print("No sequence similarity between %s and %s\n" % (pair0, pair1))
                continue

print("\n\nREDUNDANCY LIST: %s" % redundancy_list)

to_delete_pair = []
for pairs in redundancy_list:
    to_delete_pair.append(pairs[1])

# idx_list = []
# final_inputs = []
# for pair in list_of_pairs:
#     for pair2 in to_delete_pair:
#             if pair[0] == pair2[0] and pair[1] == pair2[1]:
#                 break
#             else:
#                 if pair not in final_inputs:
#                     final_inputs.append(pair)


final_inputs = copy.copy(list_of_pairs)
for pair in to_delete_pair:
    for pair2 in final_inputs:
            if pair[0] == pair2[0] and pair[1] == pair2[1]:
                idx = final_inputs.index(pair)
                del final_inputs[idx]


print("\n\nFINAL LIST: %s" % final_inputs)


print("END")
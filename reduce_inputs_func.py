# python modules
import Bio.PDB
import numpy
import os
import gzip
import re
import collections as col
import copy
import argparse
import tarfile
import ast
import sys

# our modules
from classes import *
import utilities

def save_new_pair(chain1, chain2, outputs_dir, number_file):
    """Saves the non-redundant chain pairs (.pdb) obtained from the reducing inputs process."""
    output = (outputs_dir + "/reduced_inputs/")
    if not os.path.exists(output):
        os.makedirs(output)
    name = chain1.id + chain2.id + str(number_file)
    newmodel = Bio.PDB.Model.Model("new")
    newmodel.add(chain1)
    newmodel.add(chain2)
    io = Bio.PDB.PDBIO()
    io.set_structure(newmodel)
    io.save(output + name + ".pdb")

def reduce_inputs_func(inputs_files, outputs_dir, options_verbose):
    """From a set of input pairs (.pdb) compares each chain pair with the rest. This comparision has two steps: sequence
    and structural. First of all, a pairwise alignment is performed to determine similar sequences (cut-off = 0.9).
    The score of the alignment is normalized by the length of the longest sequence. If the normalized score is higher
    than the stablished cut-off, the analysis proceeds to the second step. In this step, a superimposition is performed
    between the similar chains. These similar chains are part of two different interaction pairs, which will be refered
    as fixed and moving chains. We apply the rotran matrix to the couple of the moving chain, which will be refered as
    alternative/new chain. Finally, the distances between the CA of the comparing chain (couple of the fixed chain) and
    the new chain are computed. If the distance is lower than 9A, the interactions will be considered the same (redundant)."""

    import functions
    if options_verbose:
        sys.stderr.write("\nGetting the non-redundant interactions... Please wait.\n\n")
    PDBparser = Bio.PDB.PDBParser()

    # Creating a list of lists of the input pairs
    list_of_pairs = []
    for pair in inputs_files:
        input_pair_list = []
        pdb_code = pair.split("/")[-1].split(".")[0]
        pdb_filename = pair
        structure = PDBparser.get_structure(pdb_code, pdb_filename)
        for model in structure:
            for chain in model:
                input_pair_list.append(chain)
        list_of_pairs.append(input_pair_list)

    # Cut-off to determine similarity between chain sequences. If the chain cut-off is greater than the specified cut-off,
    # these chains will be considered equal for the superimposition.
    cutoff = 0.9
    c = 0
    redundancy_list = []
    for pair0 in list_of_pairs:
        c += 1
        for pair1 in list_of_pairs[c:]:

            for chain in pair1:
                pair00_seq = functions.get_seq_from_pdbchain(pair0[0])
                chain_seq = functions.get_seq_from_pdbchain(chain)
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

                        fixed_atoms_list = functions.get_atoms_list(fixedchain)
                        moving_atoms_list = functions.get_atoms_list(movingchain)

                        new_chain = Bio.PDB.Chain.Chain("X")
                        altchain = copy.copy(altchain)
                        for residue in altchain.get_residues():
                            new_chain.add(residue.copy())

                        if len(fixed_atoms_list) != len(moving_atoms_list):
                            chains_pattern = functions.refine_for_superimpose(fixedchain, movingchain)
                            fixed_pattern = chains_pattern[0]
                            moving_pattern = chains_pattern[1]

                            fixedchain = functions.get_chain_refined(fixedchain, fixed_pattern)
                            movingchain = functions.get_chain_refined(movingchain, moving_pattern)

                            fixed_atoms_list = functions.get_atoms_list(fixedchain)
                            moving_atoms_list = functions.get_atoms_list(movingchain)

                        super_imposer = Bio.PDB.Superimposer()
                        super_imposer.set_atoms(fixed_atoms_list, moving_atoms_list)
                        super_imposer.apply(new_chain.get_atoms())

                        new_chain_atoms_list = functions.get_atoms_list(new_chain)
                        comparing_chain_atoms_list = functions.get_atoms_list(comparingchain)

                        if len(new_chain_atoms_list) != len(comparing_chain_atoms_list):
                            chains_pattern = functions.refine_for_superimpose(new_chain, comparingchain)
                            new_pattern = chains_pattern[0]
                            comparing_pattern = chains_pattern[1]
                            new_chain = functions.get_chain_refined(new_chain, new_pattern)
                            comparingchain = functions.get_chain_refined(comparingchain, comparing_pattern)

                            new_chain_atoms_list = functions.get_atoms_list(new_chain)
                            comparing_chain_atoms_list = functions.get_atoms_list(comparingchain)

                        if len(new_chain_atoms_list) == 0 or len(comparing_chain_atoms_list) == 0:
                            break

                        distances = []
                        for new_atom, comparing_atom in zip(new_chain_atoms_list, comparing_chain_atoms_list):
                            coords_new = new_atom.get_coord()
                            coords_comparing = comparing_atom.get_coord()
                            distance = functions.get_distance(list(coords_new.tolist()), list(coords_comparing.tolist()))
                            distances.append(distance)

                        max_distance = max(distances)
                        if max_distance < 9:
                            pairs = [pair0, pair1]
                            redundancy_list.append(pairs)
                else:
                    continue

    to_delete_pair = []
    for pairs in redundancy_list:
        to_delete_pair.append(pairs[1])

    final_inputs = copy.copy(list_of_pairs)
    for pair in to_delete_pair:
        for pair2 in final_inputs:
            if pair[0] == pair2[0] and pair[1] == pair2[1]:
                idx = final_inputs.index(pair)
                del final_inputs[idx]

    number_file = 1
    for element in final_inputs:
        save_new_pair(element[0], element[1], outputs_dir, number_file)
        number_file += 1
    if options_verbose:
        sys.stderr.write("The total number of redundant interactions found is %d\n" % (len(redundancy_list)))
        sys.stderr.write("The final number of non-redundant interactions is %d\n" % (len(final_inputs)))

    return "OK"

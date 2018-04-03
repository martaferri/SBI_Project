#!/usr/bin/env python3
import sys
import Bio.PDB
import numpy
import itertools
import re
import os
import copy
import utilities

def save_pairs(outputs_dir, interactions_list, complex_struct):
    """Creates the pair files (.pdb) obtained from the get interactions step."""
    output = (outputs_dir + "/get_interactions_results/")
    if not os.path.exists(output):
        os.makedirs(output)

    for element in interactions_list:
        pdb = open(complex_struct, 'r')
        fo = open(output + element + ".pdb", "w")
        for line in pdb:
            if line.startswith('ATOM'):
                line = line.strip()
                search = re.search(r'[A-Z]*\s*\d{1,5}\s*[A-Z]{1,3}[0-9]?\'?\s*[A-Z]{1,3}\s([a-zA-Z1-9]).*', line)
                if search:
                    chain = search.group(1)
                    if chain == element[0]:
                        fo.write("%s\n" % (line))
                    elif chain == element[1]:
                        fo.write("%s\n" % (line))
        fo.close()
        pdb.close()
    return output

def calc_distances_residues(chain1, chain2):
    """Function to calculate distances between CA. Returns a tuple with: min, max, number of atoms and a counter of
    distances under 8A."""
    import functions
    chain1_redef = copy.copy(functions.adapt_chain(chain1))
    chain2_redef = copy.copy(functions.adapt_chain(chain2))
    at1 = chain1_redef.get_atoms()
    at2 = chain2_redef.get_atoms()
    list1 = []
    list2 = []
    arr = []

    for atom1 in at1:
        if atom1.id == 'CA':
            list1.append(atom1)
    for atom2 in at2:
        if atom2.id == 'CA':
            list2.append(atom2)

    for combination in itertools.product(list1, list2):
        coords1 = combination[0].get_coord()
        coords2 = combination[1].get_coord()
        dist = functions.get_distance(list(coords1.tolist()), list(coords2.tolist()))
        arr.append(dist)

    num_min = 0
    if len(arr) != 0:
        nparr = numpy.array(arr)
        for dist in nparr:
            if dist < 8:
                num_min += 1

        dist_tuple = (min(abs(nparr)), max(abs(nparr)), len(nparr), num_min)
        return dist_tuple

def get_interactions(complex_struct, outputs_dir):
    """Gets all possible interactions given a macrocomplex (.pdb). It computes the distance between chains, and pairs
    those that accomplish the following conditions: less than 8A between CA and implication of at least 8 CA in this
    interaction."""
    import functions

    # Creating the object from the complex
    PDBparser = Bio.PDB.PDBParser()
    structure = PDBparser.get_structure('prot', complex_struct)

    # Loop that compares all the chains and stores in a dict of dicts each computation of the distances between them
    main_dict = {}
    set_aachains = set()

    for model in structure:
        list_residues = model.get_residues()
        for res in list_residues:
            if res.get_resname() in utilities.res_list:
                set_aachains.add(res.get_parent())

    for chain1 in set_aachains:
        main_dict[chain1.get_id()] = {}
        for chain2 in set_aachains:
            if chain1.get_id() != chain2.get_id():
                if chain2.get_id() not in main_dict:
                    main_dict[chain1.get_id()][chain2.get_id()] = calc_distances_residues(chain1, chain2)

    # Setting a cut-off value to determine if there is or not an interaction (8A) at least between 8 CA.
    # Storing the interactions in a list.
    interactions_list = []
    for element in main_dict:
        for key, value in main_dict[element].items():
            if value[0] <= 8 and value[3] >= 8:
                interactions_list.append(element + key)

    # Saving the pairs files (.pdb)
    output = save_pairs(outputs_dir, interactions_list, complex_struct)

    return output

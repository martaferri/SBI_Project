#!/usr/bin/env python3

import Bio.PDB
import numpy
import itertools

# Opening the pdb without chains I and J because they were giving problems
complex_struct = open('3kuy_ok.pdb', 'r')

# Creating the object from the complex
p = Bio.PDB.PDBParser()
structure = p.get_structure('3kuy', complex_struct)

# Function to calculate distances between atoms, return a tuple with min, max and number of atoms
def calc_distances_residues(chain1, chain2):
	countm = 0
	countn = 0
	at1 = chain1.get_atoms()
	at2 = chain2.get_atoms()
	list1 = []
	list2 = []
	arr = []

	for atom1 in at1:
		if atom1.get_name() == 'CB':
			list1.append(atom1)
	for atom2 in at2: 
		if atom2.get_name() == 'CB':
			list2.append(atom2)
					
	for combination in itertools.product(list1, list2):
		dist = combination[0] - combination[1]
		arr.append(dist)

	nparr = numpy.array(arr)
	dist_tuple = (min(abs(nparr)), max(abs(nparr)), len(nparr))
	return dist_tuple


# Loop that compares all the chains and stores in a dict of dicts each computation of the distances between them
main_dict = {}
for model in structure:
	for chain1 in model:
		# print("1",chain1)
		main_dict[chain1.get_id()] = {}
		for chain2 in model:
			if chain1.get_id() != chain2.get_id():
				if chain2.get_id() not in main_dict:
					# print("2",chain2)
					main_dict[chain1.get_id()][chain2.get_id()] = calc_distances_residues(chain1, chain2)

print(main_dict)

# Setting a cut-off value to determine if there is or not an interaction. Storing the interactions in a list
interactions_list = []
for element in main_dict:
	for key, value in main_dict[element].items():
		if value[0] < 6: 
			interactions_list.append(element + key)
print("\n")
print(interactions_list)


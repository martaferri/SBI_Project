#!/usr/bin/env python3
import sys
import Bio.PDB
import numpy
import itertools
import re


if len(sys.argv) <= 1:
    print('USAGE: python3 file.pdb')
    exit()

# Opening the pdb without chains I and J because they were giving problems
complex_struct = open(sys.argv[1], 'r')

# Creating the object from the complex
p = Bio.PDB.PDBParser()
structure = p.get_structure('prot', complex_struct)

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

	if len(arr) != 0: 
		nparr = numpy.array(arr)
		dist_tuple = (min(abs(nparr)), max(abs(nparr)), len(nparr))
		return dist_tuple

# Loop that compares all the chains and stores in a dict of dicts each computation of the distances between them
main_dict = {}
aa_list = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]
set_aachains = set()

for model in structure: 
	list_residues = model.get_residues()
	for res in list_residues:
		if res.get_resname() in aa_list:
			set_aachains.add(res.get_parent())
			print (set_aachains)
	
	for chain1 in set_aachains:
		main_dict[chain1.get_id()] = {}
		for chain2 in set_aachains:
			if chain1.get_id() != chain2.get_id():
				if chain2.get_id() not in main_dict:
					main_dict[chain1.get_id()][chain2.get_id()] = calc_distances_residues(chain1, chain2)

# Setting a cut-off value to determine if there is or not an interaction. Storing the interactions in a list
interactions_list = []
for element in main_dict:
	for key, value in main_dict[element].items():
		if value[0] <= 8: 
			interactions_list.append(element + key)

print(interactions_list)

complex_struct.close()

# Creating the new pair files
for element in interactions_list:
	pdb = open(sys.argv[1], 'r')
	fo = open(element + ".pdb", "w")
	for line in pdb:
		if line.startswith('ATOM'):
			line = line.strip()
			search = re.search(r'[A-Z]*\s*\d{1,3}\s*[A-Z]{1,2}[0-9]?\s*[A-Z]{3}\s([a-zA-Z]).*', line)
			if search:
				chain = search.group(1)			
				if chain == element[0]:
					fo.write("%s\n"%(line))
				elif chain == element[1]:
					fo.write("%s\n"%(line))
	fo.close()
	pdb.close()	

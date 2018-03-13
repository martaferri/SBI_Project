#!/usr/bin/env python3

import Bio.PDB
import numpy
import itertools


chain1 = open('3kuy_A.pdb', 'r')
chain2 = open('3kuy_B.pdb', 'r') 

p1 = Bio.PDB.PDBParser()
struct1 = p1.get_structure('3kuy_A', chain1)
print (struct1)
p2 = Bio.PDB.PDBParser()
struct2 = p2.get_structure('3kuy_B', chain2)
print (struct2)


def calc_distances_residues(struct1, struct2):
	countm = 0
	countn = 0
	at1 = struct1.get_atoms()
	at2 = struct2.get_atoms()
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
	print (min(abs(nparr)))
	print (max(abs(nparr)))
	print (len(nparr))



calc_distances_residues(struct1, struct2)

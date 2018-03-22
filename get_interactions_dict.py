#!/usr/bin/env python3
import sys, os
import Bio.PDB
import numpy
import itertools
import re
import chain_dict as cd


if len(sys.argv) <= 1:
    print('USAGE: python3 script.py file.pdb file.fasta')
    exit()

################################
######## chain_dict.py  ########
################################

pdb_file = open(sys.argv[1], "r")

output = open("dict_and_interactions.txt", "w")


# Function to create a fasta file from the pdb file
cd.pdb_to_fasta(pdb_file)

pdb_file.close()


# Now we want to make pariwise comparison to see if there were chains in the pdb that corresponded to the same chain in the fasta
fasta = open(sys.argv[2])
fasta_from_pdb = open("fasta_from_pdb.fasta", "r")
alignments_results = open("alignments_results.txt", "w")

# Create two dict for both files to make the comparisons
dict1 = cd.FASTA_iterator(fasta_from_pdb)
dict2 = cd.FASTA_iterator(fasta)


# Function that makes pairwise comparisons an creates a resume file of the alignments and a dictionary
dict_ids = cd.dictionary_chains(dict1, dict2)

fasta_from_pdb.close()
alignments_results.close()

output.write("DICTIONARY OF IDs:\n\t%s\n" % dict_ids)

output.write("\n\n")

################################
##### get_interactions.py  #####
################################

# Opening the pdb
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

# print(main_dict)

# Setting a cut-off value to determine if there is or not an interaction. Storing the interactions in a list
interactions_list = []
for element in main_dict:
    for key, value in main_dict[element].items():
        if value[0] <= 8:
            interactions_list.append(element + key)

output.write("INTERACTIONS LIST:\n\t%s\n" % interactions_list)

complex_struct.close()

# Creating the new pair files
interactions = ("./interactions/")
if not os.path.exists(interactions):
    os.makedirs(interactions)

for element in interactions_list:
    pdb = open(sys.argv[1], 'r')
    # fo = open(sys.argv[1][0:-4] + "_" + element + ".pdb", "w")
    fo = open(interactions + element + ".pdb", "w")
    for line in pdb:
        if line.startswith('ATOM'):
            line = line.strip()
            search = re.search(r'[A-Z]*\s*\d{1,3}\s*[A-Z]{1,2}[0-9]?\s*[A-Z]{3}\s([A-Z]).*', line)
            if search:
                chain = search.group(1)
                if chain == element[0]:
                    fo.write("%s\n" % (line))
                elif chain == element[1]:
                    fo.write("%s\n" % (line))
    fo.close()
    pdb.close()

output.close()


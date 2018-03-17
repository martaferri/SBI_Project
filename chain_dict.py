import sys, re
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import utilities

# Definition of usage: It works with the pdb file of the pair and the fasta file of all the chains
if len(sys.argv) <= 1:
    print('USAGE: python3 pdb2fasta.py file.pdb file.fasta')
    exit()
    
pdb_file = open(sys.argv[1], "r")
fasta_from_pdb = open("fasta_from_pdb.fasta", "w")


# Function to create a fasta file from the pdb file
def pdb_to_fasta(pdb_file):
	main_dict = {}
	prev = "-1"
	for line in pdb_file:
	    toks = line.split()
	    if len(toks)<1: 
	        continue
	    if toks[0] != 'ATOM': 
	        continue 
	    if toks[5] != prev:        
	        aa = toks[3]

	        if toks[4] not in main_dict:
	            main_dict[toks[4]] = ""
	            if aa not in utilities.letters:
	                continue
	            else:
	                main_dict[toks[4]] += utilities.letters[aa]
	        else:
	            if aa not in utilities.letters:
	                continue
	            else:    
	                main_dict[toks[4]] += utilities.letters[aa]           
	    prev = toks[5]

	name = sys.argv[1][0:4]
	for key, value in main_dict.items():
	    name_chain = name + ":" + key 
	    fasta_from_pdb.write(">%s\n%s\n"%(name_chain, value))

pdb_to_fasta(pdb_file)

pdb_file.close()
fasta_from_pdb.close()

########################

# Now we want to make pariwise comparison to see if there were chains in the pdb that corresponded to the same chain in the fasta
fasta = open(sys.argv[2])
fasta_from_pdb = open("fasta_from_pdb.fasta", "r")
alignments_results = open("alignments_results.txt", "w")

# Create two dict for both files to make the comparisons
def FASTA_iterator(fasta_filename):
	fasta_dict = {}
	for line in fasta_filename:
		line = line.strip()
		if line.startswith(">"):
			name = line[1:]
			fasta_dict[name] = ""
		else:
			fasta_dict[name] += line
	return fasta_dict		

dict1 = FASTA_iterator(fasta_from_pdb)
dict2 = FASTA_iterator(fasta)

# Function that makes pairwise comparisons an creates a resume file of the alignments and a dictionary
def dictionary_chains(dict1, dict2):
	dict_pdb_to_fasta = {}
	for key, value in dict1.items():
		score = 0
		dict_pdb_to_fasta[key] = ""
		for key2, value2 in dict2.items():
			alignments_results.write("This is the alignment between: %s and %s: "%(key, key2[:6]))
			alignments = pairwise2.align.globalxx(value, value2)
			
			if len(alignments) > 0:
				alignments_results.write("\tThe max score is %d\n" %(alignments[0][2]))
				if alignments[0][2] > score:
					score = alignments[0][2]
					dict_pdb_to_fasta[key] = key2
			else:
				alignments_results.write("\tNO ALIGNMENT\n")

	alignments_results.write("\nThe dictionary is:\nPDB\tFASTA\n")		
	for key, value in dict_pdb_to_fasta.items():
		alignments_results.write("%s\t%s\n"%(key, value[:6]))
	return dict_pdb_to_fasta


dict_pdb_to_fasta = dictionary_chains(dict1, dict2)

fasta_from_pdb.close()
alignments_results.close()

###################

print(dict_pdb_to_fasta)

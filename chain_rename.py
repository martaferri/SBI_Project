import sys, re
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

# Definition of usage: It works with the pdb file of the pair and the fasta file of all the chains
if len(sys.argv) <= 1:
    print('USAGE: python3 pdb2fasta.py file.pdb file.fasta')
    exit()
    
pdb_file = open(sys.argv[1], "r")
fasta_from_pdb = open("fasta_from_pdb.fasta", "w")

# Dict to convert the pdb file into fasta
letters = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H',
           'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
           'TYR':'Y','VAL':'V'}


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
	            if aa not in letters:
	                continue
	            else:
	                main_dict[toks[4]] += letters[aa]
	        else:
	            if aa not in letters:
	                continue
	            else:    
	                main_dict[toks[4]] += letters[aa]           
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

# Finally, we substitute the chains of the original pdb file of the pair for the corresponding chain found in the fasta file 
pair_file = open(sys.argv[1], "r")
renamed_file = open(sys.argv[1][0:-4] + "_renamed.pdb", "w")

# Simplification of the dict to just the letter of the chain
dict_simple = {}
for key, value in dict_pdb_to_fasta.items():
    if value:
        dict_simple[key[-1:]] = value[5]

# Function that creates a new pdb file with the chains renamed
def chain_rename(pdb_file):
	for line in pair_file:
	    line = line.strip()
	    if line.startswith('ATOM'): 
	        search = re.search(r'[A-Z]*\s*\d{1,3}\s*[A-Z]{1,2}[0-9]?\s*[A-Z]{3}\s([A-Z]).*', line)
	        if search:
	            chain = search.group(1)
	            if chain in dict_simple:
	                new_chain = dict_simple[chain]
	        renamed = re.sub(r'([A-Z]*\s*\d{1,3}\s*[A-Z]{1,2}[0-9]?\s*[A-Z]{3}\s)[A-Z](.*)', r"\1" + new_chain + r"\2", line)
	        renamed_file.write("%s\n"%(renamed))  
	    elif line.startswith('TER'):
	        search = re.search(r'([A-Z]*\s*\d{1,3}\s*[A-Z]{1,3}\s)[A-Z](.*)', line)
	        if search:
	            chain = search.group(1)
	            if chain in dict_simple:
	                new_chain = dict_simple[chain]
	        renamed = re.sub(r'([A-Z]*\s*\d{1,3}\s*[A-Z]{1,3}\s)[A-Z](.*)', r"\1" + new_chain + r"\2", line)
	        renamed_file.write("%s\n"%(renamed))      
	    else:
	        renamed_file.write("%s\n"%(line))  

chain_rename(pair_file)


pair_file.close()
renamed_file.close()

# The output will be 3 files: fasta_from_pdb.fasta, alignments_results.txt, renamed_chains.pdb

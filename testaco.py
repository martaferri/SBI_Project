import Bio.PDB
import numpy
import os
import collections as col
import utilities
import copy
from functions import *
# from get_interactions_dict import dict_ids

inputs_dir = os.getcwd() + '/pairs/'
inputs_files = os.listdir(inputs_dir)

PDBparser = Bio.PDB.PDBParser()

list_of_dic = []
c=0
for pair in inputs_files:
    chain_pair_dic={}
    pdb_code = os.path.basename(os.path.splitext(pair)[0])
    pdb_filename = inputs_dir+pair
    structure = PDBparser.get_structure(pdb_code, pdb_filename)
    for model in structure:
        for chain in model:
            c+=1
            original_id = copy.copy(chain.id)
            chain.id = c
            chain_pair_dic[chain] = original_id
    list_of_dic.append(chain_pair_dic)

chain_seq_list = []
for dic in list_of_dic:
    chain_seq = tuple()
    for key in dic:
        chain_seq = (key, get_seq_from_pdbchain(key))
        chain_seq_list.append(chain_seq)


chains_eq_dic = {}
alignments_results = open("alignments_results_testaco.txt", "w")

cutoff = 0.95
for element in chain_seq_list:
    element_len = len(element[1])
    for element2 in chain_seq_list:
        element2_len = len(element2[1])
        alignments_results.write("This is the alignment between: %s and %s: " % (element, element2))
        alignments = pairwise2.align.globalxx(element[1], element2[1])
        if len(alignments) > 0:
            score_chain = alignments[0][2]
            len_max = max(element_len, element2_len)
            cutoff_chain = score_chain/len_max
            if cutoff_chain >= cutoff:
                alignments_results.write("\tThe max score is %d\n" % (alignments[0][2]))
                chains_eq_dic[element[0]] = element2[0]
                break
        else:
            alignments_results.write("\tNO ALIGNMENT\n")
alignments_results.close()

print(chains_eq_dic)









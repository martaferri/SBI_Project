#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 16:52:11 2018

@author: Aida
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 13:22:57 2018

@author: Aida
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 18 21:46:12 2018

@author: Aida
"""

import Bio.PDB
import numpy
import os
import collections as col
import utilities
import copy
from functions import *

inputs_dir = os.getcwd()+"/inputs_cm/"
inputs_files = os.listdir(inputs_dir)

PDBparser = Bio.PDB.PDBParser()


structures_list = []
chains_dic = {} #key: chain id / value: chain object
allchains = {'A': 'A', 'B': 'B', 'C':'C', 'D':'D', 'E':'A', 'F':'B', 'G':'C','H':'D'} #IMPORT
#USER!
ref_counter_user={'A':2, 'B':2, 'C':2, 'D':2}

#allchains = {'A': 'A', 'B': 'B', 'C':'C', 'D':'D', 'E':'A', 'F':'B'} #IMPORT for nucleo
#allchains = {'A': 'A', 'B': 'B', 'C': 'A', 'D': 'B'} #IMPORT for tryhemo

list_of_dic = []
for pair in inputs_files:
    chain_pair_dic={}
    pdb_code = os.path.basename(os.path.splitext(pair)[0])
    pdb_filename = inputs_dir+pair
    structure = PDBparser.get_structure(pdb_code, pdb_filename)
    for model in structure:
        for chain in model:
            for key, value in allchains.items():
                if chain.id == key:
                    chain_pair_dic[chain] = value
    list_of_dic.append(chain_pair_dic)

##############################

ref_counter_chains = copy.copy(list_of_dic[0])

current_model = [x.get_parent() for x in copy.copy(list_of_dic[0]).keys()]

n_round = 0

#define initial ref_counter
ref_counter={}
for ids, count in ref_counter_user.items():
    ref_counter[ids] = count
    for key, value in copy.copy(list_of_dic[0]).items():
        i = count
        if value == ids:
            i -= 1
            ref_counter[ids] = i



ref_chains = copy.copy(list_of_dic[0])

idx=list_of_dic.index(ref_chains)
total_list=list_of_dic[:idx]+list_of_dic[idx+1:]

run = True

while run == True:

    for alt_chains in total_list:
        print("INITIAL ref_chains: %s\t alt_chains: %s\t total_list: %s\n" % (ref_chains, alt_chains, total_list))

        print("ref_counter_chains: %s" % ref_counter_chains)
        for key, value in ref_chains.items():
            n_round += 1
            print("ROUND %d - Object: %s,  Chain: %s,  Counter: %s" % (n_round, key, value, ref_counter))

            for element in alt_chains:
                if value == alt_chains[element]:
                    print("hello")
                    print("value %s,alt_chains[element] %s"%(value,alt_chains[element]))

                    fixedchain = key
                    movingchain = element

                    alt_chains_copy = col.OrderedDict(alt_chains)
                    del alt_chains_copy[movingchain]
                    altchain = list(alt_chains_copy)[0]

                    print("\tFixed: %s, Moving: %s, Alt: %s" % (fixedchain, movingchain, altchain))

                    i_chain=alt_chains[altchain]
                    print(i_chain)

                    for ele in ref_counter:
                        print(ref_counter[i_chain])
                        if ref_counter[i_chain] != 0:
                            print("do it")

                            fixed_atoms_list = get_atoms_list(fixedchain)
                            moving_atoms_list = get_atoms_list(movingchain)
                            alt_atoms_list = get_atoms_list(altchain)

                            print("\tSUPERIMPOSING WITH:\n\t Fixed: %s, Moving: %s, Alt: %s" % (fixedchain, movingchain, altchain))

                            # creating a copy of the altchain to add it with a new id
                            list_of_ids = []
                            for element in ref_counter_chains:
                                list_of_ids.append(element.get_id())
                            for character in utilities.ascii_list:
                                if character not in list_of_ids:
                                    new_id = character
                                    break
                            new_chain = Bio.PDB.Chain.Chain(new_id)
                            altchain = copy.copy(altchain)
                            for residue in altchain.get_residues():
                                new_chain.add(residue.copy())

                            current_model[0].add(new_chain)
                            for i in current_model[0]:
                                print(i)

                            if len(fixed_atoms_list) != len(moving_atoms_list):
                                chains_pattern = refine_for_superimpose(fixedchain, movingchain)
                                fixed_pattern = chains_pattern[0]
                                moving_pattern = chains_pattern[1]

                                fixedchain = get_chain_refined(fixedchain, fixed_pattern)
                                movingchain = get_chain_refined(movingchain, moving_pattern)

                                fixed_atoms_list = get_atoms_list(fixedchain)
                                moving_atoms_list = get_atoms_list(movingchain)

                            # # find the best rotran matrices
                            # for x in range(10):
                            #     super_imposer = Bio.PDB.Superimposer()
                            #     moving_atoms_list = get_atoms_list(movingchain)
                            #     super_imposer.set_atoms(fixed_atoms_list, moving_atoms_list)
                            #     print("\t RMSD %d: %.40f" % (x, numpy.abs(super_imposer.rms)))
                            #     super_imposer.apply(movingchain.get_atoms())
                            #     super_imposer.apply(new_chain.get_atoms())

                            super_imposer = Bio.PDB.Superimposer()
                            super_imposer.set_atoms(fixed_atoms_list, moving_atoms_list)
                            print("\t RMSD : %.40f" % (numpy.abs(super_imposer.rms)))
                            super_imposer.apply(new_chain.get_atoms())

                            #TEMPORARY!!! format to save model
                            ref_chains_id = [x.id for x in current_model[0]]
                            model_name = "_".join(ref_chains_id)
                            pdb_out_filename = "%s.aligned.pdb" % model_name
                            io = Bio.PDB.PDBIO()
                            io.set_structure(current_model[0])
                            io.save(pdb_out_filename)

                            clash = check_clashes(current_model[0])
                            print("Clashes_boolean: %s" % clash)
                            if clash == True:
                                current_model[0].detach_child(new_chain.id)
                                break
                            #
                            # #format to save model
                            # ref_chains_id = [x.id for x in current_model[0]]
                            # model_name = "_".join(ref_chains_id)
                            # pdb_out_filename = "%s.aligned.pdb" % model_name
                            # io = Bio.PDB.PDBIO()
                            # io.set_structure(current_model[0])
                            # io.save(pdb_out_filename)

                            # actualize ref_conter_chains *add altchain to ref_counter_chains -> new model
                            ref_counter_chains[new_chain] = i_chain

                            print("\tRound %d ref_counter_chains: %s" % (n_round, ref_counter_chains))
                            print("\tRound %d NOT UPDATED ref_counter: %s" % (n_round, ref_counter))

                            # actualize ref_counter
                            # for key3, value3 in ref_counter_chains.items():
                            for ids, count in ref_counter.items():
                                if ids == i_chain:
                                    i = count
                                    if i != 0:
                                        i -= 1
                                        ref_counter[i_chain] = i



                            print("\tRound %d UPDATED ref_counter: %s\n" % (n_round, ref_counter))
                            break

    ref_chains = copy.copy(ref_counter_chains)
    run = False
    for ch, n in ref_counter.items():
        if n != 0:
            run = True
            break

model_filename = "1"
pdb_model_filename = "%s.aligned.pdb" % model_filename
io = Bio.PDB.PDBIO()
io.set_structure(current_model[0])
io.save(pdb_model_filename)

print("ref_counter_chains %s, ref_counter %s" % (ref_counter_chains, ref_counter))
#        ref_chains = alt_chains
#
#
print("SUMMARY\nfinal altchain: %s" % altchain)
print("final ref_counter_chains: %s" % ref_counter_chains)
print("final ref_counter: %s\n" % ref_counter)
#
print("END")

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

def get_atoms_list(chain):
    atoms = chain.get_atoms()
    atoms_list = []
    for atom in atoms:
        if atom.get_name() == 'CA':
            atoms_list.append(atom)
    return atoms_list


inputs_dir = os.getcwd()+"/inputs/"
inputs_files = os.listdir(inputs_dir)

PDBparser = Bio.PDB.PDBParser()


structures_list = []
chains_dic = {} #key: chain id / value: chain object
#allchains = {'A': 'A', 'B': 'B', 'C':'C', 'D':'D', 'E':'A', 'F':'B', 'G':'C','H': 'D'} #IMPORT
allchains = {'A': 'A', 'B': 'B', 'E':'A'} #IMPORT for nucleo
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

# print(list_of_dic)

#ref_chains = list_of_dic[0]
#ref_chains_id = [x.id for x in ref_chains]
ref_counter_chains = copy.copy(list_of_dic[0])
#
#
#ref_counter = {} # from allchains -> counter of chains
#for key, value in allchains.items():
#    i = 1
#    if value not in ref_counter:
#        ref_counter[value] = i
#    else:
#        i += 1
#        ref_counter[value] = i
#
#
current_model = [x.get_parent() for x in list_of_dic[0].keys()]
#
n_round = 0
#print("Counter: %s" % ref_counter)

ref_counter={'A':2, 'B':2}

print(list_of_dic)
for ref_chains in list_of_dic:
    idx=list_of_dic.index(ref_chains)
    total_list=list_of_dic[:idx]+list_of_dic[idx+1:]
    
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
                            for residue in altchain.get_residues():
                                new_chain.add(residue.copy())

                            current_model[0].add(new_chain)
                            for i in current_model[0]:
                                print(i)

    #get_seq etc...
    #                        if len(fixed_atoms_list) == len(moving_atoms_list):
    #                            continue
    #                        else:
                                # funcion para sacar la secuencia de cada cadena (CA!!)
                                # alinear las secuencias
                                # recorrer el alineamiento y si hay un gap, quitar el residuo de la otra cadena (por indice)
                                # pasar de 1 a 3 = fixed_seq_three
                                # mantener del pdb fixedchain original SOLO residuos de fixed_seq_three
                                # sobre fixedchain_refined hacer get_atoms_list -> para superimpose!

                            # find the best rotran matrices
                            for x in range(10):
                                super_imposer = Bio.PDB.Superimposer()
                                super_imposer.set_atoms(fixed_atoms_list, moving_atoms_list)
                                print("\t RMSD %d: %.40f" % (x, numpy.abs(super_imposer.rms)))
                            super_imposer.apply(new_chain.get_atoms())

                            ref_chains_id = [x.id for x in current_model[0]]
                            model_name = "_".join(ref_chains_id)
                            pdb_out_filename = "%s.aligned.pdb" % model_name
                            io = Bio.PDB.PDBIO()
                            io.set_structure(current_model[0])
                            io.save(pdb_out_filename)


                            # actualize ref_conter_chains *add altchain to ref_counter_chains -> new model
                            ref_counter_chains[new_chain] = i_chain

                            print("\tRound %d ref_counter_chains: %s" % (n_round, ref_counter_chains))
                            print("\tRound %d NOT UPDATED ref_counter: %s" % (n_round, ref_counter))

                            # actualize ref_counter
                            for key3, value3 in ref_counter_chains.items():
                                for ids, count in ref_counter.items():
                                    i = count
                                    if value3 == ids:
                                        if i != 0:
                                            i -= 1
                                            ref_counter[value3] = i



                            print("\tRound %d UPDATED ref_counter: %s\n" % (n_round, ref_counter))
                            break
                                    

print("ref_counter_chains %s, ref_counter %s" % (ref_counter_chains, ref_counter))
#        ref_chains = alt_chains
#
#
print("SUMMARY\nfinal altchain: %s" % altchain)
print("final ref_counter_chains: %s" % ref_counter_chains)
print("final ref_counter: %s\n" % ref_counter)
#
print("END")

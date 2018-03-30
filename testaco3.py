import Bio.PDB
import numpy
import os
import collections as col
import utilities
import copy
from functions import *
import ast

# from get_interactions_dict import dict_ids

# Setting the directory where the input files are
inputs_dir = os.getcwd() + '/test_3j7l/'
inputs_files = os.listdir(inputs_dir)

PDBparser = Bio.PDB.PDBParser()

# Creating a list of dictionaries with the pairs of the inputs. The key is the object with the id changed to numbers
# for handeling the limitation of letters (<chain1>) and the value the id from the pdb file (A)
list_of_dic_inputs = []
c = 0
for pair in inputs_files:
    chain_pair_dic = {}
    pdb_code = os.path.basename(os.path.splitext(pair)[0])
    pdb_filename = inputs_dir+pair
    structure = PDBparser.get_structure(pdb_code, pdb_filename)
    for model in structure:
        for chain in model:

            c += 1
            original_id = copy.copy(chain.id)
            chain.id = c
            chain_pair_dic[chain] = original_id
    list_of_dic_inputs.append(chain_pair_dic)

# List of tuples where the first element is the chain object and the second, its sequence
chain_seq_list = []
for dic in list_of_dic_inputs:
    chain_seq = tuple()
    for key in dic:
        chain_seq = (key, get_seq_from_pdbchain(key))
        chain_seq_list.append(chain_seq)


# Making pairwise comparisons between all the inputs to create a dictionary of equivalences
chains_eq_dic = {}
alignments_results = open("alignments_results_testaco.txt", "w")

cutoff = 0.9
cutoff_list = []
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
            info = (element[0], element2[0], cutoff_chain)
            cutoff_list.append(info)
            if cutoff_chain >= cutoff:
                alignments_results.write("\tThe max score is %d\n" % (alignments[0][2]))
                chains_eq_dic[element[0]] = element2[0]
                break
        else:
            alignments_results.write("\tNO ALIGNMENT\n")
alignments_results.close()

# Set that will contain the unique chains of the complex
unique_chains = set(chains_eq_dic.values())
print(unique_chains)

# Assignation of a letter to the unique chains to show them to the user and obtain the stoichiometry
number_to_letter = {}
for ch in unique_chains:
    number_to_letter[ch] = ""
    for character in utilities.ascii_list:
        if character not in number_to_letter.values():
            number_to_letter[ch] = character
            break

print(number_to_letter)

# Multifasta with unique chains FOR THE USER, then THE USER introduces stoichiometry
unique_chains_fasta(number_to_letter)

# Get stoichiometry from the user -> interactive or dictionary
# ref_counter_user = {'a':2, 'b':2, 'c':2, 'd':2} #'a','b','c','d' are UNIQUE chains from the values of number_to_letter dictionary (NUCLEOSOME)
# ref_counter_user = {'a': 2, 'b': 2, 'c': 2, 'd': 2, 'e': 2} # 3kuy dna
# ref_counter_user = {'a': 24} #'a','b','c','d' are UNIQUE chains from the values of number_to_letter dictionary.
# ref_counter_user = {'a': 1, 'b': 1, 'c': 1, 'd': 1, 'e': 1, 'f': 1, 'g': 1, 'h': 1, 'i': 1, 'j': 1, 'k': 1, 'l': 1, 'm': 1, 'n': 1, 'o': 1, 'p': 1, 'q': 1, 'r': 1, 's': 1, 't': 1} #4v4a
# ref_counter_user = {'a': 1, 'b': 1, 'c': 8, 'd': 1, 'e': 3, 'f': 1, 'g': 3, 'h': 3, 'i': 1, 'j': 1, 'k': 1, 'l': 1, 'm': 3, 'n': 1, 'o': 3, 'p': 1} #5vox
# ref_counter_user = {'a': 24} #2f1d
ref_counter_user = {'a': 180}  # 3j7l



# a = input("Please enter a dictionary: ")
# ref_counter_user = ast.literal_eval(a)
# print(type(ref_counter_user))
# print(ref_counter_user.keys())
#
# print("\nSTOICHIOMETRY: %s\n\n" % ref_counter_user)


list_of_dic = []
for dic in list_of_dic_inputs:
    dic_pair = {}
    for key, value in dic.items():
        for key2, value2 in chains_eq_dic.items():
            if key == key2:
                for key3, value3 in number_to_letter.items():
                    if value2 == key3:
                        dic_pair[key] = value3
    list_of_dic.append(dic_pair)

allchains = {}
for d in list_of_dic:
    for k, v in d.items():
        allchains[k] = v


#######

output_models = ("./outputs_3j7l_todo/")
if not os.path.exists(output_models):
    os.makedirs(output_models)

n_model = 0
model_number = 0


for input in list_of_dic:
    n_model += 1

    outputs = (output_models+str(n_model)+"/")
    if not os.path.exists(outputs):
        os.makedirs(outputs)

    input_idx = list_of_dic.index(input)
    ref_counter_chains = copy.copy(list_of_dic[input_idx])

    input_ids = ["a", "b"]
    i = 0
    for k, v in copy.copy(list_of_dic[input_idx]).items():
        k.id = input_ids[i]
        i += 1


    current_model = Bio.PDB.Model.Model(model_number)
    for chain in copy.copy(list_of_dic[input_idx]).keys():
        current_model.add(chain)

    list_models = [current_model]


    # current_model = [x.get_parent() for x in copy.copy(list_of_dic[input_idx]).keys()]
    # for i in current_model[model_number]:
    #     print(i)

    n_round = 0

    #define initial ref_counter
    ref_counter = {}
    for ids, count in ref_counter_user.items():
        ref_counter[ids] = count
        i = count
        for key, value in copy.copy(list_of_dic[input_idx]).items():
            if value == ids:
                i -= 1
                ref_counter[ids] = i

    ref_chains = copy.copy(list_of_dic[input_idx])
    idx = list_of_dic.index(ref_chains)
    total_list = list_of_dic[:idx]+list_of_dic[idx+1:]

    # count of ref_chains
    ref_chains_count = {}
    for v in ref_chains.values():
        if not v in ref_chains_count:
            ref_chains_count[v] = 1
        else:
            ref_chains_count[v] += 1
    print(ref_chains_count)

    # total_dict (from total_list: list with all chains except ref chains) -> total_dict_count
    total_dict = {}

    for x in total_list:
        for k, v in x.items():
            total_dict[k] = v

    total_dict_count = {}
    for ch, letter in allchains.items():
        total_dict_count[letter] = 0

    for v in total_dict.values():
        if v not in total_dict_count:
            total_dict_count[v] = 1
        else:
            total_dict_count[v] += 1

    check_input = 0
    for k, v in ref_chains.items():
        if ref_counter[v] == 0:
            check_input += 1
        elif ref_counter[v] >=1:
            if ref_chains_count[v] == 2 or total_dict_count[v] >=1:
                check_input += 1

    if check_input == 2:

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
                            print("value %s,alt_chains[element] %s" % (value, alt_chains[element]))

                            fixedchain = key
                            movingchain = element

                            alt_chains_copy = col.OrderedDict(alt_chains)
                            del alt_chains_copy[movingchain]
                            altchain = list(alt_chains_copy)[0]

                            print("\tFixed: %s, Moving: %s, Alt: %s" % (fixedchain, movingchain, altchain))

                            i_chain = alt_chains[altchain]
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
                                    used_letters_list = [x.id for x in current_model]
                                    if len(used_letters_list) == 6:
                                        list_models.append(current_model)
                                        model_number += 1
                                        current_model = Bio.PDB.Model.Model(model_number)
                                        used_letters_list = []

                                    

                                    for character in utilities.ascii_list:
                                        if character not in used_letters_list:
                                            new_id = character
                                            break
                                    new_chain = Bio.PDB.Chain.Chain(new_id)
                                    altchain = copy.copy(altchain)
                                    for residue in altchain.get_residues():
                                        new_chain.add(residue.copy())

                                    current_model.add(new_chain)
                                    for i in current_model:
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
                                    ref_chains_id = [str(x.id) for x in current_model]
                                    # model_name = "_".join(ref_chains_id)
                                    # pdb_out_filename = outputs + model_name + ".aligned.pdb"
                                    # io = Bio.PDB.PDBIO()
                                    # io.set_structure(current_model[0])
                                    # io.save(pdb_out_filename)

                                    clash = check_clashes_multi(current_model, list_models, new_chain)
                                    print("Clashes_boolean: %s" % clash)
                                    if clash == True:
                                        current_model.detach_child(new_chain.id)
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

        print("\nSAVING MODEL %d...\n\n" % n_model)

        i = 0
        for model in list_models:
            model_filename = n_model
            pdb_model_filename = outputs + str(model_filename) + "_" + str(i) + ".aligned.pdb"
            io = Bio.PDB.PDBIO()
            io.set_structure(model)
            io.save(pdb_model_filename)
            i += 1

        print("ref_counter_chains %s, ref_counter %s" % (ref_counter_chains, ref_counter))


        print("\nSUMMARY\nfinal altchain: %s" % altchain)
        print("final ref_counter_chains: %s" % ref_counter_chains)
        print("final ref_counter: %s\n" % ref_counter)
        break

print("\nEND\n")

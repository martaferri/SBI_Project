import Bio.PDB
import numpy
import os
import collections as col
import string

def get_atoms_list(chain):
    atoms = chain.get_atoms()
    atoms_list = []
    for atom in atoms:
        if atom.get_name() == 'CA':
            atoms_list.append(atom)
    return atoms_list

ascii_list = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k',
              'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', 'A', 'B', 'C', 'D', 'E', 'F',
              'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '!',
              '"', '#', '$', '%', '&', "'", '(', ')', '*', '+', ',', '-', '.', '/', ':', ';', '<', '=', '>', '?', '@',
              '[', ']', '^', '_', '`', '{', '|', '}', '~']

ascii_list2 = list(string.printable)

inputs_dir = os.getcwd()+"/inputs/"
inputs_files = os.listdir(inputs_dir)

PDBparser = Bio.PDB.PDBParser()
super_imposer = Bio.PDB.Superimposer()

structures_list = []
chains_dic = {} #key: chain id / value: chain object
allchains = {'A': 'A', 'B': 'B', 'C': 'A', 'D': 'B'} #IMPORT

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

ref_chains = list_of_dic[0]
ref_chains_id = [x.id for x in ref_chains]
ref_counter_chains = ref_chains.copy()


ref_counter = {} # from allchains -> counter of chains
for key, value in allchains.items():
    i = 1
    if value not in ref_counter:
        ref_counter[value] = i
    else:
        i += 1
        ref_counter[value] = i

# i_list=[]
# for i in current_model[0]:
#     i_list.append(i.get_parent())
current_model = [x.get_parent() for x in ref_chains.keys()]

n_round = 0
print("Counter: %s" % ref_counter)

for alt_chains in list_of_dic[1:]:
    print("INITIAL ref_chains: %s\t alt_chains: %s\n" % (ref_chains, alt_chains))

    # for each chain of ref_chains
    for key0, value0 in ref_chains.items():
        n_round += 1
        print("ROUND %d - Object: %s,  Chain: %s,  Counter: %s" % (n_round, key0, value0, ref_counter))
        print("ref_counter_chains: %s" % ref_counter_chains)
        for chain_id, counter in ref_counter.items():
            if value0 == chain_id and counter != 0:
                remaining = chain_id
                for key, value in ref_chains.items():
                    if value != remaining:
                        print("Chain %s is yet to be included by fixing %s" % (remaining, value))

                        for element in alt_chains:
                            if value == alt_chains[element]:
                                fixedchain = key
                                movingchain = element
                                break # USEFUL

                        alt_chains_copy = col.OrderedDict(alt_chains)
                        del alt_chains_copy[movingchain]
                        altchain = list(alt_chains_copy)[0]

                        print("\tFixed: %s, Moving: %s, Alt: %s" % (fixedchain, movingchain, altchain))


                        fixed_atoms_list = get_atoms_list(fixedchain)

                        moving_atoms_list = get_atoms_list(movingchain)

                        alt_atoms_list = get_atoms_list(altchain)

                        print("\tSUPERIMPOSING WITH:\n\t Fixed: %s, Moving: %s, Alt: %s" % (fixedchain, movingchain, altchain))

                        # creating a copy of the altchain to add it with a new id
                        list_of_ids = []
                        for element in ref_counter_chains:
                            list_of_ids.append(element.get_id())
                        for character in ascii_list:
                            if character not in list_of_ids:
                                new_id = character
                                break
                        new_chain = Bio.PDB.Chain.Chain(new_id)
                        for residue in altchain.get_residues():
                            new_chain.add(residue.copy())

                        current_model[0].add(new_chain)

                        if len(fixed_atoms_list) == len(moving_atoms_list):
                            continue
                        else:
                            # funcion para sacar la secuencia de cada cadena (CA!!)
                            # alinear las secuencias
                            # recorrer el alineamiento y si hay un gap, quitar el residuo de la otra cadena (por indice)
                            # actualizar las cadenas


                        # find the best rotran matrices
                        for x in range(10):
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
                        for dic in list_of_dic:
                            for key2, value2 in dic.items():
                                if altchain == key2:
                                    ref_counter_chains[new_chain] = value2


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
    ref_chains = alt_chains


print("SUMMARY\nfinal altchain: %s" % altchain)
print("final ref_counter_chains: %s" % ref_counter_chains)
print("final ref_counter: %s\n" % ref_counter)

print("END")

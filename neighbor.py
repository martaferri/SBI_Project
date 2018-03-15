
import Bio.PDB


chain_name_list = [] #will contain the names of the chains of the model
structure = Bio.PDB.PDBParser().get_structure('X', 'tot.pdb')
atom_list = Bio.PDB.Selection.unfold_entities(structure, 'A')
ns = Bio.PDB.NeighborSearch(atom_list)
chain_list = Bio.PDB.Selection.unfold_entities(structure, 'C')
for element in chain_list:
	chain_name_list.append(element.get_id()) 

for chain in chain_list:
	neighbors2 = set()
	for at in chain.get_atoms():
		center = at.get_coord()
		neighbors = ns.search(center, 2.5, level = 'C')
		for element in neighbors: 
			if element != chain:
				neighbors2.add(element)
				
	print(chain.get_id(), "\n", neighbors2)


		# print (neighbors)
# 		residue_dict[chain] = Bio.PDB.Selection.unfold_entities(neighbors, 'C') 
# 	# print (at.get_name(), residue_list[0].get_resname())
# # print(neighbors)
# print(residue_dict)

# for element in neighbors:
# 		chain = element.get_parent()
# 		print("The chain of this residue is %s\n"%(chain.get_id()))




# chain_list = Bio.PDB.Selection.unfold_entities(structure, 'C')
# print (chain_list)

# ns = Bio.PDB.NeighborSearch(chain_list)
# center = chain_list[0].get_coord()
# neighbors = ns.search(center, 5.0)
# residue_list = Bio.PDB.Selection.unfold_entities(neighbors, 'R')
# print (residue_list)

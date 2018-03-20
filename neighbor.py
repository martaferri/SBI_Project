
import Bio.PDB


chain_name_list = [] #will contain the names of the chains of the model
structure = Bio.PDB.PDBParser().get_structure('X', '1d66.pdb')
atom_list = Bio.PDB.Selection.unfold_entities(structure, 'A')
ns = Bio.PDB.NeighborSearch(atom_list)
chain_list = Bio.PDB.Selection.unfold_entities(structure, 'C')
#for element in chain_list:
#	chain_name_list.append(element.get_id()) 

for chain in chain_list:
	neighbors2 = set()
	for at in chain.get_atoms():
		center = at.get_coord()
		neighbors = ns.search(center, 1.5,  level = 'C')
		for element in neighbors: 
			if element != chain:
				neighbors2.add(element)
	if len(neighbors2) != 0:
		print(chain.get_id(), "\n", neighbors2)
	else:
		print (chain.get_id(), "\n", "No clashes were found in the chain")

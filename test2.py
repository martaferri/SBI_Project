import Bio.PDB
import numpy
import os
import collections as col

inputs_dir = os.getcwd()+"/inputs/"
inputs_files = os.listdir(inputs_dir)

PDBparser = Bio.PDB.PDBParser()
super_imposer = Bio.PDB.Superimposer()

structures_list = []
chains_dic = {} #key: chain id / value: chain object
allchains = {'A':'A', 'B':'B', 'C':'A', 'D':'B'} #IMPORT

list_of_dic=[]
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

#print(list_of_dic)  

ref_chains = list_of_dic[0]
ref_chains_id=[x.id for x in ref_chains]

ref_counter = {} #from allchains -> counter of chains
for key,value in allchains.items():
    i = 1
    if value not in ref_counter:
        ref_counter[value] = i
    else:
        i += 1
        ref_counter[value] = i

#i_list=[]
#for i in current_model[0]:
#    i_list.append(i.get_parent())

current_model = [x.get_parent() for x in ref_chains.keys()]

for alt_chains in list_of_dic[1:]:
    print("ref_chains: ",  ref_chains, "alt_chains: ",alt_chains)
    # for each chain of ref_chains
    for key, value in ref_chains.items():
#        print(key,value)
        for element in alt_chains:
#            print(element)
            if value == alt_chains[element]:
                fixedchain=key
                break
            alt_chains_copy=col.OrderedDict(alt_chains)
            movingchain = element
            del alt_chains_copy[movingchain]
            altchain=list(alt_chains_copy)[0]
            
        
            #fixedchain (from current input -> working chain)
            fixed_atoms = fixedchain.get_atoms()
            fixed_atoms_list=[]
            for atom in fixed_atoms:
                if atom.get_name() == 'CA':
                    fixed_atoms_list.append(atom)
                    
            #movingchain (same chain as working chain, but other input)
            moving_atoms = movingchain.get_atoms()
            moving_atoms_list=[]
            for atom in moving_atoms:
                if atom.get_name() == 'CA':
                    moving_atoms_list.append(atom)
                    
            #altchain (other chain in other input -> apply rotran matrix)
            alt_atoms = altchain.get_atoms()
            alt_atoms_list=[]
            for atom in alt_atoms:
                if atom.get_name() == 'CA':
                    alt_atoms_list.append(atom)
    print("fixedchain: ", fixedchain, "movingchain: ", movingchain, "altchain: ", altchain)             
    
    #actualize ref_chains *add altchain to ref_chains -> new model
    for dic in list_of_dic:
        for key2,value2 in dic.items():
            print(key2,value2)
            if altchain == key2:
                ref_chains[altchain] = value2
                
    #actualize ref_counter
    for key,value in ref_chains.items():
        for key2,value2 in ref_counter.items():
            i=value2
            if value == key2:
                i-=1
                ref_counter[value] = i
    print(ref_counter)
    
    

    
print(ref_chains)            
#        super_imposer.set_atoms(fixed_atoms_list, moving_atoms_list)
#        print(numpy.abs(super_imposer.rms))
        
#        current_model[0].add(altchain)
#        super_imposer.apply(altchain.get_atoms())
#        
#        ref_chains_id=[x.id for x in current_model[0]]
#        model_name = "_".join(ref_chains_id)
#        pdb_out_filename = "%s.aligned.pdb" % model_name
#        io=Bio.PDB.PDBIO()
#        io.set_structure(current_model[0])
#        io.save(pdb_out_filename)
        






############################################################################################
for index, ref_chains in enumerate(allchains_list):
    ref_dic={} #ref model dictionary (key: chain id / value: chain object)
    alt_dic={} #alt model dictionary (key: chain id / value: chain object)
    
    current_model = [x.get_parent() for x in ref_chains]
    ref_chains_id=[x.id for x in ref_chains]
    
    c=0
    for ch in ref_chains_id:
        for key, value in allchains.items():
            if ch == key:
                ref_dic[ref_chains[c]] = value
        c+=1

    other_chains = allchains_list.copy()
    other_chains.pop(index)

    for alt_chains in other_chains:
        alt_chains_id = [x.id for x in alt_chains] #chain_id (alt)=alt_chains_id
    i=0
    for ch in alt_chains_id:
        for key, value in allchains.items():
            if ch == key:
                alt_dic[alt_chains[i]] = value
        i+=1
    
        
    if ref_dic[ref_chains[0]] in alt_chains_id:
        fixedchain=list(ref_dic.keys())[list(ref_dic.values()).index(ref_dic[ref_chains[0]])] #retrieve key (chain object) using value (value of ref_dic -> which is equivalent to chain id)

        idx=alt_chains_id.index(ref_dic[ref_chains[0]])
        
        moving=alt_chains_id[idx]
        movingchain=list(alt_dic.keys())[list(alt_dic.values()).index(moving)]
        
        alt=alt_chains_id.copy()
        alt.pop(idx)
        alt="".join(alt)
        altchain=list(alt_dic.keys())[list(alt_dic.values()).index(alt)]
           
    elif ref_dic[ref_chains[1]] in alt_chains_id:
        fixedchain=list(ref_dic.keys())[list(ref_dic.values()).index(ref_dic[ref_chains[1]])] 
        
        idx=alt_chains_id.index(ref_dic[ref_chains[1]])
        
        moving=alt_chains_id[idx]
        movingchain=list(alt_dic.keys())[list(alt_dic.values()).index(moving)]
        
        alt=alt_chains_id.copy()
        alt.pop(idx)
        alt="".join(alt)
        altchain=list(alt_dic.keys())[list(alt_dic.values()).index(alt)]
#    print(ref_chains) 
#    print("fixed: ",fixedchain, "moving :", movingchain, "alt :", altchain)

        #fixedchain (from current input -> working chain)
        fixed_atoms = fixedchain.get_atoms()
        fixed_atoms_list=[]
        for atom in fixed_atoms:
            if atom.get_name() == 'CA':
                fixed_atoms_list.append(atom)
                
        #movingchain (same chain as working chain, but other input)
        moving_atoms = movingchain.get_atoms()
        moving_atoms_list=[]
        for atom in moving_atoms:
            if atom.get_name() == 'CA':
                moving_atoms_list.append(atom)
                
        #altchain (other chain in other input -> apply rotran matrix)
        alt_atoms = altchain.get_atoms()
        alt_atoms_list=[]
        for atom in alt_atoms:
            if atom.get_name() == 'CA':
                alt_atoms_list.append(atom)
        
        print(fixedchain == movingchain)
        super_imposer.set_atoms(fixed_atoms_list, moving_atoms_list)
        print(numpy.abs(super_imposer.rms))
        
        current_model[0].add(altchain)
        
        model_name="_".join(ref_chains_id)
        pdb_out_filename="%s_aligned.pdb" % model_name
        super_imposer.apply(altchain.get_atoms())
        
        io=Bio.PDB.PDBIO()
        io.set_structure(altchain.get_parent())
        io.save(pdb_out_filename)
        
        
        
        
        
        
        
        
            
    
    
    
    
       
    

#    allchains.setdefault(pdb_code, chain_pair)
#    allchains_ordered = col.OrderedDict(allchains)
#print(type(allchains))
#print(type(allchains_ordered))
    
#ref_chains = {}
#rc_key = list(allchains_ordered.keys())[0]
#rc_values = list(allchains_ordered.values())[0]
#ref_chains[rc_key] = rc_values
#
#
#for key, value in allchains_ordered.items():
#    print(key,value)
#    for chain in value:
#        print(chain.id)
#        if chain.id in allchains_ordered.values():
#            print ("chain in other file")
#        else:
#            print("chain not in other file")


        
            
    
    
    





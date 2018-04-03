# Biological considerations

During the development of this project, we've had to take into account some specific details in order to elaborate functions that considered the possible biological parameters and properly build the models.

### Distinguishing the type of the chain

The inputs given to the program may be of different nature: protein or nucleic acid. These types of chains are composed by different residues and atoms. To be able to work with both of them, we had to make some adaptations. 
Those are the functions involved in this matter:  
* Labeling:
    + ```check_type()```: First of all, we made a function to determine the type of the chain, this function checks if the chain has carbons alpha, if it does, we will consider it a protein, if it doesn't, a nucleic acid. Returns a string which the label.  
* Calculating distances:  
    + ```get_seq_from_pdb()```: As the residues of a protein are in three format letter, while the nucleotides are in two, or one, we had to get in a different way the sequences from the pdb file. We used a function in the protein case (```three_to_one()```), and getting the last letter by indexing in the nucleic acid case.
    + ```get_atoms_list()```: The superimposing step recquires the use of atom lists. We decided to consider just the CA atoms in protein, while P atoms in nucleic acid.
    + ```calc_distances_residues()```: To get the interactions, we compared by distances between CA atoms. As mentioned before, nucleic acids don't have this kind of atoms. We elaborated a function that renames the C1' atoms from the nucleic acid chains to CA (```adapt_chains()```). This way, we could compute the distances the same way as we would do if we only had one type of chain.  

### Superimposing

When superimposing two atoms lists, those lists have to contain the same number of atoms. Usually, arriving to this step means that the sequences of the superimposig chains are equal. However, there are cases when the chains are similar, but not exactly the same, they have different number of atoms. We handled this creating two functions:
* Comparing the chains:
    + ```refine_for_superimpose()```: The alignment of the sequences from the two chains will reveal its differences. We obtained the sequences with a function mentioned above, then made the pairwise alignment. Based on the ouput of this tool, we obtained a pattern of 1 and 0 that revealed the positions with matching or different residues. 
* Creating the new chains:
    + ```get_chain_refined()```: Once we obtained a pattern, we knew that the number of matches would be the same for both chains. This function creates a new chain that contains only the matching residues. This allowed us to finally obtain lists with the same number of atoms.
    
### Reducing the inputs

One initial step of the program is to analyse the chain interactions given by the user, and determine if there are redundant pairs to be able to reduce the inputs and build the models more fluently. This is done in the **reduce_inputs.py** script. 
For this purpose, the program makes pairwise comparisons to detect similar sequences, and if there is similarity, those chains enter a superimposition step where the structural similarity will be tested.  
To decide if the sequence similarity between the tested sequences was high enough, we set a treshold of 0.90. This was based on the distribution of the obtained scores, which showed extreme values. Setting the treshold at this point, ensured that the actual similar sequences would obtain a greater score.  

### Dealing with chain ids  

Biopython model object doesn't allow containing chains with the same id. 
To avoid this we decided to give new ids to the chains from the input to handle them without errors during the program. We used numerical annotation to deal with big numbers of chains. This solved the "same id issue" when working with objects, but to save the model in PDB format, the chain id must be of just one character.  
When a chain is added to the current model, its id is changed again. At this point, the new id is obtained from a list of ASCII characters (*ascii_list*) located in the **utilities.py** script. This last change of id allowed us to handle the saving of the created model in PDB format, but the ASCII characters list is limited, therefore if a macrocomplex is formed by more than 83 chains, we have to create a new model to continue adding chains without trouble.
To sum up, if the macrocomplex has less than 83 chains, it will be created as one model and saved in one single file, but if it doesn't, the protein will be created in more than one model and saved splitted in different files. Besides avoiding biopython errors during the program, having two or more files for one big structure avoids issues in Chimera when labeling chains.



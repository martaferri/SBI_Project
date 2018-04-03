"""Set of diverse functions to solve biological and technical problems during the analysis."""

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import Bio.PDB
import numpy
import os
import collections as col
import string
import utilities
import copy
import itertools
import classes
import re

def check_input(inputs_dir, user_dir):
    """Checks if the input given by the user is a directory."""
    if os.path.isdir(inputs_dir) == False:
        error = classes.IncorrectInputDir(inputs_dir, user_dir)
        raise error
    return "Checked."

def get_distance(list1, list2):
    """Given two lists of atom coordinates, it calculates the distance."""
    import math
    distance_result = (list1[0] - list2[0]) ** 2 + (list1[1] - list2[1]) ** 2 + (list1[2] - list2[2]) ** 2
    return math.sqrt(abs(distance_result))

def check_type(chain):
    """Function to check if the chain is a protein or a nucleic acid."""
    atoms = chain.get_atoms()
    type_chain = ""
    list_c = []
    for element in atoms:
        list_c.append(element.get_name())
    if "CA" in list_c:
        type_chain = "protein"
    else:
        type_chain = "nucleic_acid"
    return type_chain

def adapt_chain(chain):
    """Checks the type of the chain and changes the C1' to CA if it is nucleic acid.
    This will allow us to measure the distances between atoms in the same way as we do for proteins."""
    type_chain = check_type(chain)
    name = chain.id
    if type_chain == "nucleic_acid":
        new_chain = Bio.PDB.Chain.Chain(name)
        chain = copy.copy(chain)
        for residue in chain:
            new_chain.add(residue.copy())

        for residue in new_chain:
            for atom in residue:
                if atom.id == "C1'":
                    atom.id = "CA"
                    residue.add(atom.copy())
        return new_chain
    else:
        return chain

def unique_chains_fasta (number_to_letter, outputs_dir):
    """Returns a multifasta file with the sequence of the unique proteins from the total sequences of the macrocomplex."""
    fo = open(outputs_dir+"/unique_chains_fasta.mfa", "w")
    for key, value in number_to_letter.items():
        name_chain = "chain_" + value
        fo.write(">%s\n%s\n" % (name_chain, get_seq_from_pdbchain(key)))
    fo.close()
    return fo

def get_atoms_list(chain):
    """Creates a list of the atoms only taking CA or P for protein and acid nucleics, respectively.
    This list of atoms will be lately used in the superimposition process."""
    type_chain = check_type(chain)
    if type_chain == "protein":
        atom_id = "CA"
    elif type_chain == "nucleic_acid":
        atom_id = "P"
    atoms = chain.get_atoms()
    atoms_list = []
    for atom in atoms:
        if atom.id == atom_id:
            atoms_list.append(atom)
    return atoms_list

def three_to_one(three_res_list):
    """Returns sequence string (one amino acid code format). Used in the get_seq_from_pdbchain."""
    one_res_list = []

    for res in three_res_list:
        try:
            one = utilities.three_to_one[res]
            one_res_list.append(one)
        except KeyError:
            return False
    return "".join(one_res_list)

def get_seq_from_pdbchain(chain):
    """Returns the pdb amino acid or nucleic acid sequence. In case it is a protein, the sequence is returned in one
    amino acid code format. In case it is a nucleic acid, the sequence is composed by the nucleotides."""
    type_chain = check_type(chain)
    if type_chain == "protein":
        three_res_list = []
        for res in chain:
            residues_atoms = res.get_atoms()
            for atom in residues_atoms:
                if atom.get_name() == 'CA':
                    residue = atom.get_parent()
                    three_res_list.append(residue.get_resname())
        return three_to_one(three_res_list)  # three_to_one function
    else:
        nucleic_acid_res = []
        for res in chain:
            residues_atoms = res.get_atoms()
            for atom in residues_atoms:
                if atom.get_name() == 'P':
                    residue = atom.get_parent()
                    nucleic_acid_res.append(residue.get_resname())
        nucleic_acid_seq = [x[2] for x in nucleic_acid_res]
        return "".join(nucleic_acid_seq)

def refine_for_superimpose(fixedchain, movingchain):
    """Given two chains with different sequence length, returns a pattern of 0 and 1 (for each chain) based on the pairwise alignment, to
     include or exclude the residue of the original chains. This is needed for the superimposition."""
    fixedchain_seq = get_seq_from_pdbchain(fixedchain)
    movingchain_seq = get_seq_from_pdbchain(movingchain)

    alignments = pairwise2.align.globalxx(fixedchain_seq, movingchain_seq)

    fixedchain_seq = alignments[0][0]
    movingchain_seq = alignments[0][1]

    fixedchain_pattern = []
    movingchain_pattern = []

    for res1, res2 in zip(fixedchain_seq, movingchain_seq):
        if res1 == res2:
            fixedchain_pattern.append(1)
            movingchain_pattern.append(1)
        elif res1 == '-':
            movingchain_pattern.append(0)
        elif res2 == '-':
            fixedchain_pattern.append(0)

    chains_pattern = (fixedchain_pattern, movingchain_pattern)
    return chains_pattern

def get_chain_refined(chain_original, chain_pattern):
    """Creates new chain objects filtering the residues of the original chains, to get their atoms later for the
    superimposition process."""
    new_chain = Bio.PDB.Chain.Chain('X')

    for residue, pattern in zip(chain_original.get_residues(), chain_pattern):
        if pattern == 1:
            new_chain.add(residue.copy())
    return new_chain

def one_to_three(chain_refined):
    """Returns amino acid sequence of a chain (three amino acid code format) from a sequence in one amino acid code format.
    Used in refine_for_superimpose function."""
    three_res_list = []

    for res in chain_refined:
        three = utilities.one_to_three[res]
        three_res_list.append(three)
    return three_res_list

def check_clashes_multi(current_model, list_models, newchain, clash_distance_cutoff):
    """Checks clashes between the last added chain and the current model depending on the clash_distance_cutoff, which
    is set by the user (default 1.2A)."""
    neighbors2 = set()

    atom_list = Bio.PDB.Selection.unfold_entities(current_model, 'A')

    for model in list_models:
            atom_list2 = Bio.PDB.Selection.unfold_entities(model, 'A')
            for element in atom_list2:
                atom_list.append(element)

    ns = Bio.PDB.NeighborSearch(atom_list)

    for at in newchain.get_atoms():
        center = at.get_coord()
        neighbors = ns.search(center, clash_distance_cutoff, level='C')
        for element in neighbors:
            if element != newchain:
                neighbors2.add(element)

    if len(neighbors2) != 0:
        return True
    else:
        return False

def create_tempPDB(outputs):
    """Erases the acid nucleic chains of the pdb files that will be used in the refining process, once the model
    has been built (to avoid Modeller errors)."""
    temp_dir = outputs + "temp/"
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    files_list = []
    models_list = os.listdir(outputs)
    for file in models_list:
        if file.endswith(".pdb") or file.endswith(".B"):
            files_list.append(file)

    for element in files_list:
        pdb = open(outputs + element, 'r')
        fo = open(temp_dir + "mod_" + element, "w")
        for line in pdb:
            if line.startswith('ATOM'):
                line = line.strip()
                search = re.search(r'[A-Z]*\s*\d{1,5}\s*[A-Z]{1,3}[0-9]?\'?\s*([A-Z]{1,3}).*', line)
                residue = search.group(1)
                if residue in utilities.nucleic_list:
                    continue
                else:
                    fo.write("%s\n" % (line))
        fo.close()
        pdb.close()
    return temp_dir




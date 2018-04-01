import os
import sys
import re

def blockPrint(file = os.devnull):
    sys.stdout = open(file, 'w')

def enablePrint():
    sys.stdout = sys.__stdout__

blockPrint()

from modeller import *
from modeller.scripts import complete_pdb
from modeller.optimizers import conjugate_gradients, molecular_dynamics, actions

files_list = []
dict_finalen = {} 
list_finalen = []

models_dir = "./models/"
dl = os.listdir(models_dir)

for file in dl:
	if file.endswith(".pdb"):
		files_list.append(file)

enablePrint()
print ("Analyzing", len(files_list), "models...")
blockPrint()

if not os.path.exists("./models/optimization_results/refined_models/"):
    os.makedirs("./models/optimization_results/refined_models/")
if not os.path.exists("./models/optimization_results/stats/"):
    os.makedirs("./models/optimization_results/stats/")
if not os.path.exists("./models/optimization_results/dope_profile/"):
    os.makedirs("./models/optimization_results/dope_profile/")
if not os.path.exists("./models/optimization_results/log_files/"):
    os.makedirs("./models/optimization_results/log_files/")  
if not os.path.exists("./models/optimization_results/restraints/"):
    os.makedirs("./models/optimization_results/restraints/")

env = environ()
env.io.atom_files_directory = ['../atom_files']
env.edat.dynamic_sphere = True

env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

for file in files_list:
    blockPrint(file = "./models/optimization_results/log_files/" + str(file[:-4]) + '.log')
    print ("\n", "##########", file, "##########", "\n")
    code = "./models/" + str(file)
    mdl = complete_pdb(env, code)
    atmsel = selection(mdl)
    mdl.restraints.make(atmsel, restraint_type='improper', spline_on_site=False)
    mdl.restraints.write(file= "./models/optimization_results/restraints/" + str(file)+'.rsr')
    cg = conjugate_gradients(output='REPORT')
    
    trcfil = open("./models/optimization_results/stats/" + str(file) +'.stats', 'w')
    
# Run CG on the all-atom selection; write stats every 5 steps
    cg.optimize(atmsel, max_iterations=20, actions=actions.trace(10, trcfil))
    mpdf2 = atmsel.energy()
    mdl.write(file="./models/optimization_results/refined_models/" + str(file)+'.B')
    print ("Final energy: ", mpdf2[0])
    dict_finalen[file] = mpdf2[0]
    list_finalen.append(mpdf2[0])
	
enablePrint()

print ("The minimum final energy for all the models is:", min(list_finalen), "and corresponds to:", list(dict_finalen.keys())[list(dict_finalen.values()).index(min(list_finalen))] + '.B')

blockPrint()

import pylab
import modeller

def get_profile(profile_file, seq):
    """Read `profile_file` into a Python array, and add gaps corresponding to
       the alignment sequence `seq`."""
    # Read all non-comment and non-blank lines from the file:
    f = open(profile_file, 'r')
    vals = []
    for line in f:
        if not line.startswith('#') and len(line) > 10:
            spl = line.split()
            vals.append(float(spl[-1]))
    # Insert gaps into the profile corresponding to those in seq:
    for n, res in enumerate(seq.residues):
        for gap in range(res.get_leading_gaps()):
            vals.insert(n, None)
    # Add a gap at position '0', so that we effectively count from 1:
    vals.insert(0, None)
    return vals

#Create a PDB files without nucleic acids. Stored in temp. 
files_list = []
dict_nonres = ["A", "G", "T", "U", "C", "DA", "DG", "DT", "DC"]
models_dir = "./models/"
dl = os.listdir(models_dir)
for file in dl:
    if file.endswith(".pdb"):
        files_list.append(file)


for file in files_list: 
    i_pdb = open("./models/" + file, 'r')
    fo = open("./models/temp/" + "new" + file + ".B", "w")
    for line in i_pdb:
        if line.startswith('ATOM'):
            line = line.strip()
            search = re.search(r'[A-Z]*\s*\d{1,5}\s*[A-Z]{1,3}[0-9]?\'?\s*([A-Z]{1,3})\s[a-zA-Z1-9].*', line)
            residue = search.group(1)
            if residue in dict_nonres: 
                continue
            else: 
                fo.write("%s\n" % (line))
    fo.close()

#Obtention of files necessary to plot -- models refined
models_dir = "./models/temp/"
env.io.atom_files_directory = [models_dir]
mdl_list = []
aln = modeller.alignment(env)
code_list = []

for file in os.listdir(models_dir):
    if file.endswith(".pdb.B"):
        print ("hola")
        mdl = modeller.model(env)
        mdl.read(file = file)
        code = str(file)
        code_list.append(code)
        s = selection(mdl)
        s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file= "./models/temp/" + code + '.profile',
              normalize_profile=True, smoothing_window=15)
        mdl_list.append(mdl)
        aln.append_model(mdl, align_codes = code, atom_files = code)
        aln.write(file='build_profile_ref.ali', alignment_format='PIR')


models_dir = "./models/temp/"
env.io.atom_files_directory = [models_dir]
mdl_nr_list = []
aln_nr = modeller.alignment(env)
code_list_nr = []


for file in os.listdir(models_dir): 
    if file.endswith(".pdb"):
        mdl_nr = modeller.model(env)
        mdl_nr.read(file = file)
        code = str(file)
        print (code)
        code_list_nr.append(code)
        t = selection(mdl_nr)
        t.assess_dope(output='ENERGY_PROFILE NO_REPORT', file= "./models/temp/" + code + '.profile',
              normalize_profile=True, smoothing_window=15)
        mdl_nr_list.append(mdl_nr)
        aln_nr.append_model(mdl_nr, align_codes = code, atom_files = code)
        aln_nr.write(file='build_profile_notref.ali', alignment_format='PIR')


if len(mdl_nr_list) == len(mdl_nr_list):
    for a, b, c, d in zip(mdl_nr_list, mdl_list, code_list_nr, code_list):
        model1 = get_profile("./models/temp/" + c + ".profile", aln_nr[str(c)])
        model2 = get_profile("./models/temp/" + d + ".profile", aln[str(d)])
        pylab.figure(1, figsize=(30,18))
        pylab.xlabel('Alignment position', fontsize = 20)
        pylab.ylabel('DOPE per-residue score', fontsize = 20)
        pylab.plot(model1, color='red', linewidth=2, label='Optimized model')
        pylab.plot(model2, color='green', linewidth=2, label='Model')
        pylab.legend(fontsize = 20)
        pylab.savefig("./models/optimization_results/dope_profile/" + c + 'dope_profile.jpg', dpi=100)
        pylab.close()


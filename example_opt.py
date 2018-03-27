import os
import sys

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
env = environ()
env.io.atom_files_directory = ['../atom_files']
env.edat.dynamic_sphere = True

env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

models_dir = "./outputs/"
dl = os.listdir(models_dir)

for file in dl:
	if file.endswith(".pdb"):
		files_list.append(file)

enablePrint()
print ("Analyzing", len(files_list), "models...")
blockPrint()

if not os.path.exists("./outputs/optimization_results/refined_models/"):
    os.makedirs("./outputs/optimization_results/refined_models/")
if not os.path.exists("./outputs/optimization_results/stats/"):
    os.makedirs("./outputs/optimization_results/stats/")
if not os.path.exists("./outputs/optimization_results/restraints/"):
    os.makedirs("./outputs/optimization_results/restraints/")
if not os.path.exists("./outputs/optimization_results/log_files/"):
    os.makedirs("./outputs/optimization_results/log_files/")  


for file in files_list:
	blockPrint(file = "./outputs/optimization_results/log_files/" + str(file[:-4]) + '.log')
	print ("\n", "##########", file, "##########", "\n")
	code = "./outputs/" + str(file)
	mdl = complete_pdb(env, code)
#	mdl.write(file="./outputs/optimization_results/" + str(file)+'.ini')

# Select all atoms:
	atmsel = selection(mdl)
#	mpdf0 = atmsel.energy()
	
# Generate the restraints:
	mdl.restraints.make(atmsel, restraint_type='improper', spline_on_site=False)
	mdl.restraints.write(file= "./outputs/optimization_results/restraints/" + str(file)+'.rsr')
	
# Create optimizer objects and set defaults for all further optimizations
	cg = conjugate_gradients(output='REPORT')
# md = molecular_dynamics(output='REPORT')

# Open a file to get basic stats on each optimization
	trcfil = open("./outputs/optimization_results/stats/" + str(file) +'.stats', 'w')

# Run CG on the all-atom selection; write stats every 5 steps
	cg.optimize(atmsel, max_iterations=20, actions=actions.trace(10, trcfil))
	mpdf2 = atmsel.energy()
	mdl.write(file="./outputs/optimization_results/refined_models/" + str(file)+'.B')
	# print ("Initial energy: ", mpdf0[0])
	print ("Final energy: ", mpdf2[0])
	dict_finalen[file] = mpdf2[0]
	list_finalen.append(mpdf2[0])

enablePrint()

print ("The minimum final energy for all the models is:", min(list_finalen), "and corresponds to:", list(dict_finalen.keys())[list(dict_finalen.values()).index(min(list_finalen))] + '.B')

# Run MD; write out a PDB structure (called '1fas.D9999xxxx.pdb') every
# 10 steps during the run, and write stats every 10 steps
# md.optimize(atmsel, temperature=300, max_iterations=50,
#             actions=[actions.write_structure(10, code+'.D9999%04d.pdb'),
#                      actions.trace(10, trcfil)])
# Finish off with some more CG, and write stats every 5 steps
# cg.optimize(atmsel, max_iterations=20, actions=[actions.trace(10, trcfil)])

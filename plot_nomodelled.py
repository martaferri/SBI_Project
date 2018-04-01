import os
import sys

def blockPrint(file = os.devnull):
	"""Blocks the printing in the terminal window"""
	sys.stdout = open(file, 'w')

def enablePrint():
	"""Enables the printing in the terminal window -- normally blocked by blockPrint() function"""
	sys.stdout = sys.__stdout__

blockPrint()

import modeller
import re
import pylab

from modeller import *
from modeller.scripts import complete_pdb
from modeller.optimizers import conjugate_gradients, molecular_dynamics, actions

def get_profile(profile_file, seq):
    """Read `profile_file` into a Python array, and add gaps corresponding to
       the alignment sequence `seq`."""
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

if not os.path.exists("./models/temp/"):
    os.makedirs("./models/temp/")
if not os.path.exists("./models/dope_profiles"):
	os.makedirs("./models/dope_profiles/")

files_list = []
dict_nonres = ["A", "G", "T", "U", "C", "DA", "DG", "DT", "DC"]
models_dir = "./models/"
dl = os.listdir(models_dir)
for file in dl:
	if file.endswith(".pdb"):
		files_list.append(file)

#Create a PDB files without nucleic acids. Stored in temp. 
for file in files_list: 
	i_pdb = open("./models/" + file, 'r')
	fo = open("./models/temp/" + "new" + file, "w")
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

# Obtention of DOPE profiles and alignments required for the plot. 
flist = []
mdir = "./models/"
models_dir = "./models/temp"
env = environ()
env.io.atom_files_directory = [models_dir]
dl = os.listdir(models_dir)

for file in dl:
	if file.endswith(".pdb") & file.startswith("new"):
		flist.append(file)

aln = modeller.alignment(env)

for file in flist:
	mdl = modeller.model(env)
	code = str(file)
	print (code)
	mdl.read(file = code, model_segment=('FIRST:@', 'END:'))
	aln.append_model(mdl, align_codes = code, atom_files = code)
	t = selection(mdl)
	t.assess_dope(output='ENERGY_PROFILE NO_REPORT', file= "./models/" + code + '.profile', normalize_profile=True, smoothing_window=15)
	model = get_profile("./models/" + file + ".profile", aln[str(file)])
	

	pylab.figure(1, figsize=(20,12))
	pylab.xlabel('Alignment position', fontsize = 20)
	pylab.ylabel('DOPE per-residue score', fontsize = 20)
	pylab.plot(model, color='green', linewidth=3, label= file[3:-4])
	pylab.legend(fontsize = 15)
	pylab.savefig("./models/dope_profiles/" + file + 'dope_profile.jpg', dpi=100)
	pylab.close()





import os
import sys
import modeller
import re
import pylab

from modeller import *
from modeller.scripts import complete_pdb
from modeller.optimizers import conjugate_gradients, molecular_dynamics, actions


def blockPrint(file=os.devnull):
    """Blocks the printing in the terminal window."""
    sys.stdout = open(file, 'w')


def enablePrint():
    """Enables the printing in the terminal window -- normally blocked by blockPrint() function."""
    sys.stdout = sys.__stdout__


def get_profile(profile_file, seq):
    """Read 'file_dope' into a Python array, and add gaps corresponding to the alignment sequence 'seq'."""
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


def DOPE_profiles_maker(temp_dir, outputs):
    """Creates a DOPE profile plot (.jpg) from a macrocomplex (.pdb), which has no acid nucleic chains using Modeller."""
    flist = []
    env = environ()
    env.io.atom_files_directory = [temp_dir]
    dl = os.listdir(temp_dir)

    for file in dl:
        if file.startswith("mod"):
            flist.append(file)

    aln = modeller.alignment(env)

    for file in flist:
        mdl = modeller.model(env)
        code = str(file)
        mdl.read(file=code, model_segment=('FIRST:@', 'END:'))
        aln.append_model(mdl, align_codes=code, atom_files=code)
        t = selection(mdl)
        file_dope = outputs + code + '.profile'
        t.assess_dope(output='ENERGY_PROFILE NO_REPORT', file=file_dope, normalize_profile=True, smoothing_window=15)
        model = get_profile(file_dope, aln[str(file)])
        pylab.figure(1, figsize=(20, 12))
        pylab.xlabel('Alignment position', fontsize=20)
        pylab.ylabel('DOPE per-residue score', fontsize=20)
        pylab.plot(model, color='green', linewidth=3, label=file[3:-4])
        pylab.savefig(outputs + file[:-4] + '.dope_profile.jpg', dpi=100)
        pylab.close()

    path_img = outputs + file[:-4] + '.dope_profile.jpg'
    return("DOPE profile plot for model created here:\n  %s\n" % (path_img))

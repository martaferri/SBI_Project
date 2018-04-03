"""Error class"""

# python modules
import Bio.PDB
import numpy
import os
import gzip
import re
import collections as col
import copy
import argparse
import tarfile
import ast

# our modules
from functions import *
from classes import *
from reduce_inputs_func import *
import utilities

class IncorrectInputDir(ValueError):
    """Checks if the input variable is a directory or a compressed directory (.tar.gz)."""
    def __init__(self, input_dir, user_dir):
        self.input_dir = input_dir
        self.user_dir = user_dir

    def __str__(self):
        return "%s is not neither a directory nor a compressed directory (.tar.gz). Please, try again with a proper directory." % (self.user_dir)


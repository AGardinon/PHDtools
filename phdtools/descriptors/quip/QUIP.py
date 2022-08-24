# -------------------------------------------------- #
# Descriptors tools 
#
#
# AUTHOR: Andrea Gardin
# -------------------------------------------------- #

import numpy as np
import quippy
from phdtools.computes.misc import parse_config
from .QUIPtools import quipSOAP_string

class Descriptors:
    def __init__(self, type, params_file):
        if type(params_file) == dict:
            self.params_dict = self.params_file
        elif type(params_file) == str:
            self.params_dict == parse_config(filename=self.params_file)[self.type]

    def calc_desc(self, ase_db, frames, save_file=None):
        
        pass

    def _calc(self):
        print(f"Computing the descriptor ...\n")
        print()
        dscr_obj = quippy.descriptors.Descriptor()
        pass
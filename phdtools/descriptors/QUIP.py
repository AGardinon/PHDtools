# -------------------------------------------------- #
# Descriptors tools 
#
#
# AUTHOR: Andrea Gardin
# -------------------------------------------------- #

import numpy as np
import quippy
from phdtools.computes.misc import parse_config
from .quipTools import quipSOAP_string

class Descriptors:
    def __init__(self, desc_type, params_file):
        self.desc_type = desc_type
        if type(params_file) == dict:
            self.params_dict = params_file
        elif type(params_file) == str:
            print('Loading the file.')
            self.params_dict = parse_config(filename=params_file)[self.desc_type]
        self.params_str = quipSOAP_string(self.params_dict, desc_type=self.desc_type)
        self.Obj = quippy.descriptors.Descriptor(self.params_str)

    def calc_desc(self, ase_db, save_file=None):
        print(f"Computing the descriptor ...\n")
        print(f"{self.params_str}")
        desc_db = self._calc(ase_db)
        if type(save_file) == str:
            np.save(save_file, desc_db)
        else:
            raise ValueError("You need to insert a string type name.")
        return 

    def _calc(self, ase_db):
        return self.Obj.calc_descriptor(ase_db)

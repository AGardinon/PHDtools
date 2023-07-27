# -------------------------------------------------- #
# ASE tools - trajectories handler tool in ASE
#
#
# AUTHOR: Andrea Gardin
# -------------------------------------------------- #

import os
import json
import numpy as np
from ase.io import read, write
from .atomstools import *
from phdtools.computes import misc, traj

# -------------------------------------------------- #
# --- trajtools 

class ASEProjectManager:

    def func():
        pass

class Universe:

    @misc.my_timer
    def __init__(self, 
                 projectName: str, 
                 trajPath: str,
                 **kwargs):
        """ASE universe tool (similar to the MDA one).
        Optional kwarg is rcutCorrection (dict), but can be set later.

        :param projectName: Name of the project (aka the system).
        :type projectName: str
        :param trajPath: Path for the traj file.
        :type trajPath: str
        """
        # getting project info
        self.projectName = projectName
        self.trajPath = trajPath
        self.rcutCorrection = None
        self.mol_info_dict = dict()
        # ---
        # init the project
        self._get_config
        try:
            self.rcutCorrection = kwargs['rcutCorrection']
            print(f"rcut correction: {self.rcutCorrection}\n")
        except:
            print("!!! Warning !!!\n"
                  "`rcutCorrection` not set, default values will be used!"
                  )
            # creates the default one
            self.rcutCorrection = {k:1.0 for k in self.symbols}
            print(f"rcutCorrection = {self.rcutCorrection}")
            
    @property
    def _get_config(self):
        """Reads the first frame to get
        information about the system, 
        effectively building the Universe.
        """
        print("Gathering the Universe ...\n")
        self._at0 = read(filename=self.trajPath, index='0')
        self.symbols = np.unique(self._at0.symbols)
        print(f"Total atoms: {len(self._at0.symbols)}\n"
              f"Atom types: {self.symbols}\n"
              )


    def set_rcut_correction(self, rcut_dict: dict) -> dict:
        """Allow the setting of a rcut dictionary for
        the atoms types in the system.

        :param rcut_dict: LJ radius correction parameter (1.0 = unchanged).
        :type rcut_dict: dict
        :raises ValueError: Only defined if the atoms are in the system.
        :return: Dictionary containing kw: atoms and args: rcut correction.
        :rtype: dict
        """
        if any(rcut_dict.keys()) == any(rcut_dict.keys()):
            self.rcutCorrection = rcut_dict
            return print(f"Updated rcutCorrection = {self.rcutCorrection}\n")
        else:
            raise ValueError("Given atom types does not match with the system types.\n")

    @misc.my_timer
    def find_molecs(self) -> None:
        self.mol_dict = get_chemFormulas(self._at0, 
                                         fct=self.rcutCorrection)
        print(f"Molecules found: {self.mol_dict}")
        pass
    

    def set_mol_names(self, mol_name_list: list) -> dict:
        try:
            keys = list(self.mol_dict.keys())
        except:
            raise ValueError("The molecule dict is not set.\n"
                             "(see `Universe.find_molecs` function)")
        self.mol_names = {k:f for k,f in zip(keys, mol_name_list)}
        # updating the self.mol_info_dict
        self.mol_info_dict['names'] = self.mol_names
        return self.mol_names

    @misc.my_timer
    def get_mol_info(self):
        # """
        # Get the molecules IDs for the frame configuration.
        # (i.e., an array with N_mol index for each different molecule)
        # The molIDs is also added in place to the at.array dict
        # """
        # self.molIDs = get_molIDs(at=self._at0,
        #                          fct=self.rcutCorrection)
        # return self.molIDs
        pass
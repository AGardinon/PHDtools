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
        self.info_dict = dict()
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
    def find_molecs(self, 
                    mol_name: list = None) -> None:
        """Finds the molecules as whole defined by the 
        current LJ rcut.
        """
        self.mol_dict = get_chemFormulas(self._at0, 
                                         fct=self.rcutCorrection)
        print(f"Uniques molecules found: {len(self.mol_dict)}")
        if mol_name:
            self.mol_dict = self.set_mol_names(mol_name=mol_name)
        print(f"Molecules found: {self.mol_dict}")
        # update the self.info_dict()
        self.info_dict['names'] = self.mol_dict
        pass
    
    # @property
    def set_mol_names(self, 
                      mol_name: list) -> dict:
        """Set new specific names for molecules found by
        Universe.find_molecs().

        :param mol_name_list: List of names in the same order of appearence.
        :type mol_name_list: list
        :raises ValueError: The molecules needs to be defined before this operation.
        :return: Dictionary with as args: custom mol names
        :rtype: dict
        """
        if self.mol_dict:
            if len(self.mol_dict) == len(mol_name):
                for key,val in zip(self.mol_dict.keys(),mol_name):
                    self.mol_dict[key] = val
                self.info_dict['names'] = self.mol_dict
                return self.mol_dict
            else:
                raise ValueError("The list of molecules provided is not compatible,"
                                 "(the system has {len(self.mol_dict)} molecs)")
        else:
            raise ValueError("The molecules have to be found before setting the name.")

    @property
    def get_mol_info(self):
        print("Computing MolIDs\t...")
        self.molIDs = get_molIDs(db=self._at0,
                                 fct=self.rcutCorrection)       
        print("Computing MolSymbols\t...")
        for mol in np.unique(self.molIDs):
            pass
        pass
# -------------------------------------------------- #
# ASE tools - trajectories handler tool in ASE
#
#
# AUTHOR: Andrea Gardin
# -------------------------------------------------- #

import os
import json
import numpy as np
from typing import Union, List
from ase.io import read, write
from .atomstools import *
from ..computes import misc, traj

# -------------------------------------------------- #
# --- Universe 

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
        Optional kwargs are:
        - rcutCorrection: dict (Default=1.0)
        - moleculeNames: list (Default='molXXX')

        :param projectName: Name of the project (aka the system).
        :type projectName: str
        :param trajPath: Path for the traj file.
        :type trajPath: str
        """
        # --- Universe info
        # -
        self.projectName = projectName
        self.trajPath = trajPath
        self.rcutCorrection = None
        self.moleculeNames = None
        # --- init the project
        # -
        self._get_config
        try:
            self.rcutCorrection = kwargs['rcutCorrection']
            print(f"rcut correction: {self.rcutCorrection}\n")
        except:
            print("!!! Warning: `rcutCorrection` not set !!!\n"
                  "Default values will be used, "
                  "this might cause artifacts in the molecules detection."
                  )
            # creates the default one
            self.rcutCorrection = {k:1.0 for k in self.symbols}
            print(f"rcutCorrection = {self.rcutCorrection}\n")
        # -
        print("Searching for molecules in the system ...\n")
        try:
            self.moleculeNames = kwargs['moleculeNames']
            print(self.moleculeNames)
        except:
            print("!!! Warning: `moleculeNames` not set !!!\n"
                  "Default molecules names will be used\n")
        self.find_molecs(mol_name=self.moleculeNames)
        self.get_mol_info
        # - Update the dictionary
        self.projectDictionary = dict(
            projectName = self.projectName,
            trajPath = self.trajPath,
            rcutCorrection = self.rcutCorrection,
            moleculeNames = self.moleculeNames,
            moleculeFormulas = self.moleculeFormulas,
            molIDs = self.molIDs,
            molSym = self.molSym
        )
        # -
        print("<end>")
            
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
        pass


    def set_rcut_correction(self, 
                            rcut_dict: dict) -> dict:
        """Allow the setting of a rcut dictionary for
        the atoms types in the system.

        :param rcut_dict: LJ radius correction parameter (1.0 = unchanged).
        :type rcut_dict: dict
        :raises ValueError: Only defined if the atoms are in the system.
        :return: Dictionary containing kw: atoms and args: rcut correction.
        :rtype: dict
        """
        if set(rcut_dict.keys()) == set(self.symbols):
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
        self.moleculeFormulas = get_chemFormulas(self._at0, 
                                                 fct=self.rcutCorrection)
        self.moleculeNames = list(self.moleculeFormulas.values())
        print(f"Uniques molecules found: {len(self.moleculeFormulas)}")
        if mol_name:
            self.moleculeFormulas = self.set_mol_names(mol_name=mol_name)
        print(f"Molecules found: {self.moleculeFormulas}")
        pass
    

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
        self.moleculeNames = mol_name
        if self.moleculeFormulas:
            if len(self.moleculeFormulas) == len(mol_name):
                for key,val in zip(self.moleculeFormulas.keys(),mol_name):
                    self.moleculeFormulas[key] = val
                return self.moleculeFormulas
            else:
                raise ValueError("The list of molecules provided is not compatible,"
                                 "(the system has {len(self.moleculeFormulas)} molecs)")
        else:
            raise ValueError("The molecules have to be found before setting the name.")

    @property
    def get_mol_info(self):
        """TODO
        """
        # ---
        print("Computing MolIDs ...")
        self.molIDs = get_molIDs(at=self._at0,
                                 fct=self.rcutCorrection)
        # ---
        print("Computing MolSymbols ...")
        self.molSym = get_molSym(at=self._at0,
                                 molIDs=self.molIDs,
                                 mol_name=self.moleculeFormulas)
        pass


# -------------------------------------------------- #
# --- trajtools

# should be improved with proper exception creation!
def frame_type_checker(frame_value : Union[tuple, list]) -> None:
    try:
        if isinstance(frame_value, (tuple, str)):
            if frame_value == str and not frame_value == 'all':
                raise ValueError("ValueError: string values supported is only 'all'.")
        else:
            raise ValueError("ValueError: only tuple and string values are accepted.")
        
    except ValueError as error:
        print("!!! " + repr(error)) 


class ASEtraj(Universe):
    """Child class for handling trajectories and doing operation whitin them.

    :param Universe: sets the Ase Universe as a baseline to compute trajectory 
    related quantities.
    :type Universe: class
    """


    def __init__(self, 
                 projectName: str, 
                 trajPath: str, 
                 frames : Union[tuple, list] = 'all',
                 **kwargs):
        super().__init__(projectName, trajPath, **kwargs)

        # - specific attributes
        _ = frame_type_checker(frame_value=frames)
        self._frames = frames
        pass

    @property
    def frames(self) -> Union[tuple, list]:
        """Get the frames used for the analysis.

        :return: Frame value range.
        :rtype: Union[tuple, list]
        """
        return self._frames
    
    @frames.setter
    def frames(self, 
               value: Union[tuple, list]) -> None:
        """Allow to set different frames values for the analysis.

        :param value: Frame value range.
        :type value: Union[tuple, list]
        """
        _ = frame_type_checker(frame_value=value)
        self._frames = value
        pass
    
    @property
    def _read(self) -> List[ase.ase.Atoms]:
        if isinstance(self._frames, tuple):
            try:
                b,e,s = self._frames
            except:
                s = 1
                b,e = self._frames
            print("Reading traj:\n"
                    f"Begin: {b} | End: {e} | Stride: {s}")
            return read(self.trajPath, index=f'{b}:{e}:{s}')

        else:
            return read(self.trajPath, index=f':')
            
        
    # @property
    # def output_dir(self):
    #     return self._outputDir
    
    # @output_dir.setter
    # def output_dir(self, dir):
    #     self._outputDir = dir
    #     pass
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
from phdtools.computes import misc

# -------------------------------------------------- #
# --- trajtools 

# -- base class
class BaseUni:

    DBFILE = os.path.dirname(__file__) + "/chemFormulaToName.json"

    @property
    def _project_db(self):
        with open(self.DBFILE) as jsonfile:
            project_db = json.load(jsonfile)
        return project_db

    def __init__(self, 
                 projectName, 
                 trajPath, 
                 rcut_correction=1.):
        self.projectName = projectName
        self.trajPath = trajPath
        self.rcutCorrection = rcut_correction
        # little check
        if self.projectName in self._project_db.keys():
            self.mol_chemName_db = self._project_db[self.projectName]
        elif not self.projectName in self._project_db.keys():
            print(f"*Warning*: {self.projectName} not an avaliable project.")
            print(f"You need to ad it to the project database.\n(default @ {self.DBFILE})")

    def addNewProject(self, save_to_default=False, save_name=None):
        if self.projectName in self._project_db.keys() and self.mol_chemName_db:
            print(f"Project already exists:\n{self.mol_chemName_db}")
        elif not self.projectName in self._project_db.keys():
            print("Adding new project.")
            self.mol_chemName_db = self.__add_new_project()
        # once created it is possible to add it
        if save_to_default:
            new_project_dict = self._project_db.copy()
            if save_name:
                new_project_dict[save_name] = self.mol_chemName_db
            elif not save_name:
                new_project_name = "NewProject"+misc.todayDate()
                new_project_dict[new_project_name] = self.mol_chemName_db
            # update .json file
            with open(self.DBFILE,'w+') as jsonfile:
                json.dump(new_project_dict, jsonfile)
            print(f"New project added.\n{self._project_db}")
        elif not save_to_default:
            print("New project not added to the defaul database.")

    def __add_new_project(self):
        # care for list of ase.Atoms instead of just the Atoms obj
        at0 = read(self.trajPath, index='0:1')[0]
        return get_chemFormulas(at=at0, fct=self.rcutCorrection)

    def deleteProject(self, project_name):
        try:
            new_project_dict = self._project_db.copy()
            del new_project_dict[project_name]
            # update .json file
            with open(self.DBFILE,'w+') as jsonfile:
                json.dump(new_project_dict, jsonfile)
            print(f"Project {project_name} deleted.\n{self._project_db}")
        except:
            raise NameError("Can not erase the project because is not there.")


# -------------------------------------------------- #
# -- main class
# Attributes:
# + readd: frame, chunk, whole
# + get COM-CG traj for the syst mol
# + unwrap traj
# + Z number shift (to cheat the atoms difference)

class Universe(BaseUni):
    """ASE traj handles class"""

    @property
    def _get_info(self):
        print("\nGathering the Universe ...")
        traj_info = dict()
        # general info on the atomistic traj
        dummy_frame_atoms = self._read(frames=0)
        traj_info['atmsZ'] = dummy_frame_atoms[0].numbers
        # general infro on the molecular traj
        dummy_frame_mols = extract_molInfo(db=dummy_frame_atoms,
                                           mol_chemName_db=self.mol_chemName_db, 
                                           fct=self.rcutCorrection)
        traj_info['molSym'] = dummy_frame_mols[0].arrays['molSym']
        traj_info['molID'] = dummy_frame_atoms[0].arrays['molID']
        return traj_info

    def __init__(self, projectName, trajPath, rcut_correction=1.):
        super().__init__(projectName, trajPath, rcut_correction)
        self.infoDict = self._get_info

    def read(self, frames, Zshift_tuple=None, save_file=None):
        print("\nReading trajectory ...")
        ase_db = self._read(frames=frames)
        if Zshift_tuple:
            new_Znumbers = self._shift_Znumbers(Zshift_tuple=Zshift_tuple)
            for snap in tqdm(ase_db, desc='Applying Z shift'):
                snap.numbers = new_Znumbers
        if save_file and type(save_file) == str:
            self.save_traj(ase_db=ase_db, file_name=save_file, frames=frames)
        return ase_db

    def read_COM(self, frames, save_traj=False):
        print("\nReading trajectory and getting COM ...")
        ase_COM_db = self._get_molCOM(ase_db=self._read(frames=frames),
                                      molID=self.infoDict['molID'], 
                                      molSym=self.infoDict['molSym'])
        if save_traj:
            self._save_traj(ase_db=ase_COM_db, name='COM_traj', frames=frames)
        return ase_COM_db


    def _read(self, frames):
        if type(frames) == int:
            ase_dbtraj_list = read(self.trajPath, index=f'{frames}:{frames+1}')
        elif type(frames) == tuple:
            if len(frames) == 3:
                b,e,s = frames
                ase_dbtraj_list = read(self.trajPath, index=f'{b}:{e}:{s}')
            else:
                raise ValueError("Frame tupe must be of len=3 (begin, end, stide)")
        elif type(frames) == str:
            if frames == 'all':
                ase_dbtraj_list = read(self.trajPath, index=f':')
            else:
                raise ValueError("Only accepted value is 'all' to read all the frames.")
        return ase_dbtraj_list

    def _get_molCOM(self, ase_db, molID, molSym):
        # how to set their names ??
        ase_COMdb = list()
        for at in tqdm(ase_db, desc='Computing mol COM'):
            molCOM_tmp = list()
            for m in np.unique(molID):
                mol = at[molID==m] #copy by value
                mass = mol.get_masses()
                cm = np.sum(mol.positions*mass.reshape(-1,1), axis=0)/np.sum(mass)
                molCOM_tmp.append(cm)
            newmol = Atoms(positions=np.array(molCOM_tmp), pbc=True, cell=at.cell)
            newmol.arrays['molSym'] = np.array(molSym)
            ase_COMdb.append(newmol)
        return ase_COMdb

    def _shift_Znumbers(self, Zshift_tuple):
        mol_toshift, Z_toshift = Zshift_tuple
        toshift = np.max(self.infoDict['atmsZ'])
        new_Znumbers = self.infoDict['atmsZ'].copy()
        # loop inside the moltypes and select the right type
        for idx, m in enumerate(self.infoDict['molSym']):
            if m == mol_toshift:
                # create a mask
                mask = self.infoDict['molID'] == idx
                for i,ture in enumerate(mask):
                    if ture:
                        if new_Znumbers[i] in Z_toshift:
                            new_Znumbers[i] += toshift
                        else:
                            pass
        return new_Znumbers

    def save_traj(self, ase_db, file_name, frames=None, format=None):
        # -
        if format and format[0] == '.':
            format = format[1:]
        else:
            format = 'extxyz'
        # -  
        if frames:
            if type(frames) == tuple:
                add_str = '-'.join(map(str,frames))
                save_name = file_name + add_str
            elif type(frames) == int:
                add_str = 'frame'+str(frames)
                save_name = file_name + add_str
            elif type(frames) == str:
                save_name = file_name + frames
            save_name = save_name+'.'+format
        else:
            save_name = file_name+'.'+format
        print(f"Saving trajectory to file {save_name} ...")
        write(save_name, ase_db, format=format)


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

class _BaseTraj:

    DBFILE = os.path.dirname(__file__) + "/chemFormulaToName.json"

    @property
    def _project_db(self):
        with open(self.DBFILE) as jsonfile:
            project_db = json.load(jsonfile)
        return project_db

    def __init__(self, 
                 projectName, 
                 trajPath, 
                 rcut_correction=1.0):
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

# --- helper functions

# def _project_dict(project_dict_name):
#     if type(project_dict_name) == str:
#         assert '.json' in project_dict_name
#         with open(project_dict_name) as jsonfile:
#             name_dict = json.load(jsonfile)
#     elif type(project_dict_name) == dict:
#         name_dict = project_dict_name
#     return name_dict
# -------------------------------------------------- #
# ASE tools - atoms tools module
# 
#
# AUTHOR: Andrea Gardin, Ioan-Bogdan Magdau
# -> I make use of just a few of ibm tools, the
# -> complete tools list can be foun in ibm repo
# -> @ https://github.com/imagdau/Python-Atoms
# -------------------------------------------------- #

import numpy as np
from tqdm import tqdm
import ase
from ase import Atoms, neighborlist
from scipy import sparse
from typing import Union
from phdtools.computes import misc

# -------------------------------------------------- #
# --- Atom tools

# -- extract molecules features and identity into a new database (i.e. db, ase traj obj)
def extract_molInfo(db, mol_chemName_db, fct=1):
    """
    Extract all valuable information from a database (db)
    db: -> list(ase.Atoms)
    Returns a db with COM-CG trajectory
    """
    moldb = list()
    for at in tqdm(db, desc='Extracting mol information'):
        if 'molID' not in at.arrays.keys():
            get_molIDs([at], fct=fct)
        molID = at.arrays['molID']
        molCM = list()
        molSym = list()
        molQ = list()
        molF = list()
        molT = list()
        for m in np.unique(molID):
            mol = at[molID==m] #copy by value
            mass = mol.get_masses()
            cm = np.sum(mol.positions*mass.reshape(-1,1), axis=0)/np.sum(mass)
            molCM.append(cm)
            molSym.append(mol_chemName_db[mol.symbols.get_chemical_formula()])
            if 'initial_charges' in at.arrays:
                molQ.append(np.sum(mol.arrays['initial_charges']))
            if 'forces' in at.arrays:
                molF.append(np.sum(mol.arrays['forces'], axis=0))
                molT.append(np.sum(np.cross(mol.positions-cm, mol.arrays['forces'], axis=1), axis=0))
        newmol = Atoms(positions=np.array(molCM), pbc=True, cell=at.cell)
        newmol.arrays['molID'] = molID
        newmol.arrays['molSym'] = np.array(molSym)
        if molQ:
            newmol.arrays['initial_charges'] = np.array(molQ)
        if molF:
            newmol.arrays['forces'] = np.array(molF)
            newmol.arrays['torques'] = np.array(molT)
        moldb.append(newmol)
    return moldb

@misc.my_timer
def get_molIDs(db: Union[list, ase.ase.Atoms], 
               fct: Union[float, ase.ase.Atoms]) -> list:
    if type(db) == list:
        for at in db:
            # doc @ https://wiki.fysik.dtu.dk/ase/ase/neighborlist.html
            _, molID = get_connected_atoms(at, fct)
            at.arrays['molID'] = molID
    else:
        _, molID = get_connected_atoms(db, fct)
        db.arrays['molID'] = molID
    print(f"Total numner of molecule: {len(molID)}")
    return molID


def get_chemFormulas(at: ase.ase.Atoms, 
                     fct: Union[float, dict] = 1.0) -> dict:
    """Returns a dictionary with the whole molecules inside the
    system (dependent on the rcut correction, fct parameter).

    :param at: Atomic configuration in ase format.
    :type at: ase.ase.Atoms
    :param fct: Rcut correction, defaults to 1.0 (i.e., no correction)
    :type fct: Union[float, dict], optional
    :return: Whole molecules dictionary.
    :rtype: dict
    """
    _, molID = get_connected_atoms(at, fct)
    chemFormulas_list = list()
    for m in np.unique(molID):
        mol = at[molID==m]
        chemFormulas_list.append(mol.symbols.get_chemical_formula())
    chemFormulas_dict = dict()
    for i,chem in enumerate(np.unique(chemFormulas_list)):
        chemFormulas_dict[chem] = f'mol{i+1}'
    return chemFormulas_dict
    
# - computes molID for single config, not adding molID to atoms.arrays
def find_molecules(at, fct=1.0):
    # doc @ https://wiki.fysik.dtu.dk/ase/ase/neighborlist.html
    _, molID = get_connected_atoms(at, fct)
    Natoms, Nmols = np.unique(np.unique(molID, return_counts=True)[1], return_counts=True)
    return Nmols, Natoms


def get_connected_atoms(at: ase.ase.Atoms, 
                        fct: Union[float, dict]) -> tuple[int, np.ndarray]:
    cutOff = modif_natural_cutoffs(at, fct)
    nbLst = neighborlist.NeighborList(cutOff, 
                                      self_interaction=False, 
                                      bothways=True)
    nbLst.update(at)
    conMat = nbLst.get_connectivity_matrix(sparse=True)
    Nmol, molID = sparse.csgraph.connected_components(conMat)
    return Nmol, molID

# --------------------------------------------------

# -- natural cutoff modifier
def modif_natural_cutoffs(at: ase.ase.Atoms,
                          fct: Union[float, dict]):
    if type(fct) is int or type(fct) is float:
        return neighborlist.natural_cutoffs(at, mult=fct)
    elif type(fct) is dict:
        cutOff = neighborlist.natural_cutoffs(at, mult=1)
        return [ctf*fct[el] for ctf, el in zip(cutOff, at.get_chemical_symbols())]
    else:
        raise NameError('Unknown fct type '+str(type(fct)))



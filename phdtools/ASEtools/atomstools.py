# -------------------------------------------------- #
# ASE tools - atoms tools module
# 
#
# AUTHOR: Andrea Gardin, Ioan-Bogdan Magdau
# -> I made use of just a few of ibm tools;
# -> the complete tools list can be foun in 
# -> ibm repo @ https://github.com/imagdau/Python-Atoms
# -------------------------------------------------- #

import numpy as np
from tqdm import tqdm
import ase
from ase import Atoms, neighborlist
from scipy import sparse
from typing import Union, Tuple, List
from ..computes import misc

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
def get_molIDs(at: ase.ase.Atoms, 
               fct: Union[float, ase.ase.Atoms]) -> list:
    # doc @ https://wiki.fysik.dtu.dk/ase/ase/neighborlist.html
    _, molID = get_connected_atoms(at=at, 
                                   fct=fct)
    at.arrays['molID'] = molID
    return molID

@misc.my_timer
def get_molSym(at: ase.ase.Atoms,
               molIDs: np.ndarray,
               mol_name: dict) -> list:
    """Get the molecules symbols as they appear in the system xyz configuration.

    :param at: Atoms configuration.
    :type at: ase.ase.Atoms
    :param molIDs: Molecules IDs.
    :type molIDs: np.ndarray
    :param mol_name: Names association to chemical formulas.
    :type mol_name: dict
    :return: List of molecules of the atomic configuration.
    :rtype: list
    """
    molSym = list()
    for m in np.unique(molIDs):
        mol = at[molIDs == m]
        molSym.append(mol_name[mol.symbols.get_chemical_formula()])
    print(f"Total numner of molecules: {len(molSym)}")
    return molSym


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
def find_molecules(at: ase.ase.Atoms, 
                   fct: Union[float, dict]):
    """Computes the whole molecules based on the LJ cutoff values of each 
    atoms in the configuration.

    :param at: ase atom configuration.
    :type at: ase.ase.Atoms
    :param fct: scaling parameters for the LJ cutoffs.
    :type fct: Union[float, dict]
    :return: _description_
    :rtype: _type_
    """
    _, molID = get_connected_atoms(at, fct)
    Natoms, Nmols = np.unique(np.unique(molID, return_counts=True)[1], return_counts=True)
    return Nmols, Natoms


def get_connected_atoms(at: ase.ase.Atoms, 
                        fct: Union[float, dict]) -> Tuple[int, np.ndarray]:
    """Computes connected atoms based on the natural LJ cutoff range.
    Doc @ https://wiki.fysik.dtu.dk/ase/ase/neighborlist.html

    :param at: ase atom configuration.
    :type at: ase.ase.Atoms
    :param fct: scaling parameters for the LJ cutoffs.
    :type fct: Union[float, dict]
    :return: _description_
    :rtype: Tuple[int, np.ndarray]
    """
    cutOff = modif_natural_cutoffs(at, fct)
    nbLst = neighborlist.NeighborList(cutOff, 
                                      self_interaction=False, 
                                      bothways=True)
    nbLst.update(at)
    conMat = nbLst.get_connectivity_matrix(sparse=True)
    Nmol, molID = sparse.csgraph.connected_components(conMat)
    return Nmol, molID


def modif_natural_cutoffs(at: ase.ase.Atoms,
                          fct: Union[float, dict]) -> dict:
    """Modifies the natural cutoff of the LJ interactions.

    :param at: ase atom configuration.
    :type at: ase.ase.Atoms
    :param fct: newly defined scaling parameters for the LJ cutoffs.
    :type fct: Union[float, dict]
    :raises NameError: only accepts int and dictionary values.
    :return: newly defined dictionary containing the LJ scaling values for each atoms.
    :rtype: dict
    """
    if type(fct) is int or type(fct) is float:
        return neighborlist.natural_cutoffs(at, mult=fct)
    elif type(fct) is dict:
        cutOff = neighborlist.natural_cutoffs(at, mult=1)
        return [ctf*fct[el] for ctf, el in zip(cutOff, at.get_chemical_symbols())]
    else:
        raise NameError('Unknown fct type '+str(type(fct)))


def ZnumberShift(Znumbers: np.ndarray, 
                 molSymbols: list,
                 molIDs: list,
                 to_shift: Tuple[str, list]) -> np.ndarray:
    
    mol_to_shift, Z_to_shift = to_shift
    dummy = np.max(Znumbers)
    shifted_Znumbers = Znumbers.copy()
    #
    for idx, mol in enumerate(molSymbols):
        if mol == mol_to_shift:
            mask = molIDs == idx
            #
            for i, tgt in enumerate(mask):
                if tgt:
                    if shifted_Znumbers[i] in Z_to_shift:
                        shifted_Znumbers[i] += dummy
                    else:
                        pass
    return shifted_Znumbers


def center_of_mass(ase_db: List[ase.ase.Atoms],
                   molSymbols: list,
                   molIDs: list) -> List[ase.ase.Atoms]:
    """Computes the COM of a given ase atoms databas of frames.

    :param ase_db: ase atoms database.
    :type ase_db: List[ase.ase.Atoms]
    :param molSymbols: list of molecule-wise symbols as they appear in the frame configuration.
    :type molSymbols: list
    :param molIDs: list of molecule-wise ids as they appear in the frame configuration.
    :type molIDs: list
    :return: ase atoms database containitng the COM position.
    :rtype: List[ase.ase.Atoms]
    """
    ase_db_com_list = list()
    for at in tqdm(ase_db, desc='Computing COM:'):
        mol_com_tmp = list()
        for m in np.unique(molIDs):
            mol = at[molIDs==m] #copy by value
            mass = mol.get_masses()
            cm = np.sum(mol.positions*mass.reshape(-1,1), axis=0)/np.sum(mass)
            mol_com_tmp.append(cm)
        new_com_at = Atoms(positions=np.array(mol_com_tmp), pbc=True, cell=at.cell)
        new_com_at.arrays['molSym'] = np.array(molSymbols)
        ase_db_com_list.append(new_com_at)
    return ase_db_com_list
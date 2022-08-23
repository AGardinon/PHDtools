#!
from textwrap import indent
import numpy as np
from ase.io import read
from phdtools.ASEtools import atomstools



trajpath = 'data/traj_2.1_0-100-1.xyz'
rcut_correction = {'H':1,'C':1,'O':1,'Li':0.1,'P':1,'F':1}

asetraj = read(trajpath, index='0:1')

print(type(asetraj) == list)

dummy = atomstools.extract_molInfo(asetraj[:1], 
mol_chemName_db={"C3H4O3": "EC", "C4H8O3": "EMC", "Li": "Li", "F6P": "PF6"},
fct=rcut_correction
)

print(dummy[0].arrays.keys())
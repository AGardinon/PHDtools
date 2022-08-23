#!
from textwrap import indent
import numpy as np
from ase.io import read
from phdtools.ASEtools import atomstools

trajpath = 'data/traj_2.1_0-100-1.xyz'
rcut_correction = {'H':1,'C':1,'O':1,'Li':0.1,'P':1,'F':1}

asetraj = read(trajpath, index='0:1')

print(type(asetraj) == list)

atomstools.modif_natural_cutoffs(asetraj[0], rcut_correction)
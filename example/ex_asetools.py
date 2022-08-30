# -------------------------------------------------- #
# 
# Example on how to use ASEtools
# 
# -------------------------------------------------- #

import numpy as np
from phdtools.ASEtools import asetools

# 
input_dict = dict(
    projectName = "LiIonSolutions",
    trajPath = "data/traj_2.1_0-100-1.xyz",
    rcut_correction = {'H':1,'C':1,'O':1,'Li':0.1,'P':1,'F':1}
)

# init the traj-tool object
ase_uni = asetools.Universe(**input_dict)
print(ase_uni.mol_chemName_db)
print(ase_uni.infoDict.keys())

# reading traj and saving traj
traj_output = './example_outputs/traj_example_'
ex_frames = (0,10,2)

# atomistic traj with Z number shift
ase_traj_db = ase_uni.read(frames=ex_frames, 
                           Zshift_tuple=('EC', [6, 8]),
                           save_file='example_outputs/test_output_1_')

# molecules COM traj
ase_trajCOM_db = ase_uni.read_COM(frames=ex_frames, 
                                  save_file='example_outputs/test_output_2_')

print(np.shape(ase_trajCOM_db))
# problema !!!

unwrap_traj_dict = ase_uni.MolUnwrapper(ase_mol_db=ase_trajCOM_db, 
                                        mol_species=['Li', 'PF6'], 
                                        save_file='example_outputs/')

print(unwrap_traj_dict.keys())
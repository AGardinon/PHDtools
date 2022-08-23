#!

from phdtools.ASEtools import asetools

input_dict = dict(
    projectName = "LiIonSolutions",
    trajPath = "data/traj_2.1_0-100-1.xyz",
    rcut_correction = {'H':1,'C':1,'O':1,'Li':0.1,'P':1,'F':1}
)

aseTraj = asetools.Traj(**input_dict)
print(aseTraj.mol_chemName_db)

print(aseTraj.infoDict.keys())
print(aseTraj.infoDict['molID'])

traj_db = aseTraj.read(frames=(0,10,1))
assert len(traj_db) == 10

traj_COM_db = aseTraj.read_COM(frames=(0,10,1))
assert len(traj_COM_db) == 10
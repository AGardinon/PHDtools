#!

from phdtools.ASEtools import asetools

# test input
projname = 'NewProject23Aug2022'
trajpath = 'data/traj_2.1_0-100-1.xyz'
rcut_correction = {'H':1,'C':1,'O':1,'Li':0.1,'P':1,'F':1}

asetraj = asetools._BaseTraj(projectName=projname,
                             trajPath=trajpath,
                             rcut_correction=rcut_correction)

#print("prima",hasattr(asetraj, 'mol_chemName_db'))

#asetraj.addNewProject(save_to_default=True)

#print("dopo",hasattr(asetraj, 'mol_chemName_db'))
#print(asetraj.mol_chemName_db)

asetraj.deleteProject(project_name="NewProject23Aug2022")
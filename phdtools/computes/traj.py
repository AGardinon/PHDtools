# -------------------------------------------------- #
# Computes - trajectory module
#
#
# AUTHOR: Andrea Gardin
# -------------------------------------------------- #

import numpy as np
from tqdm import tqdm

# -------------------------------------------------- #
# --- Coordinate unwrapper

class XYZunwrapper:
    """Class to Unwrap trajectories."""
    METHODS = ("heuristic", "displacement", "hybrid")

    def __init__(self, xyz, box):
        self.xyz = xyz
        self.box = box

    def __call__(self):
        methods_function_dict = dict(
            heuristic = heuristic_unwrapping,
            displacement = heuristic_unwrapping,
            hybrid = hybrid_unwrapping,
        )
        self._methods_functions_dict = methods_function_dict

    def run(self, method):
        return self._run(method)

    def _run(self, method):
        if not method in XYZunwrapper.METHODS:
            raise ValueError(f'Method not recognised, available: {XYZunwrapper.METHODS}')
        print(f'Method: {method}')
        unwrappedXYZ = self._methods_functions_dict[method](w=self.xyz, 
                                                            box=self.box)
        return unwrappedXYZ

# -- Helper functions
def heuristic_unwrapping(w,box):
    # init the coordinates
    u = np.empty((len(w),3))
    # set the same starting point
    u[0] = w[0]
    # computing the diff of the w traj
    difw = np.diff(w, axis=0)
    # ---
    # Eq. 1 Heuristic method
    for i,dw in tqdm(enumerate(difw), desc='Unwrapping'):
        u[i+1] = w[i+1]-np.floor((w[i+1]-u[i])/box[i+1]+0.5)*box[i+1]
    return u

def displacement_unwrapping(w,box):
    # init the coordinates
    u = np.empty((len(w),3))
    # set the same starting point
    u[0] = w[0]
    # computing the diff of the w traj
    difw = np.diff(w, axis=0)
    # ---
    # Eq 2. Displacemnet method
    for i,dw in tqdm(enumerate(difw), desc='Unwrapping'):
        u[i+1] = u[i]+(w[i+1]-w[i])-np.floor((w[i+1]-w[i])/box[i+1]+0.5)*box[i+1]
    return u
        
def hybrid_unwrapping(w,box):
    # init the coordinates
    u = np.empty((len(w),3))
    # set the same starting point
    u[0] = w[0]
    # computing the diff of the w traj
    difw = np.diff(w, axis=0)
    # ---
    # Eq. 12 Hybrid mehtod
    for i,dw in tqdm(enumerate(difw), desc='Unwrapping'):
        u[i+1] = u[i]+(w[i+1]-w[i])-np.floor((w[i+1]-w[i])/box[i+1]+0.5)*box[i+1]-np.floor((w[i]-u[i])/box[i]+0.5)*(box[i+1]-box[i])
    return u

# -------------------------------------------------- #
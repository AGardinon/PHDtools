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

    @property
    def _methods_dict(self):
        methods_function_dict = dict(
            heuristic = heuristic_unwrapping,
            displacement = heuristic_unwrapping,
            hybrid = hybrid_unwrapping,
        )
        return methods_function_dict

    def __init__(self, method):
        self.method = method
        self.methdos_list = list(self._methods_dict.keys())
        if self.method not in self.methdos_list:
            raise ValueError(f"Chosen method not available, choose from {self.methdos_list}")

    def run(self, xyz, box, verbose=False):
        if verbose:
            print(f"Unwrapping the trajectory with {self.method} method.")
        return self._run(xyz, box)

    def _run(self, xyz, box):
        return self._methods_dict[self.method](w=xyz, box=box)

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
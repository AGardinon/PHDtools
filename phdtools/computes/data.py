# -------------------------------------------------- #
# Computes - sampling stuffs
#
#
# AUTHOR: Andrea Gardin
# -------------------------------------------------- #

import numpy as np
import random
from typing import Union

def random_shuffle(X: np.ndarray, 
                   Y: np.ndarray=None, 
                   n: int=None):
    """Random shuffler for a dataset.
    Include also the possibility to shuffle a
    properties array in the same way.
    """
    l = np.arange(X.shape[0])
    random.shuffle(l)
    if Y is None:
        return X[l[:n], :]
    elif Y is None and n is None:
        return X[l, :]
    elif n is None:
        return X[l, :], Y[l]
    else:
        return X[l[:n], :], Y[l[:n]]


def FPS(X: np.ndarray,
        n: int=-1, 
        ndx: int=None, 
        retDist: bool=False):
    """Does Farthest Point Selection on a set of points X
    X is in the form of {X_i} with X_i(x_0,x_1,...,x_N)
    where N are the features or dimensions and i are the
    data sample size.
    """
    N = X.shape[0]
    # if no n points are selected it takes all of them
    if n <= 0:
        n = N

    # init the arrays for the ndxs and distances
    fps_ndxs = np.zeros(n, dtype=np.int)
    D = np.zeros(n)
    if ndx is None:
        ndx = np.random.randint(0, N)
    fps_ndxs[0] = ndx

    # compute the distance from selected point
    dist1 = np.linalg.norm(X - X[ndx], axis=1)

    # loop over the remaining points
    for i in range(1, n):
        # get and store the index for the max dist from the point chosen
        fps_ndxs[i] = np.argmax(dist1)
        D[i - 1] = np.amax(dist1)

        # compute the dists from the newly selected point
        dist2 = np.linalg.norm(X - X[fps_ndxs[i]], axis=1)
        # takes the min from the two arrays dist1 2
        dist1 = np.minimum(dist1, dist2)

        # little stopping condition
        if np.abs(dist1).max() == 0.0:
            print(f"Only {i} iteration possible")
            return fps_ndxs[:i], D[:i]

    if retDist:
        return X[fps_ndxs], fps_ndxs, D
    else:
        return X[fps_ndxs], fps_ndxs

# there is a bug, so it has been removed
# def normalizeData(X: np.ndarray, 
#                   mean: Union[int, np.ndarray, ]=None, 
#                   std: np.ndarray=None):
#     """Normalize a dataset with default mean=0 and std=1
#     """
#     if mean is None:
#         mean = np.mean(X, axis=0)
#     Xcenter = X - mean

#     if std is None:
#         std = np.linalg.norm(X) / np.sqrt(len(Xcenter))
#     Xnorm = Xcenter / std

#     return Xnorm
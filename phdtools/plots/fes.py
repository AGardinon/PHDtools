# -------------------------------------------------- #
# Plot tools -General tools
#
#
# AUTHOR: Andrea Gardin
# -------------------------------------------------- #

import numpy as np
from .general import *


def plot_fes_1d(fes: np.ndarray,
                grid: np.ndarray,
                axes,
                kwargs: dict =None):
    """Plot 1D fes.
    """
    try:
        axes.plot(grid[:-1], fes, **kwargs)
    except:
        axes.plot(grid[:-1], fes)
    pass


def plot_fes_2d(fes: np.ndarray,
                grid: np.ndarray,
                levels: int,
                axes,
                cont_kwargs: dict =None,
                cmap = 'coolwarm_r'):
    """Plot 2D fes.
    """
    default_cont_kwargs = dict(
        colors = 'k',
        linewidths = .5,
        alpha= 0.0
    )

    X, Y = grid

    # contour object
    try:
        contour = axes.contour(X, Y, fes, 
                               levels,
                               **cont_kwargs,
                               zorder=2)
    except:
        contour = axes.contour(X, Y, fes, 
                               levels,
                               **default_cont_kwargs, 
                               zorder=2)

    # surface object
    surface = axes.contourf(X, Y, fes, 
                            levels,
                            cmap=cmap, 
                            zorder=1)
    
    return surface, contour

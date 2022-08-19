# -------------------------------------------------- #
# Plot tools - FES
#
#
# AUTHOR: Andrea Gardin
# -------------------------------------------------- #

import numpy as np
import matplotlib.pyplot as plt

# -------------------------------------------------- #
# --- FES

class _BaseFES:
    """Base FES class to init the computation"""

    kbUNITS = dict(
        kJ = 0.00831446261815324,
        kcal = 0.00198720425864083,
        unit =  1.0
    )

    def __init__(self, units, temp):
        self.units = units
        self.temp = temp

    @property
    def kb(self):
        return _BaseFES.kbUNITS[self.units]

    @property
    def kbT(self):
        if self.units == 'unit':
            return self.kb
        else:
            return self.kb * self.temp


class FES(_BaseFES):
    """Class to compute and plot the pseudo-FES"""
    def __init__(self, units, temp):
        super().__init__(units, temp)

    def fit_plot(self):
        pass

    def fit(self, X, Y, bins,
            range=None, fillEmpty=True):
        return self._fit(X, Y, bins, range=None, fillEmpty=True)

    def plot(self):
        print('plotting time')
        pass

    def _fit(self, X, Y, bins, range=None, fillEmpty=True):
        # <-> Double variables, 2D fes <->
        if isinstance(Y, (list, np.ndarray)):
            H, xedges, yedges = np.histogram2d(X, Y, bins=bins, 
                                               range=range, density=True)
            Z = self._histo_to_fes(H, fillEmpty)
            XX, YY = self._2DmeshGrid(xedges, yedges, bins)
            compute_dict = dict(fes = Z.T, 
                                mesh_grid = (XX, YY))
            self.fesDim = '2D'
            return compute_dict

        # <-> Single variable, 1D fes <->
        elif not Y:
            # 1. compute Histograms
            H, edges = np.histogram(X, bins=bins, range=range, density=True)
            # 2. compute 1D pseudofes
            Z = self._histo_to_fes(H, fillEmpty)
            # 3. compute minvalue
            edges_ = edges[:-1]
            minVal = edges_[Z == np.min(Z)]
            compute_dict = dict(fes = Z.T, 
                                edges = edges, 
                                min_value = minVal)
            self.fesDim = '1D'
            return compute_dict

    def _histo_to_fes(self, H, fillEmpty):
        # inversion
        Z = -1 * self.kbT * np.log(H)
        # min to zero
        Z = Z - np.min(Z)
        if fillEmpty:
            max_ = np.max(Z[Z != np.inf])
            Z[Z == np.inf] = max_
        return Z

    def _2DmeshGrid(self, x, y, bins):
        Xmin = x.min()
        Xmax = x.max()
        Ymin = y.min()
        Ymax = y.max()
        xx = np.arange(Xmin, Xmax, ((Xmax - Xmin) / bins))
        yy = np.arange(Ymin, Ymax, ((Ymax - Ymin) / bins))
        XX, YY = np.meshgrid(xx, yy)
        return XX, YY
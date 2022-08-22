# -------------------------------------------------- #
# Plot tools - FES
#
#
# AUTHOR: Andrea Gardin
# -------------------------------------------------- #

import numpy as np
import phdtools.plots as phdplot

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
    """Class to compute and plot the pseudo-FES."""
    def __init__(self, units, temp):
        super().__init__(units, temp)

    def fit_plot(self, X, Y, bins, range=None, levels=None):
        """
        Returns a Dictionary with the fit data.
        According if the fit is 1D or 2D it changes.
        Quickly plot the data as well.
        """
        if not levels:
            levels = 10
        compute_dict = self._fit(X, Y, bins, range=range)
        self._plot(levels=levels,**compute_dict)
        return compute_dict

    def fit(self, X, Y, bins,
            range=None, fillEmpty=True):
        """
        Returns a Dictionary with the fit data.
        According if the fit is 1D or 2D it changes.
        """
        return self._fit(X, Y, bins, range=None, fillEmpty=True)

    def plot(self):
        print("Plotting time")
        print("Not yet implemented.")
        # probably best to put the plotting
        # function in plt tools as stand alone
        # in addition to this.
        pass

    def _fit(self, X, Y, bins, range=None, fillEmpty=True):
        # -> Double variables, 2D fes
        if isinstance(Y, (list, np.ndarray)):
            H, xedges, yedges = np.histogram2d(X, Y, bins=bins, 
                                               range=range, density=True)
            Z = self._histo_to_fes(H, fillEmpty)
            XX, YY = self._2DmeshGrid(xedges, yedges, bins)
            compute_dict = dict(fes = Z.T, 
                                grid = (XX, YY))
            self._fesDim = '2D'
            return compute_dict

        # -> Single variable, 1D fes
        elif not Y:
            # 1. compute Histograms
            H, edges = np.histogram(X, bins=bins, range=range, density=True)
            # 2. compute 1D pseudofes
            Z = self._histo_to_fes(H, fillEmpty)
            # 3. compute minvalue
            edges_ = edges[:-1]
            minVal = edges_[Z == np.min(Z)]
            compute_dict = dict(fes = Z.T, 
                                grid = edges, 
                                min_value = minVal)
            self._fesDim = '1D'
            return compute_dict

    def _plot(self, grid, fes, levels,
              #fesArgs, # -> needed?
              figure=None, axes=None, 
              ghost=False,
              contlabels=True,
              cbar=True,
              cbar_label=None):
        # Def fig and axes if not defined
        if not figure and not axes:
            figure, axes = phdplot.get_axes(1,1)

        # Plot either 1D or 2D FES
        # -> 1D
        if self._fesDim == '1D':
            print(f"Plotting {self._fesDim} FES.")
            axes.plot(grid[:-1],fes)

        # -> 2D
        elif self._fesDim == '2D':
            print(f"Plotting {self._fesDim} FES.")
            X,Y = grid
            cont = axes.contour(X, Y, fes, levels,
                                colors='k', 
                                linewidths=0.5, 
                                zorder=2)
            if not ghost:
                surf = axes.contourf(X, Y, fes, levels,
                                     cmap='coolwarm_r', 
                                     zorder=1)
                if cbar:
                    cbar = figure.colorbar(surf,ax=axes)
                    if cbar_label:
                        cbar.set_label(cbar_label)
            if contlabels:
                axes.clabel(cont, inline=True, 
                            colors='k', fontsize=8, 
                            fmt='%1.1f', zorder=3)
        
        figure.savefig(f"{self._fesDim}_pseudo_fes.png")
        


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
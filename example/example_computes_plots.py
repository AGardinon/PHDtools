# -------------------------------------------------- #
# 
# Example on how to use Plots tools
# 
# -------------------------------------------------- #

import numpy as np
from phdtools import plots
from phdtools import computes

# - gen dummy dataset
from sklearn import datasets

dummy_dict = dict(
    n_samples = 10000,
    factor = 0.4,
    noise = 0.1,
    random_state = 73,
)
X, y = datasets.make_circles(**dummy_dict)
print(X.shape, y.shape)

# - setting the image
fig_settings = dict(
    P = 3,
    max_col = 2,
    fig_frame = (5,4),
    res = 200,
)

fig, ax = plots.get_axes(**fig_settings)

ax[0].scatter(*X.T, c=y, s=5)

# - computes and draw pseudo-fes
FES = computes.fes.FES(units='kcal', temp=300)
print(FES.kbT)

fes_compute = dict(
    X = X[:,0],
    Y = X[:,1],
    bins = 42,
    range=((-1.3,1.3),(-1.3,1.3))
)
fes_plot = dict(
    figure=fig, axes=ax[1],
    levels = 4,
    contlabels = False,
    cbar_label = r'FES [Kcal mol$^{-1}$]'
)

fes_output = FES.fit_plot(**fes_compute, plotArgs=fes_plot)


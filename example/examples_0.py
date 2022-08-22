# -
#
# Examples
#
# -

import numpy as np
import matplotlib.pyplot as plt
from sklearn import datasets

# phd tools
import phdtools as phd

print("# -\t Generating Dummy dataset")
dummy_dict = dict(
    n_samples = 10000,
    factor = 0.1,
    noise = 0.5,
    random_state = 73,
)
X, y = datasets.make_circles(**dummy_dict)
print(X.shape)

print("# -\t Plots examples")
fig_dict = dict(
    P = 3,
    max_col = 2,
    fig_frame = (5,4),
    res = 100,
)
print(fig_dict)

fig, ax = phd.plots.get_axes(**fig_dict)
ax[0].scatter(*X.T)
ax[1].scatter(*X.T)
phd.plots.axarrows(axes=ax[1], figure=fig)

fes_dict = dict(
    units = 'kcal',
    temp = 300,
)
fes = phd.computes.FES(**fes_dict)
print(fes.kbT)

fes_compute = dict(
    X = X[:,0],
    Y = X[:,1],
    bins = 42,
)

compute_output = fes.fit(**fes_compute)
print(compute_output)
print(len(compute_output['grid']))
#print(fes.compute_dict.keys())

compute_plot_output = fes.fit_plot(**fes_compute)
print(compute_plot_output)
# -------------------------------------------------- #
# Plot tools
#
#
# AUTHOR: Andrea Gardin
# -------------------------------------------------- #

import matplotlib.pyplot as plt

# -------------------------------------------------- #
# --- FIG & AXES

def get_axes(L: int, max_col: int =2, fig_frame: tuple =(3.3,3.), res: int =200):
    """
    Define Fig and Axes objects.

    """
    # cols and rows definitions
    cols = L if L <= max_col else max_col
    rows = int(L / max_col) + int(L % max_col != 0)

    fig, axes = plt.subplots(rows,
                             cols,
                             figsize=(cols * fig_frame[0], rows * fig_frame[1]),
                             dpi=res)
    # beauty
    if L > 1:
        axes = axes.flatten()
        for i in range(L, max_col*rows):
            remove_frame(axes[i])
    elif L == 1:
        pass
    
    return fig, axes


def remove_frame(axes):
    for side in ['bottom', 'right', 'top', 'left']:
        axes.spines[side].set_visible(False)
    axes.set_yticks([])
    axes.set_xticks([])
    axes.xaxis.set_ticks_position('none')
    axes.yaxis.set_ticks_position('none')
    pass

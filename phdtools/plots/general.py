# -------------------------------------------------- #
# Plot tools -General tools
#
#
# AUTHOR: Andrea Gardin
# -------------------------------------------------- #

import matplotlib.pyplot as plt

# -------------------------------------------------- #
# --- FIG & AXES

# -- Set Fig and Axes
def get_axes(plots: int, 
             max_col: int =2, 
             fig_frame: tuple =(3.3,3.), 
             res: int =200):
    """Define Fig and Axes objects.

    :param plots: number of plots frames in a fig.
    :type plots: int
    :param max_col: number of columns to arrange the frames, defaults to 2
    :type max_col: int, optional
    :param fig_frame: frame size, defaults to (3.3,3.)
    :type fig_frame: tuple, optional
    :param res: resolution, defaults to 200
    :type res: int, optional
    :return: fig and axes object from matplotlib.
    :rtype: _type_
    """
    # cols and rows definitions
    cols = plots if plots <= max_col else max_col
    rows = int(plots / max_col) + int(plots % max_col != 0)

    fig, axes = plt.subplots(rows,
                             cols,
                             figsize=(cols * fig_frame[0], rows * fig_frame[1]),
                             dpi=res)
    # beauty
    if plots > 1:
        axes = axes.flatten()
        for i in range(plots, max_col*rows):
            remove_frame(axes[i])
    elif plots == 1:
        pass
    
    return fig, axes


# -- Arrowed frame axes
def axarrows(axes, figure):
    # - 1
    xmin, xmax = axes.get_xlim()
    ymin, ymax = axes.get_ylim()
    remove_frame(axes)
    # - 2
    # get width and height of axes object to compute
    # matching arrowhead length and width
    dps = figure.dpi_scale_trans.inverted()
    bbox = axes.get_window_extent().transformed(dps)
    width, height = bbox.width, bbox.height
    # manual arrowhead width and length
    hw = 1./20.*(ymax-ymin)
    hl = 1./20.*(xmax-xmin)
    lw = 1. # axis line width
    ohg = 0.3 # arrow overhang
    # compute matching arrowhead length and width
    yhw = hw/(ymax-ymin)*(xmax-xmin)* height/width
    yhl = hl/(xmax-xmin)*(ymax-ymin)* width/height
    # - 3
    # draw x and y axis
    axes.arrow(xmin, ymin, xmax-xmin, 0., fc='k', ec='k', lw=lw,
             head_width=hw, head_length=hl, overhang=ohg,
             length_includes_head=True, clip_on=False)
    axes.arrow(xmin, ymin, 0., ymax-ymin, fc='k', ec='k', lw=lw,
             head_width=yhw, head_length=yhl, overhang=ohg,
             length_includes_head=True, clip_on=False)

    pass


# -- Remove frame
def remove_frame(axes):
    for side in ['bottom', 'right', 'top', 'left']:
        axes.spines[side].set_visible(False)
    axes.set_yticks([])
    axes.set_xticks([])
    axes.xaxis.set_ticks_position('none')
    axes.yaxis.set_ticks_position('none')
    pass

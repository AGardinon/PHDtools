# -------------------------------------------------- #
# MSM tools - markov state model tools
#
# AUTHOR: Andrea Gardin 
# -------------------------------------------------- #

import numpy as np

# -------------------------------------------------- #
# --- Time series util
from deeptime.util.data import TrajectoryDataset, TrajectoriesDataset

# - init timeseries
def build_timeseries(trajfile, lagtime=1):
    # if trajfile is a file string location
    if isinstance(trajfile, str):
        data = np.load(trajfile)
    # if trajfile is a data tensor
    else:
        data = trajfile
    # for multiparticles
    # the expected shape of data is "static"
    # (N_frames x N_at x N_features)
    if len(np.shape(data)) > 2:
        print(f"Multi particles data - {np.shape(data)}")
        traj = [np.vstack([ts[at] for ts in data]) for at in np.arange(0,data.shape[1])]
        # create time series
        traj_list = [TrajectoryDataset(lagtime=lagtime,
                                       trajectory=ts.astype(np.float64))
                     for ts in traj]
        dataset = TrajectoriesDataset(traj_list)
    else:
        dataset = TrajectoryDataset(lagtime=lagtime, trajectory=data.astype(np.float64))
    return dataset
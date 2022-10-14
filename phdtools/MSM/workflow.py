# -------------------------------------------------- #
# MSM tools - markov state model tools
#
# AUTHOR: Andrea Gardin 
# -------------------------------------------------- #
#
# Wokkiflow:
# - timeseries
# - tica
# - discretization
# - MSM

import json
import numpy as np
import phdtools.MSM.msmtools as msmtools

class baseMSM:

    def __init__(self, input_dict):
        if type(input_dict) == str:
            with open(input_dict) as jsonfile:
                workflowInput = json.load(jsonfile)
        elif type(input_dict) == dict:
            workflowInput = input_dict
        else:
            raise ValueError("Input must be a .json file or a py dictonary")

        self.workflowInput = workflowInput
        self.data = np.load(self.workflowInput["data"])


    def timeSeries(self, lag):
        return msmtools.build_timeseries(trajfile=self.data, lagtime=lag)


class myMSM(baseMSM):

    def __init__(self, input_dict):
        super().__init__(input_dict)


    def runWorflow(self):
        print("Workflow as input:\n")
        for k in self.workflowInput.key:
            print(k)
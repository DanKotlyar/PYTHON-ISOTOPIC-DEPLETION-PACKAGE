# -*- coding: utf-8 -*-
"""example_depletion

A case that shows how the depletion is carried out.

Created on Mon Oct 11 21:30:00 2021 @author: Dan Kotlyar
Last updated on Sun Oct 17 18:15:00 2021 @author: Dan Kotlyar

"""

import numpy as np

from pyIsoDep.functions.maindepletionsolver import MainDepletion
from pyIsoDep.functions.generatedata import TransmutationData
from pyIsoDep.functions.postprocessresults import Results
# import pre-generated data
from pyIsoDep.tests.pregenerated_decay_isotopes import ID, N0,\
    timepoints, timeUnits, vol

from datetime import datetime

start_time = datetime.now()

# -----------------------------------------------------------------------------
#                            DATA GENERATION
# -----------------------------------------------------------------------------
# Reset the data container
data = TransmutationData(libraryFlag=True, wgtFY=1.0)

# -----------------------------------------------------------------------------
#                            DECAY CALCULATIONS
# -----------------------------------------------------------------------------
dep = MainDepletion(0.0, data)
# define metadata (steps, flux, and so on)
dep.SetDepScenario(timeUnits=timeUnits, timepoints=timepoints)
# set initial composition
dep.SetInitialComposition(ID, N0, vol=vol)
# solve the Bateman equations
dep.SolveDecay(method="cram")
# Post depletion analysis
dep.DecayHeat()
dep.Radiotoxicity()
dep.Activity()
dep.Mass()

time_elapsed = datetime.now() - start_time
print('Time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))


# -----------------------------------------------------------------------------
#                  POST-PROCESS RESULTS
# -----------------------------------------------------------------------------
res = Results(dep)
res.plot("totalQt", timeUnits="hours", norm=1E+6, ylabel="Decay heat, MW",
         pltType="loglog", markers="--r^")
res.plot("Nt", timeUnits="hours", isotopes=[390900, 942380],
         ylabel="Atomic density",  pltType="loglog")

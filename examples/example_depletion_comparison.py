# -*- coding: utf-8 -*-
"""example_depletion

Example on how to execute depletion calculations

Created on Sun Oct 11 07:30:00 2021 @author: Dan Kotlyar
Last updated on Sat Oct 16 15:30:00 2021 @author: Dan Kotlyar

"""

import numpy as np

from pyIsoDep.functions.maindepletionsolver import MainDepletion
from pyIsoDep.functions.generatedata import TransmutationData
from pyIsoDep.functions.postprocessresults import Results

from pyIsoDep.tests.pregenerated_xs import ID, N0, sig_c,\
    sig_c2m, sig_n2n, sig_n3n, sig_f

import timeit


# -----------------------------------------------------------------------------
#                            DATA GENERATION
# -----------------------------------------------------------------------------
tic = timeit.timeit()
# Reset the data container
data = TransmutationData(libraryFlag=True, wgtFY=1.0)
# Feed cross sections into the container
data.ReadData(ID, sig_f=sig_f, sig_c=sig_c, sig_c2m=sig_c2m,
              sig_n2n=sig_n2n, sig_n3n=sig_n3n, flagBarns=False)
toc = timeit.timeit()
print("Time elapsed for reading data = {} seconds ".format(toc-tic))

# -----------------------------------------------------------------------------
#                            DEPLETION
# -----------------------------------------------------------------------------
depCram = MainDepletion(0.0, data)
nsteps = 10
timeDays = 25*np.ones(nsteps)
power = 348E+6*np.ones(nsteps)
volume = 332097.750  # volume in cm**3
# define metadata (steps, flux, and so on)
depCram.SetDepScenario(power=power, timeUnits="days", timesteps=timeDays)
# set initial composition
depCram.SetInitialComposition(ID, N0, vol=volume)
# solve the Bateman equations
depCram.SolveDepletion(method="cram")
# Post depletion analysis
depCram.DecayHeat()
depCram.Radiotoxicity()
depCram.Activity()
depCram.Mass()


depExpm = MainDepletion(0.0, data)
nsteps = 10
timeDays = 25*np.ones(nsteps)
power = 348E+6*np.ones(nsteps)
volume = 332097.750  # volume in cm**3
# define metadata (steps, flux, and so on)
depCram.SetDepScenario(power=power, timeUnits="days", timesteps=timeDays)
# set initial composition
depCram.SetInitialComposition(ID, N0, vol=volume)
# solve the Bateman equations
depCram.SolveDepletion(method="expm")
# Post depletion analysis
depCram.DecayHeat()
depCram.Radiotoxicity()
depCram.Activity()
depCram.Mass()


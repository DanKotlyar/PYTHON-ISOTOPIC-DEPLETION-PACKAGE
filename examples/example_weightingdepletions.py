# -*- coding: utf-8 -*-
"""example_weighting

Example on how to weight diffent depletion objects

Created on Wed Feb 09 15:30:00 2022 @author: Dan Kotlyar
Last updated on Wed Feb 09 15:30:00 2022 @author: Dan Kotlyar

"""

import numpy as np

from pyIsoDep.functions.maindepletionsolver import MainDepletion
from pyIsoDep.functions.generatedata import TransmutationData
from pyIsoDep.functions.weightdepletionobjects import WeightDepObjects
from pyIsoDep.functions.postprocessresults import Results

from pyIsoDep.tests.pregenerated_xs import ID, N0, sig_c,\
    sig_c2m, sig_n2n, sig_n3n, sig_f


# -----------------------------------------------------------------------------
#                            DATA GENERATION
# -----------------------------------------------------------------------------
# Reset the data container
data = TransmutationData(libraryFlag=True, wgtFY=1.0)
# Feed cross sections into the container
data.ReadData(ID, sig_f=sig_f, sig_c=sig_c, sig_c2m=sig_c2m,
              sig_n2n=sig_n2n, sig_n3n=sig_n3n, flagBarns=False)

# -----------------------------------------------------------------------------
#                            DEPLETION-1
# -----------------------------------------------------------------------------
dep1 = MainDepletion(0.0, data)
nsteps = 10
timeDays = 25*np.ones(nsteps)
power = 348E+6*np.ones(nsteps)
volume = 332097.750  # volume in cm**3
# define metadata (steps, flux, and so on)
dep1.SetDepScenario(power=power, timeUnits="days", timesteps=timeDays)
# set initial composition
dep1.SetInitialComposition(ID, N0, vol=volume)
# solve the Bateman equations
dep1.SolveDepletion(method="cram")
# Post depletion analysis
dep1.DecayHeat()
dep1.Radiotoxicity()
dep1.Activity()
dep1.Mass()
dep1.Reactivity()

# -----------------------------------------------------------------------------
#                            DEPLETION-2 (different volume)
# -----------------------------------------------------------------------------
dep2 = MainDepletion(0.0, data)
nsteps = 10
timeDays = 25*np.ones(nsteps)
power = 348E+6*np.ones(nsteps)
volume = 55555.750  # volume in cm**3
# define metadata (steps, flux, and so on)
dep2.SetDepScenario(power=power, timeUnits="days", timesteps=timeDays)
# set initial composition
dep2.SetInitialComposition(ID, N0, vol=volume)
# solve the Bateman equations
dep2.NoDepletion()
# Post depletion analysis
dep2.DecayHeat()
dep2.Radiotoxicity()
dep2.Activity()
dep2.Mass()
dep2.Reactivity()


# -----------------------------------------------------------------------------
#                            DEPLETION-weighted
# -----------------------------------------------------------------------------

wgtDep = WeightDepObjects(dep1, dep2)

# -----------------------------------------------------------------------------
#                  POST-PROCESS RESULTS
# -----------------------------------------------------------------------------
res = Results(wgtDep)
res.plot("Nt", timeUnits="hours", isotopes=[922350],
         ylabel="Atomic density")
res.getvalues("totalQt")
res.plot("totalQt", timeUnits="days", norm=1E+6, ylabel="Decay heat, MW")
res.plot("Nt", timeUnits="hours", isotopes=[531350, 541350],
         ylabel="Atomic density")
res.plot("Qt", timeUnits="hours", isotopes=[531350],
         ylabel="Decay heat, Watts")

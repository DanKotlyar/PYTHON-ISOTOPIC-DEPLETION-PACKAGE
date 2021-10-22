# -*- coding: utf-8 -*-
"""example_depletion

A case that shows how the depletion is carried out.

Created on Mon Oct 11 21:30:00 2021 @author: Dan Kotlyar
Last updated on Mon Oct 11 21:45:00 2021 @author: Dan Kotlyar

"""

import numpy as np

from pyIsoDep.functions.maindepletionsolver import MainDepletion
from pyIsoDep.functions.generatedata import TransmutationData
from pyIsoDep.functions.read_csv import ReadCsv
from pyIsoDep.functions.postprocessresults import Results

from datetime import datetime

start_time = datetime.now()


FY_WGT = 0.0  # determines the fission yield wieghting
VOL = 332097.8  # cm^3


# -----------------------------------------------------------------------------
#                            DATA GENERATION
# -----------------------------------------------------------------------------
xsTable, fields = ReadCsv("./bootstrap.csv")
bootstrap = TransmutationData(libraryFlag=True, wgtFY=1.0)
bootstrap.ReadData(ID=xsTable[:, 0], sig_f=xsTable[:, 3], sig_c=xsTable[:, 2],
                   sig_alpha=xsTable[:, 4], flagBarns=True)

xsTable, fields = ReadCsv("./tempramp.csv")
tempramp = TransmutationData(libraryFlag=True, wgtFY=1.0)
tempramp.ReadData(ID=xsTable[:, 0], sig_f=xsTable[:, 3], sig_c=xsTable[:, 2],
                  sig_alpha=xsTable[:, 4], flagBarns=True)

xsTable, fields = ReadCsv("./fullthrust.csv")
fullthrust = TransmutationData(libraryFlag=True, wgtFY=1.0)
fullthrust.ReadData(ID=xsTable[:, 0], sig_f=xsTable[:, 3], sig_c=xsTable[:, 2],
                    sig_alpha=xsTable[:, 4], flagBarns=True)


# -------------------------------------------------------------------------
#                            DEPLETION
# -------------------------------------------------------------------------
dep = MainDepletion([0.0, 5.5, 30.0], bootstrap, tempramp, fullthrust)
# define metadata (steps, flux, and so on)
power = 1E+6*np.array([16.545, 118.49, 272.52, 330.22, 272.52, 214.82, 118.49])
dt = np.array([5.5, 24.5, 7., 1800., 7., 180., 40.])
dep.SetDepScenario(power=power, timeUnits="seconds", timesteps=dt)
# set initial composition
dep.SetInitialComposition(xsTable[:, 0], xsTable[:, 1], vol=VOL)
# solve the Bateman equations
dep.SolveDepletion(method="adaptive", xsinterp=False)
# Post depletion analysis
dep.DecayHeat()
dep.Radiotoxicity()
dep.Activity()
dep.Mass()


# Post-process results
# -----------------------------------------------------------------------------
res = Results(dep)
res.plot("Nt", isotopes=[541350])
res.plot("totalQt", norm=1E+6, ylabel="Total Decay Heat, MW",
         pltType="semilogx")

res.plot("Qt", isotopes=[531350, 541350], norm=1E+6, ylabel="Total Decay Heat, MW")

res.plot("flux", ylabel="Flux, n/cm2/s", pltType="semilogx")
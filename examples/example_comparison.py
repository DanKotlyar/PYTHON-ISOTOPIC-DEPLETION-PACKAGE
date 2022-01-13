# -*- coding: utf-8 -*-
"""simplexenonchains

Use custom depletion chains for xenon and samarium and perform
depletion calculations.

Created on Sun Oct 11 07:30:00 2021 @author: Dan Kotlyar
Last updated on Sat Oct 16 15:30:00 2021 @author: Dan Kotlyar

"""
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

from pyIsoDep.functions.maindepletionsolver import MainDepletion
from pyIsoDep.functions.generatedata import TransmutationData
from pyIsoDep.functions.postprocessresults import Results

ID = [531350, 541350, 611490, 621490, 922350, 922380]
sig_c = [6.8, 250537.62, 132.47, 6968.75, 5.0, 8.0]
sig_f = [0.0, 0.0000000, 0.000, 0.00000, 97., 3.8]
kappa = [0.0, 0.0000000, 0.000, 0.00000, 202.44, 202.44]
N0 = [0.0, 0.0000000, 0.000, 0.00000, 6.43230E-04, 2.58062E-03]

#    531350, 541350, 611490, 621490, 922350, 922380
mtxFY = [
    [0.0000, 0.0000, 0.0000, 0.0000, 0.06306, 0.06306],  # 531350
    [0.0000, 0.0000, 0.0000, 0.0000, 0.00248, 0.00248],  # 541350
    [0.0000, 0.0000, 0.0000, 0.0000, 0.01100, 0.01100],  # 611490
    [0.0000, 0.0000, 0.0000, 0.0000, 0.00000, 0],  # 621490
    [0.0000, 0.0000, 0.0000, 0.0000, 0.00000, 0],  # 922350
    [0.0000, 0.0000, 0.0000, 0.0000, 0.00000, 0],  # 922380
    ]

volume = 332097.750  # volume in cm**3

timepoints = np.linspace(0, 48.0, 24)
power = 330000000.*np.ones(len(timepoints)-1)  # Watts

start_time = datetime.now()


# -----------------------------------------------------------------------------
#                            DATA GENERATION
# -----------------------------------------------------------------------------
# Reset the data container
data = TransmutationData(libraryFlag=True, wgtFY=1.0)
# Feed cross sections into the container
data.ReadData(ID, sig_f=sig_f, sig_c=sig_c, fymtx=mtxFY, EfissMeV=kappa)
# Condense the data only to specific set of isotopes
data.Condense(ID)

# -----------------------------------------------------------------------------
#                            DEPLETION
# -----------------------------------------------------------------------------
dep = MainDepletion(0.0, data)
dep.SetDepScenario(power=power, timeUnits="hours", timepoints=timepoints)
# set initial composition
dep.SetInitialComposition(ID, N0, vol=volume)
# solve the Bateman equations
dep.SolveDepletion(method="cram")
# Post depletion analysis
dep.DecayHeat()
dep.Radiotoxicity()
dep.Activity()
dep.Mass()
print("cram {:0.3e}".format(dep._solveTime))
resCram = Results(dep)

dep = MainDepletion(0.0, data)
dep.SetDepScenario(power=power, timeUnits="hours", timepoints=timepoints)
# set initial composition
dep.SetInitialComposition(ID, N0, vol=volume)
# solve the Bateman equations
dep.SolveDepletion(method="expm")
# Post depletion analysis
dep.DecayHeat()
dep.Radiotoxicity()
dep.Activity()
dep.Mass()
print("expm {:0.3e}".format(dep._solveTime))
resExpm = Results(dep)

dep = MainDepletion(0.0, data)
dep.SetDepScenario(power=power, timeUnits="hours", timepoints=timepoints)
# set initial composition
dep.SetInitialComposition(ID, N0, vol=volume)
# solve the Bateman equations
dep.SolveDepletion(method="odeint")
# Post depletion analysis
dep.DecayHeat()
dep.Radiotoxicity()
dep.Activity()
dep.Mass()
print("oden {:0.3e}".format(dep._solveTime))
resOdeint = Results(dep)


# Reset the data container
timepoints = np.linspace(0, 48.0, 2)
power = 330000000.*np.ones(len(timepoints)-1)  # Watts
data = TransmutationData(libraryFlag=True, wgtFY=1.0)
# Feed cross sections into the container
data.ReadData(ID, sig_f=sig_f, sig_c=sig_c, fymtx=mtxFY, EfissMeV=kappa)
# Condense the data only to specific set of isotopes
data.Condense(ID)
dep = MainDepletion(0.0, data)
dep.SetDepScenario(power=power, timeUnits="hours", timepoints=timepoints)
# set initial composition
dep.SetInitialComposition(ID, N0, vol=volume)
# solve the Bateman equations
dep.SolveDepletion(method="adaptive")
# Post depletion analysis
dep.DecayHeat()
dep.Radiotoxicity()
dep.Activity()
dep.Mass()
print("adpt {:0.3e}".format(dep._solveTime))
resAdapt = Results(dep)



plt.plot(resCram.timepoints, resCram.Nt[1,:], label="Cram")
plt.plot(resExpm.timepoints, resExpm.Nt[1,:], label="Expm")
plt.plot(resOdeint.timepoints, resOdeint.Nt[1,:], label="Odeint")
plt.scatter(resAdapt.timepoints, resAdapt.Nt[1,:], label="Adaptive", s=10, marker="s")
plt.ylabel("135Xe Conc, atoms/b-cm")
plt.xlabel("time, hours")
plt.legend()
plt.show()
plt.close()


plt.plot(resCram.timepoints, resCram.Nt[0,:], label="Cram")
plt.plot(resExpm.timepoints, resExpm.Nt[0,:], label="Expm")
plt.plot(resOdeint.timepoints, resOdeint.Nt[0,:], label="Odeint")
plt.scatter(resAdapt.timepoints, resAdapt.Nt[0,:], label="Adaptive", s=10, marker="s")
plt.ylabel("135I Conc, atoms/b-cm")
plt.xlabel("time, hours")
plt.legend()
plt.show()
plt.close()

plt.plot(resCram.timepoints, resCram.Nt[3,:], label="Cram")
plt.plot(resExpm.timepoints, resExpm.Nt[3,:], label="Expm")
plt.plot(resOdeint.timepoints, resOdeint.Nt[3,:], label="Odeint")
plt.scatter(resAdapt.timepoints, resAdapt.Nt[3,:], label="Adaptive", s=10, marker="s")
plt.ylabel("149Sm Conc, atoms/b-cm")
plt.xlabel("time, hours")
plt.legend()
plt.show()
plt.close()

plt.plot(resCram.timepoints, resCram.Nt[4,:], label="Cram")
plt.plot(resExpm.timepoints, resExpm.Nt[4,:], label="Expm")
plt.plot(resOdeint.timepoints, resOdeint.Nt[4,:], label="Odeint")
plt.scatter(resAdapt.timepoints, resAdapt.Nt[4,:], label="Adaptive", s=10, marker="s")
plt.ylabel("235U Conc, atoms/b-cm")
plt.xlabel("time, hours")
plt.legend()
plt.show()
plt.close()

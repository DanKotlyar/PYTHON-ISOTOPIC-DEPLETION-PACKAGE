# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 16:19:50 2021

@author: mkrecicki3
"""
import pytest
import numpy as np
import matplotlib.pyplot as plt
from pyIsoDep.functions.maindepletionsolver import MainDepletion
from pyIsoDep.functions.generatedata import TransmutationData
from pyIsoDep.functions.postprocessresults import Results

def _xenonAnalyticMatt(time, I0, X0, sigaX, SIG_F, flx, lmdaI, lmdaX, fyI, 
                       fyX):
    B2CM = 10E-24
    sigaX = sigaX * B2CM
    
    It = ((fyI * SIG_F * flx)/lmdaI)*(1-np.exp(-1*lmdaI*time)) + (I0*np.exp(-1*lmdaI*time))
    
    C1 = ((((lmdaI+lmdaX)*SIG_F*flx)) / (lmdaX + sigaX*flx) )*(1 - np.exp(-1*(lmdaX + sigaX*flx)*time))
    C2 = ( (lmdaI*SIG_F*flx - lmdaI*I0) / (lmdaX - lmdaI + sigaX*flx ) ) * (np.exp(-1*(lmdaX+sigaX*flx)*time) - np.exp(-1*lmdaI*time))
    C3 = X0*np.exp(-1*(lmdaX + sigaX*flx)*time)
    
    Xt = C1 + C2 + C3
    
    return It, Xt


def analyticalIodine():
    pass


def analyticalXenon():
    pass


def _xenonAnalytic(timepoints, I0, X0, sigaX, SIG_F, flx, lmdaI, lmdaX, fyI,
                   fyX):
    """analytical solution for Iodine -> Xenon chain"""
    B2CM = 10E-24  # conversion from barns to cm**2
        
    I135 = ((fyI*SIG_F*flx) / lmdaI) * (1-np.exp(-lmdaI*timepoints)) 
    # + I0*np.exp(-lmdaI*timepoints)
        
    Xt = (fyI+fyX)*SIG_F*flx/(lmdaX+sigaX*flx*B2CM)*(
        1-np.exp((-lmdaX-sigaX*flx*B2CM)*timepoints)) +\
        (fyI*SIG_F*flx-lmdaI*I0)/(lmdaX-lmdaI+sigaX*flx*B2CM)*(
            np.exp((-lmdaX-sigaX*B2CM*flx)*timepoints)-np.exp(-lmdaI*timepoints)) +\
        X0*np.exp((-lmdaX-sigaX*flx*B2CM)*timepoints)
        
    return I135*B2CM, Xt*B2CM

   
def test_analytical_xenon():
    """compare analytical and numerical solutions"""
    
    timepoints = np.array([   0,  300,  600,  900, 1200, 1500, 1800, 2100, 2400]) #time in seconds
    volume = 1.0  # volume in cm**3
    flux = 4.91024e+14*np.ones(len(timepoints)-1)  #flux
        
    ID =    [531350, 541350,     922350]
    sig_c = [0.0,    250537.62,  0.0]
    sig_f = [0.0,    0.0,        97.]
    kappa = [0.0,    0.0,        202.44]
    N0 =    [0.0,    0.0,        6.43230E-04]
    
    #    531350, 541350, 922350
    mtxFY = [
        [0.0000, 0.0000, 0.06306],  # 531350
        [0.0000, 0.0000, 0.00248],  # 541350
        [0.0000, 0.0000, 0.00000]]  # 922350
      
    #    531350,         541350,      922350
    mtxDY = [
        [-2.93061e-05,   0.0000,      0.0000],  # 531350
        [2.93061e-05,   -2.10657e-05, 0.0000],  # 541350
        [0.0000,         0.0000,      0.00000],  # 922350
        ]
    
    data = TransmutationData(libraryFlag=True, wgtFY=1.0) # Reset the data container
    data.ReadData(ID, sig_f=sig_f, sig_c=sig_c, fymtx=mtxFY, EfissMeV=kappa,
                  decaymtx=mtxDY) # Feed cross sections into the container
    data.Condense(ID) # Condense the data only to specific set of isotopes
    dep = MainDepletion(0.0, data)
    
    
    timeUnits="seconds"
    dep.SetInitialComposition(ID, N0, vol=volume)
    dep.SetDepScenario(flux=flux, timeUnits=timeUnits, timepoints=timepoints)
    dep.SolveDepletion(method="cram") # solve the Bateman equations
    dep.DecayHeat()
    dep.Radiotoxicity()
    dep.Activity()
    dep.Mass()
    res = Results(dep)
    
    ICalc = res.getvalues(isotopes=[531350], attribute="Nt")
    XeCalc = res.getvalues(isotopes=[541350], attribute="Nt")
    
    I135, XeRef = _xenonAnalyticMatt(timepoints,
                                 float(N0[0]),
                                 float(N0[1]),
                                 float(sig_c[1]),
                                 float(N0[2]*sig_f[2]),
                                 float(flux[0]),
                                 float(abs(mtxDY[0][0])),
                                 float(abs(mtxDY[1][1])),
                                 float(mtxFY[0][2]), 
                                 float(mtxFY[1][2]))
    
    plt.plot(timepoints/60, I135*1E-24, label="Analytical I")
    plt.plot(timepoints/60, XeRef*1E-24, label="Analytical Xe")
    plt.scatter(timepoints/60, ICalc, label="Cram I", s=12, marker="o")
    plt.scatter(timepoints/60, XeCalc, label="Cram Xe", s=12, marker="x")
    plt.xlabel("Time, minutes")
    plt.ylabel("Conc., atoms/b-cm")
    #plt.yscale("log")
    plt.legend()
    plt.show()
    plt.close()
    
    return res, dep, data

    
def test_analytical_xenon_deacy():
    
    I0 = 4.47736495160465e-08
    X0 = 1.2989627590704352e-09
    
    ID = [531350, 541350]
    N0 = [I0, X0]
    
    time = np.linspace(0, 80*60*60, 60)
    
    flx = 0.0 
    lmdaI = 2.93061e-05
    lmdaX = 2.10657e-05
    fyI = 0.0
    fyX = 0.0
    sigaX = 0.0
    SIG_F = 0.0
    
    vol = 1.0
    
    #    531350, 541350
    mtxDY = [
            [ -lmdaI, 0.0],  # 531350
            [ lmdaI, -lmdaX]]  # 541350   
    
    #obtain analytical solution
    IRef, XeRef = _xenonAnalyticMatt(time, I0, X0, sigaX, SIG_F, flx, lmdaI,
                                     lmdaX, fyI, fyX)
        
    #create numerical solution
    data = TransmutationData(libraryFlag=True)
    data.ReadData(ID, sig_f=np.array([0.0, 0.0]), sig_c=np.array([0.0, 0.0]),
                  decaymtx=mtxDY)
    data.Condense(ID)
        
    dep = MainDepletion(0.0, data)
    dep.SetDepScenario(timeUnits="seconds", timepoints=time)
    dep.SetInitialComposition(ID, N0, vol=vol)
    dep.SolveDecay(method="cram")
    resCram = Results(dep)
    
    dep = MainDepletion(0.0, data)
    dep.SetDepScenario(timeUnits="seconds", timepoints=time)
    dep.SetInitialComposition(ID, N0, vol=vol)
    dep.SolveDecay(method="expm")
    resExpm = Results(dep)
    
    dep = MainDepletion(0.0, data)
    dep.SetDepScenario(timeUnits="seconds", timepoints=time)
    dep.SetInitialComposition(ID, N0, vol=vol)
    dep.SolveDecay(method="odeint")
    resOde = Results(dep)
    
    ICram = resCram.getvalues(isotopes=[531350], attribute="Nt").flatten()
    XeCram = resCram.getvalues(isotopes=[541350], attribute="Nt").flatten()
    IExpm = resExpm.getvalues(isotopes=[531350], attribute="Nt").flatten()
    XeExpm = resExpm.getvalues(isotopes=[541350], attribute="Nt").flatten()
    IOde = resOde.getvalues(isotopes=[531350], attribute="Nt").flatten()
    XeOde = resOde.getvalues(isotopes=[541350], attribute="Nt").flatten()
    
    #plot solution 
    plt.plot(time/(60*60), IRef, label="Analytical I")
    plt.plot(time/(60*60), XeRef, label="Analytical Xe")
    
    plt.scatter(time/(60*60), ICram, marker="<", label="Cram I", c="C2")
    plt.scatter(time/(60*60), XeCram, marker="<", label="Cram Xe", c="C3")
    
    #plt.scatter(time/(60*60), IExpm, marker="x", label="Expm I", c="C4")
    #plt.scatter(time/(60*60), XeExpm, marker="x", label="Expm Xe", c="C5")
    
    plt.scatter(time/(60*60), IOde, marker="v", label="Odeint I", c="C6")
    plt.scatter(time/(60*60), XeOde, marker="v", label="Odeint Xe", c="C7")
    
    plt.xlabel("Time, hours")
    plt.ylabel("Conc., atoms/b-cm")
    plt.legend()
    plt.show()
    plt.close()


#test_analytical_xenon_deacy()
res, dep, data = test_analytical_xenon()
# -*- coding: utf-8 -*-
"""test_resultsreconstruction

Tests that a single depletion step is carried out properly.
The entire sequence from cross section generation to depletion execution is
tested. Results are compared against pre-generated data using a different code.

Created on Thu Oct 28 08:59:44 2021 @author: Matt Krecicki
Last updated onThu Oct 28 08:59:44 2021 @author: Matt Krecicki

"""
import pytest
import numpy as np
from pyIsoDep.functions.maindepletionsolver import MainDepletion
from pyIsoDep.functions.generatedata import TransmutationData
from pyIsoDep.functions.postprocessresults import Results

from pyIsoDep.tests.pregenerated_xs import flux, ID, N0, sig_c,\
    sig_c2m, sig_n2n, sig_n3n, sig_f, compareNt
    
    
# -----------------------------------------------------------------------------
#                            DATA GENERATION
# -----------------------------------------------------------------------------

# Reset the data container
data = TransmutationData(libraryFlag=True, wgtFY=1.0)
# Feed cross sections into the container
data.ReadData(ID, sig_f=sig_f, sig_c=sig_c, sig_c2m=sig_c2m,
              sig_n2n=sig_n2n, sig_n3n=sig_n3n, flagBarns=False)

def test_reconstruction():
    """Test that ensure results reconstruction is correct"""
    # -------------------------------------------------------------------------
    #                            DEPLETION
    # -------------------------------------------------------------------------
    dep = MainDepletion(0.0, data)
    # define metadata (steps, flux, and so on)
    dep.SetDepScenario(power=None, flux=[flux], timeUnits="seconds",
                       timesteps=[6.630851880276299780234694480896E+05],
                       timepoints=None)
    # set initial composition
    dep.SetInitialComposition(ID, N0, vol=1.0)
    # solve the Bateman equations
    dep.SolveDepletion(method="cram")
    # Post depletion analysis
    dep.DecayHeat()
    dep.Radiotoxicity()
    dep.Activity()
    dep.Mass()
    
    #export results to hdf5 file
    res = Results(dep)
    res.export("test.h5")
    
    #reconstruct results from hdf5 file
    res2 = Results("test.h5")
    
    #compare exported and reconstructed results
    assert res.N0 == pytest.approx(res2.N0, rel=0.001)
    assert res.flagPower == res2.flagPower
    assert res.flux == pytest.approx(res2.flux, rel=0.001)
    assert res.totalQt == pytest.approx(res2.totalQt, rel=0.001)
    assert res.decaymtx[25,:] == pytest.approx(res2.decaymtx[25,:], rel=0.001)
    assert res.Nt[30,:] == pytest.approx(res2.Nt[30,:], rel=0.001)  
    assert res._xsDataSets[0.0][50,:] ==\
        pytest.approx(res2._xsDataSets[0.0][50,:], rel=0.001) 
       
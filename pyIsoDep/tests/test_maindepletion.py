# -*- coding: utf-8 -*-
"""test_maindepletion

Tests that a single depletion step is carried out properly.
The entire sequence from cross section generation to depletion execution is
tested. Results are compared against pre-generated data using a different code.

Created on Sun Oct 11 07:30:00 2021 @author: Dan Kotlyar
Last updated on Sun Oct 17 17:30:00 2021 @author: Dan Kotlyar

"""

import pytest

from pyIsoDep.functions.maindepletionsolver import MainDepletion
from pyIsoDep.functions.generatedata import TransmutationData

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


def test_depletion():
    """Test that depletion is carried out properly"""
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

    assert dep.Nt[1656, 1] == pytest.approx(compareNt[1656], rel=0.001)


def test_badMainDepletion():
    """Errors for the main depletion definitions"""

    with pytest.raises(ValueError, match="Time Frames*"):
        MainDepletion(0.0, data, data)

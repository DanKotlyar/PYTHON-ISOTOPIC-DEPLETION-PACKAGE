# -*- coding: utf-8 -*-
"""test_xsinterface

Tests that cross sections are properly stored on the interface and that
the interpolation is carried out properly.

Created on Wed Feb 02 09:10:00 2022 @author: Dan Kotlyar
Last updated on Wed Feb 02 09:10:00 2022 @author: Dan Kotlyar

"""

import numpy as np

import pytest

from pyIsoDep.functions.generatedata import TransmutationData
from pyIsoDep.functions.xsinterface import XsInterface
from pyIsoDep.tests.interp2d import _interp2D


def test_2d_interpolation():
    """Test that 2-dim interpolation works properly"""

    # -------------------------------------------------------------------------
    #      Build a fake fission and capture cross sections only fot U235
    # -------------------------------------------------------------------------
    pressures = np.array([2, 3, 4, 5])
    temperatures = np.array([500, 1000, 1500, 2000, 2500])
    npert = len(pressures) * len(temperatures)
    xssets = []  # store all cross section sets
    states = []  # store all operational states

    sigf0 = 35.0  # fission cross section for U235
    counter = 0
    dsig = 0.5
    sigVec = np.zeros(npert)
    for pressure in pressures:
        for temperature in temperatures:
            sigVec[counter] = sigf0+dsig*counter
            data = TransmutationData(libraryFlag=True, wgtFY=1.0)
            # Feed cross sections into the container
            data.ReadData(ID=[922350], sig_f=[sigf0+dsig*counter],
                          sig_c=[sigf0+dsig*counter])
            xssets.append(data)
            states.append([pressure, temperature])
            counter += 1

    # -------------------------------------------------------------------------
    #                  XS Interface (BiLinear Interpolation/Extrp.)
    # -------------------------------------------------------------------------

    xs0 = XsInterface(numdepn=2, numpert=len(xssets), states=states,
                      xssets=xssets, extrpFlag=True)

    xTrace = [3.25, 3.75, 4.0, 4.25, 4.75]
    yTrace = [925, 1230, 1625, 2210, 2395]

    timepoints, intrXs = xs0.setTimeTrace([0, 1, 2, 3, 4], xTrace, yTrace)

    intrXs = list(intrXs)
    predVals = np.zeros(5)
    expVals = np.zeros(5)

    states = np.array((states))
    for idx in range(5):
        predVals[idx] = 1E+24*intrXs[idx].xsData[1656, 2]
        expVals[idx] = _interp2D(xTrace[idx], states[:, 0],
                                 yTrace[idx], states[:, 1], sigVec)

    assert predVals == pytest.approx(expVals, rel=0.001)


def test_1d_interpolation():
    """Test that 1-dim interpolation works properly"""

    temperatures = np.array([500, 1000, 1500, 2000, 2500])
    npert = len(temperatures)
    xssets = []  # store all cross section sets
    states = []  # store all operational states

    sigf0 = 35.0  # fission cross section for U235
    counter = 0
    dsig = 0.5
    sigVec = np.zeros(npert)
    for temperature in temperatures:
        sigVec[counter] = sigf0+dsig*counter
        data = TransmutationData(libraryFlag=True, wgtFY=1.0)
        # Feed cross sections into the container
        data.ReadData(ID=[922350], sig_f=[sigf0+dsig*counter],
                      sig_c=[sigf0+dsig*counter])
        xssets.append(data)
        states.append([temperature])
        counter += 1

    # -------------------------------------------------------------------------
    #                  XS Interface (BiLinear Interpolation/Extrp.)
    # -------------------------------------------------------------------------

    xs0 = XsInterface(numdepn=1, numpert=len(xssets), states=states,
                      xssets=xssets, extrpFlag=True)
    xTrace = [750]

    timepoints, intrXs = xs0.setTimeTrace([0], xTrace)

    intrXs = list(intrXs)
    states = np.array((states))
    predVal = 1E+24*intrXs[0].xsData[1656, 2]
    expVal = 0.5*(sigVec[0]+sigVec[1])

    assert predVal == pytest.approx(expVal, rel=0.001)

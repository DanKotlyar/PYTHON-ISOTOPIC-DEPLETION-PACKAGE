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
from pyIsoDep.functions.xsinterface import XsInterface

from datetime import datetime

start_time = datetime.now()


FY_WGT = 1.0  # determines the fission yield wieghting
VOL = 332097.8  # cm^3


# -----------------------------------------------------------------------------
#                            DATA GENERATION
# -----------------------------------------------------------------------------
xsTable, fields = ReadCsv("./bootstrap.csv")
bootstrap = TransmutationData(libraryFlag=True, wgtFY=1.0)
bootstrap.ReadData(ID=xsTable[:, 0], sig_f=xsTable[:, 3], sig_c=xsTable[:, 2],
                   sig_a=xsTable[:, 4], flagBarns=True)

xsTable, fields = ReadCsv("./tempramp.csv")
tempramp = TransmutationData(libraryFlag=True, wgtFY=1.0)
tempramp.ReadData(ID=xsTable[:, 0], sig_f=xsTable[:, 3], sig_c=xsTable[:, 2],
                  sig_a=xsTable[:, 4], flagBarns=True)

xsTable, fields = ReadCsv("./fullthrust.csv")
fullthrust = TransmutationData(libraryFlag=True, wgtFY=1.0)
fullthrust.ReadData(ID=xsTable[:, 0], sig_f=xsTable[:, 3], sig_c=xsTable[:, 2],
                    sig_a=xsTable[:, 4], flagBarns=True)


# -------------------------------------------------------------------------
#                            XS Interface (Linear Interpolation/Extrp.)
# -------------------------------------------------------------------------


xs = XsInterface(numdepn=1, numpert=3, states=[[500], [600], [700]],
                 xssets=[bootstrap, tempramp, fullthrust])

timepoints, xsTimeSets = xs.setTimeTrace([0, 3.5, 7.0], [525, 550, 575])


# -------------------------------------------------------------------------
#                            XS Interface (BiLinear Interpolation/Extrp.)
# -------------------------------------------------------------------------

xs = XsInterface(numdepn=2, numpert=6, states=[[500, 2], [500, 3], [500, 4],
                                                [600, 2], [600, 3], [600, 4]],
                  xssets=[bootstrap, bootstrap, bootstrap, bootstrap, bootstrap,
                          bootstrap])

timepoints, xsTimeSets = xs.setTimeTrace([0, 3.5], [500, 550], [3.0, 3.5])


a = 1

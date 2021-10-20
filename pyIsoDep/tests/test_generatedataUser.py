# -*- coding: utf-8 -*-
"""test_generatedata

Tests that the transmutation matrix is built properly.

Created on Sun Oct 11 07:30:00 2021 @author: Dan Kotlyar
Last updated on Fri Oct 15 14:30:00 2021 @author: Dan Kotlyar

"""

import numpy as np
import timeit

from pyIsoDep.functions.generatedata import TransmutationData

# -----------------------------------------------------------------------------
# User defined data
# -----------------------------------------------------------------------------
ID = [531350, 541350, 611490, 621490, 922350, 922380]
sig_c = [6.8, 250537.62, 132.47, 6968.75, 5.0, 8.0]
sig_f = [0.0, 0.0000000, 0.000, 0.00000, 97., 3.8]
kappa = [0.0, 0.0000000, 0.000, 0.00000, 202.44, 202.44]

#    531350, 541350, 611490, 621490, 922350, 922380
fymtx = [
    [0.0000, 0.0000, 0.0000, 0.0000, 0.06306, 0.06306],  # 531350
    [0.0000, 0.0000, 0.0000, 0.0000, 0.00248, 0.00248],  # 541350
    [0.0000, 0.0000, 0.0000, 0.0000, 0.01100, 0.01100],  # 611490
    [0.0000, 0.0000, 0.0000, 0.0000, 0.00000, 0],  # 621490
    [0.0000, 0.0000, 0.0000, 0.0000, 0.00000, 0],  # 922350
    [0.0000, 0.0000, 0.0000, 0.0000, 0.00000, 0],  # 922380
    ]

decaymtx = np.zeros((6, 6))

tic = timeit.timeit()

# -----------------------------------------------------------------------------
# Reset the data container
# -----------------------------------------------------------------------------
data = TransmutationData(libraryFlag=False, h5path=datapath)

# -----------------------------------------------------------------------------
# Feed cross sections into the container
# -----------------------------------------------------------------------------
data.ReadData(ID, sig_f=sig_f, sig_c=sig_c, fymtx=fymtx, decaymtx=decaymtx)


toc = timeit.timeit()

print("Reading data = {} seconds ".format(toc-tic))

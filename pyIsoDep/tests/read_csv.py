# -*- coding: utf-8 -*-
"""read_csv

Read the different csv files

Created on Mon Oct 11 21:30:00 2021 @author: Dan Kotlyar
Last updated on Mon Oct 11 21:45:00 2021 @author: Dan Kotlyar

"""

import numpy as np
import pandas as pd


def ReadCsv(csvFile):

    data = pd.read_csv('bootstrap.csv')
    ID = np.array(data['ZAID'], dtype=int)
    xsTypes = np.array(data['MT'], dtype=int)
    xsVals = np.array(data["XS [barns]"], dtype=float)
    N0 = np.array(data["N0 [atoms/b-cm]"], dtype=float)

    fullID = np.unique(ID)  # unique isotopes
    nIsotopes = len(fullID)
    # 1-ID, 2-ND, 3-cap, 4-fiss, 5-(n,alpha)
    xsTable = np.zeros((nIsotopes, 5))
    xsTable[:, 0] = fullID

    # obtain all the cross section types
    numMTs = np.array([102, 18, 107])

    for idx, numMT in enumerate(numMTs):
        vals, idxFull, idx0 =\
            np.intersect1d(fullID, ID[xsTypes == numMT], assume_unique=False,
                           return_indices=True)
        if idx == 0:
            xsTable[idxFull, 1] = N0[xsTypes == numMT][idx0]
        xsTable[idxFull, idx+2] = xsVals[xsTypes == numMT][idx0]

    idxFields = {"ID": 0, "N0": 1, "sig_c": 2, "sig_alpha": 3, "sig_f": 4}

    return xsTable, idxFields

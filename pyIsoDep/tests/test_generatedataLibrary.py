# -*- coding: utf-8 -*-
"""test_generatedata

Tests that the transmutation matrix is built properly.
Data is read directly from a pre-generated default library.

Created on Sun Oct 11 07:30:00 2021 @author: Dan Kotlyar
Last updated on Sat Oct 16 00:30:00 2021 @author: Dan Kotlyar

"""

from pyIsoDep.functions.generatedata import TransmutationData

# import pre-generated data containing set of cross sections
from depletionPackage.tests.pregenerated_xs import ID, sig_c, sig_c2m,\
    sig_n2n, sig_n3n, sig_f

import timeit

tic = timeit.timeit()
# -----------------------------------------------------------------------------
# Reset the data container
# -----------------------------------------------------------------------------
data = TransmutationData(libraryFlag=True, wgtFY=1.0)

# -----------------------------------------------------------------------------
# Feed cross sections into the container
# -----------------------------------------------------------------------------
data.ReadData(ID, sig_f=sig_f, sig_c=sig_c, sig_c2m=sig_c2m,
              sig_n2n=sig_n2n, sig_n3n=sig_n3n, flagBarns=False)


toc = timeit.timeit()
print("Time elapsed for reading data = {} seconds ".format(toc-tic))

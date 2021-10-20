# -*- coding: utf-8 -*-
"""test_loaddecaydata

Tests that data is loaded properly

Created on Sat Oct 10 05:30:00 2021 @author: Dan Kotlyar
Last updated on Sat Oct 10 05:30:00 2021 @author: Dan Kotlyar

"""

from linkdirs import datapath

from loaddecaydata import DecayData

import timeit

tic = timeit.timeit()

# load the decay data
table = DecayData(datapath)

# list all the properties in the data file
properties = table.properties()

# get description for a specific propery
ptydesc = table.description("AW")

# Obtain values for a specific property
ID = table.getvalues("IDlist")
AW = table.getvalues("AW")
Q = table.getvalues("Q")
lmbda = table.getvalues("lambda")
lmbdaMtx = table.getvalues("decayMatrix")
fastFY = table.getvalues("fastFY")
thermalFY = table.getvalues("thermalFY")
ingestion = table.getvalues("ingestion")
inhalation = table.getvalues("inhalation")

toc = timeit.timeit()

print("{} seconds ".format(toc-tic))

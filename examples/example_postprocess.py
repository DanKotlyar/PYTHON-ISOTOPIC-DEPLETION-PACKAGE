# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 19:45:41 2021

@author: matt krecicki

"""
from pyIsoDep.functions.postprocessresults import Results


res = Results("scenario.h5")

rankQt = res.rank(parameter="Qt")
rankIngst = res.rank(parameter="toxicityIngestion")
rankInhal = res.rank(parameter="toxicityInhalation")
rankAt = res.rank(parameter="At")
rankRho = res.rank(parameter="reactivity")
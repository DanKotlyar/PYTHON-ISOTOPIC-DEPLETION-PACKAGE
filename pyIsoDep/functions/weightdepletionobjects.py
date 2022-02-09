"""weightdepletionobjects

Allows to weight different depletion objects so that volume averages
can be created.


Created on Wed Feb 09 15:30:00 2022 @author: Dan Kotlyar
Last updated on Wed Feb 09 15:30:00 2022 @author: Dan Kotlyar

"""

import copy

from pyIsoDep.functions.header import NOT_WEIGHT_ATTR


def WeightDepObjects(*argv):
    """weight different depletion objects"""

    volT = 0  # total volume
    # calculate the total volume
    for idx, objDep in enumerate(argv):
        volT += objDep.volume

    for idx, objDep in enumerate(argv):
        if idx == 0:
            # create a weighted depletion set
            wgtDepObj = copy.copy(objDep)
            wgtDepObj.volume = volT

        for key, values in objDep.__dict__.items():
            if key not in NOT_WEIGHT_ATTR:
                setvals = values * objDep.volume / volT
                if idx > 0:
                    setvals += getattr(wgtDepObj, key)
                setattr(wgtDepObj, key, setvals)

    return wgtDepObj

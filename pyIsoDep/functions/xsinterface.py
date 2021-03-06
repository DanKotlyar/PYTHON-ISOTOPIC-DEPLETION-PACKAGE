"""xsinterface

An interface to store all the cross-section sets for a range of operational
conditions.
The interface also contains interpolation techniques to allow constructing
a new cross section set for given operational conditions.


Created on Tue Nov 16 18:00:00 2021 @author: Dan Kotlyar
Last updated on Sat Jan 15 12:00:00 2021 @author: Dan Kotlyar

"""

import copy

import numpy as np

from pyIsoDep.functions.checkerrors import _isint, _exp2dshape,\
    _isequallength, _anynegative, _inrange, _islist

from pyIsoDep.functions.header import DATA_ATTR, INTRP_ATTR


class XsInterface:
    """A class to interact with pre-generated cross section sets.

    This interface is used to link between cross sections generated in
    branch-off calculations for multiple operational points and the actual
    operational regime. The interface allows to generate on-the-fly cross
    section sets that represent different time-dependent scenarios.

    Parameters
    ----------
    numdepn : int
        Number of dependencies, e.g., 2 for when only pressure chamber and
        temperatures are considered.
    numpert : int
        Total number of perturbation points for which cross sections were
        generated.
    states : 2-dim array
        The actual operational states. Rows represent different operational
        state, and columns contain the actual operational values, e.g.,
        [[P1, T1],[P2, T2],[P3, T3],...]. The states must all be unique.
    xssets : tuple of TransmutationData objects
        The number of objects must correspond to the number of perturbations

    Attributes
    ----------
    _xsDataSets : dictionary
        Keys represent the timeframes and values the TransmutationData objects
    _timeframes : array
        Time frames for all the stored objects (can be used to interpolate xs)

    Raises
    ------
    TypeError
        If any of ``numdepn`` or ``numpert`` is not int.
        If ``states`` is not an array.
        If ``xssets`` is not a tuple.
        If ``extrpFlag`` is not boolean.
    ValueError
        If ``numdepn`` is above 3.
        If the number of perturbabtions within ``states`` is incorrect.

    Examples
    --------
    >>> xs = XsInterface(numdepn=2, numpert=3,
                         states=[[500, 5],[600, 5],[700, 6]],
                         xssets=(xsset1, xsset2, xsset3))

    """

    def __init__(self, numdepn, numpert, states, xssets):
        """reset values with a complete list of all the nuclides"""

        # Check for potential errors
        # ---------------------------------------------------------------------
        _isint(numdepn, "Number of dependencies")
        _isint(numpert, "Number of perturbations")
        states = np.array((states))
        _exp2dshape(states, (numpert, numdepn), "States")
        _islist(xssets, "XS Sets")
        _isequallength(xssets, numpert, "XS Sets List")
        # Check values
        if numdepn > 2:
            raise ValueError("Current version supports up to 2 dependencies "
                             "and not {}".format(numdepn))
        # Check that each data set contains the required attributes
        # ---------------------------------------------------------------------
        # Loop over all the required fields/attributes
        for attr in DATA_ATTR.keys():
            for idx, data in enumerate(xssets):
                refId = data.fullId
                if not hasattr(data, attr):
                    print("No attribute <{}> in data".format(attr))
                Id0 = data.fullId
                if (refId != Id0).any():
                    raise ValueError(
                        "fullId for index set {} is not identical to the 0th "
                        "set.".format(idx))

        # Check that data is complete/squared
        # ---------------------------------------------------------------------
        uniqStates = 1
        for i in range(numdepn):
            uniqStates *= len(np.unique(states[:, i]))
        if uniqStates != numpert:
            raise ValueError("Complete square matrix must be provided for "
                             "states. Number of unique states={}, however {} "
                             "perturbations were provided".format(uniqStates,
                                                                  numpert))

        # Store the user-provided data unsorted data
        # ---------------------------------------------------------------------
        self._states = states
        self._xssets = xssets

        # Sort data
        # ---------------------------------------------------------------------
        uniqStates, lenStates, xssetsMtx =\
            self._sortsets(numdepn, numpert, states, xssets)

        # Store the data
        # ---------------------------------------------------------------------
        self.numdepn = numdepn
        self.numpert = numpert
        self.states = uniqStates
        self.xssets = xssetsMtx
        self._lenStates = lenStates

    def setTimeTrace(self, timepoints, *argv):
        """Feed in operational trace.

        The use can feed an operational trace for each of the dependencies
        defined in the problem. If only one dependency exists (e.g. pressure),
        only a pressure trace as a function needs to be provided.
        If more than one dependnecies exist then separate vectors for each
        dependnecy will need to be provided.

        Parameters
        ----------
        timepoints : array
            the time points at which interpolated data will be created.
        argv : arrays
            time-dependent traces for each dependency (e.g. pressure).


        Attributes
        ----------
        xsIntrpSets : dict
            all the interpolated xs sets corresponding to ``timepoints``.
            Keys are the time points and values are xs objects.

        Returns
        -------
        timepoints : array
            the time points at which interpolated data will be created.
        *xsTimeSets : objects
            the cross sections objects corresponding to each time-point

        Raises
        ------
        TypeError
            If ``timepoints`` is not an array.
            If any of the ``argv`` is not an array
        TypeValue
            If the vectors are not of the same size.
            If any of the traces contain negative values.

        Examples
        --------
        >>> timepoints, xsTimeSets = xs.setTimeTrace(
                                     [0, 3.5], [500, 550], [3.0, 3.5])

        """

        timepts = np.array(timepoints)
        npts = len(timepts)  # number of time points

        traces = {}  # time-dependent dependencies
        for idx, arg in enumerate(argv):
            arg = np.array(arg)  # convert to an array if a list is provided
            _isequallength(arg, npts, "Dep<{}>".format(idx))
            _anynegative(arg, "Dep<{}>".format(idx))  # negatives not allowed
            traces[idx] = arg

        # Check that number of traces correspond to the number of dependencies
        _isequallength(traces, self.numdepn, "Number of dependencies")

        # Find required number of points based on the number of dependencies
        if self.numdepn == 1:  # linear interpolation
            self._LinearSplineInterp(traces, timepts)
        elif self.numdepn == 2:  # bilinear
            self._BiLinearSplineInterp(traces, timepts)
        elif self.numdepn == 3:  # trilinear
            self._TriLinearSplineInterp(traces, timepts)
        else:
            raise ValueError("Only up to 3 dependencies are supported and not"
                             " {}".format(self.numdepn))

        self.timepoints = timepts
        return self.timepoints, self.xsTimeSets.values()

    def _LinearSplineInterp(self, traces, timepts):
        """Interpolate the data for a single dependency"""

        if (traces[0] > max(self.states[0])).any() or (
                traces[0] < min(self.states[0])).any():
            raise ValueError("Traces must be within defined range of data "
                             "[{},{}]".format(max(self.states[0]),
                                              min(self.states[0])))
        xsTimeSets = {}
        x = self.states[0]  # create a 1-dim array
        # loop over all the values within the trace and interpolate
        for it, tval in enumerate(timepts):
            xval = traces[0][it]  # time-dependent value of the trace
            idx, = np.where(x == xval)
            if idx.size:  # if the dep. exists in states
                # No need to change the xs-set
                xsTimeSets[it] = self.xssets[idx[0]]
                continue
            else:
                idx = np.where(x > xval)[0][0]

                # Obtain the corresponding xs sets
                xsR = self.xssets[idx]
                xsL = self.xssets[idx-1]
                xsIntp = copy.copy(xsL)  # create a new interp. xs set

            # loop over all the required interpolated attributes
            for attr in INTRP_ATTR:
                vals0 = getattr(xsL, attr)
                vals1 = getattr(xsR, attr)
                vals = LinearInterp(x[idx-1], x[idx], vals0, vals1, xval)
                setattr(xsIntp, attr, vals)  # assign interpolated attribute
            xsTimeSets[it] = xsIntp
        # store all the interpolated xs sets on the object
        self.xsTimeSets = xsTimeSets

    def _BiLinearSplineInterp(self, traces, timepts):
        """Interpolate the data for two dependencies"""

        if (traces[0] > max(self.states[0])).any() or (
                traces[0] < min(self.states[0])).any():
            raise ValueError("Traces must be within defined range of data "
                             "[{},{}]".format(max(self.states[0]),
                                              min(self.states[0])))
        if (traces[1] > max(self.states[1])).any() or (
                traces[1] < min(self.states[1])).any():
            raise ValueError("Traces must be within defined range of data "
                             "[{},{}]".format(max(self.states[1]),
                                              min(self.states[1])))
        xsTimeSets = {}
        x = self.states[0]  # dep. 1
        y = self.states[1]  # dep. 2
        # loop over all the values within the trace and interpolate
        for it, tval in enumerate(timepts):
            # specific time-dependent (xval, yval) for interpolation
            xval = traces[0][it]
            yval = traces[1][it]
            idx, = np.where(x == xval)
            idy, = np.where(y == yval)
            if idx.size and idy.size:  # if the dep. pair exists
                # No need to change the xs-set
                xsTimeSets[it] = self.xssets[idx, idy]
                continue
            else:
                # Find the right and left bounds for x and y dep.
                idx = np.where(x > xval)[0][0]
                idy = np.where(y > yval)[0][0]

                # Obtain the corresponding xs sets
                xs00 = self.xssets[idx-1, idy-1]
                xs10 = self.xssets[idx, idy-1]
                xs01 = self.xssets[idx-1, idy]
                xs11 = self.xssets[idx, idy]
                xsIntp = copy.copy(xs00)  # create a new interp. xs set
            # loop over all the required interpolated attributes
            for attr in INTRP_ATTR:
                vals00 = getattr(xs00, attr)
                vals10 = getattr(xs10, attr)
                vals01 = getattr(xs01, attr)
                vals11 = getattr(xs11, attr)
                # Get and assign interpolated values
                vals = BiLinearInterp(x[idx-1], x[idx], y[idy-1], y[idy],
                                      vals00, vals10, vals01, vals11,
                                      xval, yval)
                setattr(xsIntp, attr, vals)  # assign
            xsTimeSets[it] = xsIntp
        # store all the interpolated xs sets on the object
        self.xsTimeSets = xsTimeSets

    def _TriLinearSplineInterp(self, traces, timepts):
        """Interpolate the data for two dependencies"""
        # function will be completed in the future
        pass

    @staticmethod
    def _sortsets(numdepn, numpert, states, xssets):
        """sorts the data into multidimensional arrays"""
        # store the unique values foe each dependency
        uniqDep = np.empty(numdepn, dtype=np.object)
        lenDep = np.empty(numdepn, dtype=int)

        for idx in range(numdepn):
            uniqDep[idx] = np.unique(states[:, idx])
            lenDep[idx] = len(uniqDep[idx])

        # A matrix to store all the objects/xs sets for all dependencies
        xssetsMtx = np.empty(lenDep, dtype=np.object)

        if numdepn == 1:
            for ix, xval in enumerate(uniqDep[0]):
                cond, = np.where(states[:, 0] == xval)
                if cond.size:
                    xssetsMtx[ix] = xssets[cond[0]]
                else:
                    raise ValueError("Value {} for dep.1 does not exist"
                                     .format(xval))
        if numdepn == 2:
            for ix, xval in enumerate(uniqDep[0]):
                for iy, yval in enumerate(uniqDep[1]):
                    cond, = np.where((states[:, 0] == xval) & (
                        states[:, 1] == yval))
                    if cond.size:
                        xssetsMtx[ix, iy] = xssets[cond[0]]
                    else:
                        raise ValueError("Value {},{} for dep.1-2 do not "
                                         "exist".format(xval, yval))

        if numdepn == 3:
            for ix, xval in enumerate(uniqDep[0]):
                for iy, yval in enumerate(uniqDep[1]):
                    for iz, zval in enumerate(uniqDep[2]):
                        cond, = np.where((states[:, 0] == xval) & (
                            states[:, 1] == yval) & (states[:, 2] == zval))
                        if cond.size:
                            xssetsMtx[ix, iy, iz] = xssets[cond[0]]
                        else:
                            raise ValueError("Value {},{},{} for dep.1-3 do "
                                             "not exist"
                                             .format(xval, yval, zval))

        return uniqDep, lenDep, xssetsMtx


def LinearInterp(x0, x1, vals0, vals1, x):
    """Linear interpolation

    Given the dependnecies x0 (left) and x1 (right) a linear interpolation is
    going to be performed on vectors/floats vals0 and vals1 for a specific x.

    Parameters
    ----------
    x0 : float
        a value of the dependnecy (e.g. temperature)
    x1 : float
        a value of the dependnecy
    vals0 : float or array
        values corresponding to a dependnecy `x0`
    vals1 : float or array
        values corresponding to a dependnecy `x1`
    x : float
        the value of the dependency for the interpolation
    extrpFlag : boolean, default True
        extrapolation flag to indicate whether extrapolation is allowed.

    Raises
    ------
    ValueError
        If x is outside the range [x0, x1] and extrapolation flag is False

    Examples
    --------
    >>> LinearInterp(500., 600., 0.01, 0.02, 550, False)
    ... 0.015
    >>> LinearInterp(500., 600., 0.01, 0.02, 490, True)
    ... 0.009

    """

    _inrange(x, "value of x", [x0, x1])
    xd = (x - x0) / (x1 - x0)
    vals = (1 - xd)*vals0 + xd*vals1

    return vals


def BiLinearInterp(x0, x1, y0, y1, vals00, vals10, vals01, vals11, x, y):
    """Bi-Linear interpolation on a x and y grid

    The grid and values are provided in the schematics below:

        (x0, y0)--------(x1, y0)           vals00 -------- vals10
            |   (x,y)   |                    |     vals       |
            |           |         ==>        |                |
        (x0, y1)-------(x1, y1)            vals01 -------- vals11

     (x, y) are interpolated from knowing x0, x1, y0, and y1, and the weights
     are used to obtain the interpolated `vals` from known values (e.g. vals00)
     that correpond to the grid structure.

    Parameters
    ----------
    x0 : float
        a value of the x-dependnecy (e.g. temperature)
    x1 : float
        a value of the x-dependnecy (e.g. temperature)
    y0 : float
        a value of the y-dependnecy (e.g. pressure)
    y1 : float
        a value of the y-dependnecy (e.g. pressure)
    vals00 : float or array
        values corresponding to dependnecy pair (x0, y0)
    vals10 : float or array
        values corresponding to dependnecy pair (x0, y0)
    vals01 : float or array
        values corresponding to dependnecy pair (x0, y0)
    vals11 : float or array
        values corresponding to dependnecy pair (x0, y0)
    x : float
        the value of the x-dependency for the interpolation
    y : float
        the value of the y-dependency for the interpolation
    extrpFlag : boolean, default True
        extrapolation flag to indicate whether extrapolation is allowed.

    Raises
    ------
    ValueError
        If x is outside the range [x0, x1] or y is outside the range [y0, y1]
        and extrapolation flag is False (i.e., extrapolation is not allowed).

    Examples
    --------
    >>> BiLinearInterp(500., 600., 4, 5, 0.01, 0.02, 0.03, 0.04, 550, 4.5)
    ... 0.025
    >>> BiLinearInterp(500., 600., 4, 5, 0.01, 0.02, 0.03, 0.04, 500, 4)
    ... 0.01
    >>> BiLinearInterp(500., 600., 4, 5, 0.01, 0.02, 0.03, 0.04, 600, 5)
    ... 0.04
    >>> BiLinearInterp(500., 600., 4, 5, 0.01, 0.02, 0.03, 0.04, 450, 3.9)
    ... 0.0029999999999999975

    """

    _inrange(x, "value of x", [x0, x1])
    _inrange(y, "value of y", [y0, y1])

    xd = (x - x0) / (x1 - x0)
    yd = (y - y0) / (y1 - y0)

    # Interpolation over the first x-dependency for different y-values
    vals0 = (1 - xd)*vals00 + xd*vals10
    vals1 = (1 - xd)*vals01 + xd*vals11

    # Interpolation for the second y-dep. based on the x-interpolated values
    vals = (1 - yd)*vals0 + yd*vals1

    return vals


def TriLinearInterp(x0, x1, y0, y1, z0, z1, vals000, vals100, vals010, vals110,
                    vals001, vals101, vals011, vals111, x, y, z):
    """Tri-Linear interpolation on (x, y, z) grid

    The grid and values are provided in the schematics below:


  vals011 *------------* vals111
         /|           /|
        / |          / |
vals001* -----------*  |vals101
       |  |         |  |
       |  |         |  |
   010 |  *---------|--* vals110
       | /          | /
       |/           |/
       * -----------*
     vals000        vals100

     (x, y, z) are interpolated from knowing [x0, x1], [y0, y1], [z1, z2]
     and the weights are used to obtain the interpolated `vals`
     from known values (e.g. vals001, vals111, ...)
     that correpond to the parallelepiped/cubic-based grid structure.

    Parameters
    ----------
    x0 : float
        a value of the x-dependnecy (e.g. temperature)
    x1 : float
        a value of the x-dependnecy (e.g. temperature)
    y0 : float
        a value of the y-dependnecy (e.g. pressure)
    y1 : float
        a value of the y-dependnecy (e.g. pressure)
    z0 : float
        a value of the z-dependnecy (e.g. burnup)
    z1 : float
        a value of the z-dependnecy (e.g. burnup)
    vals000 : float or array
        values corresponding to f(x0, y0, z0)
    vals100 : float or array
        values corresponding to f(x1, y0, z0)
    vals010 : float or array
        values corresponding to f(x0, y1, z0)
    vals110 : float or array
        values corresponding to f(x1, y1, z0)
    vals001 : float or array
        values corresponding to f(x0, y0, z1)
    vals101 : float or array
        values corresponding to f(x1, y0, z1)
    vals011 : float or array
        values corresponding to f(x0, y1, z1)
    vals111 : float or array
        values corresponding to f(x1, y1, z1)
    x : float
        the value of the x-dependency for the interpolation
    y : float
        the value of the y-dependency for the interpolation
    z : float
        the value of the z-dependency for the interpolation
    extrpFlag : boolean, default True
        extrapolation flag to indicate whether extrapolation is allowed.

    Raises
    ------
    ValueError
        If x is outside the range [x0, x1] or y is outside the range [y0, y1]
        and extrapolation flag is False (i.e., extrapolation is not allowed).

    Examples
    --------
    >>> TriLinearInterp(500., 600., 4, 5, 0.01, 0.02, 0.03, 0.04, 550, 4.5)
    ... 0.025
    >>> TriLinearInterp(500, 600, 4, 5, 0, 10, 0.01, 0.02, 0.03, 0.04,
                    0.05, 0.06, 0.07, 0.08, 550, 4.5, 5,
                    extrpFlag=True)
    ... 0.045

    """

    _inrange(x, "value of x", [x0, x1])
    _inrange(y, "value of y", [y0, y1])
    _inrange(z, "value of z", [z0, z1])

    xd = (x - x0) / (x1 - x0)
    yd = (y - y0) / (y1 - y0)
    zd = (z - z0) / (z1 - z0)

    # Interpolation over the first x-dependency for diff y and z values
    vals00 = (1 - xd)*vals000 + xd*vals100
    vals01 = (1 - xd)*vals001 + xd*vals101
    vals10 = (1 - xd)*vals010 + xd*vals110
    vals11 = (1 - xd)*vals011 + xd*vals111

    # Interpolation over the y-dep. for diff. z-values
    vals0 = (1 - yd)*vals00 + yd*vals10
    vals1 = (1 - yd)*vals01 + yd*vals11

    # Interpolation for the z-dep.
    vals = (1 - zd)*vals0 + zd*vals1

    return vals

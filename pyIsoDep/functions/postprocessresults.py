"""postprocessresults

Post processing tool to get specific values and plot results.
Data can be obtain for individual isotopes.

Created on Sat Oct 16 01:00:00 2021 @author: Dan Kotlyar
Last updated on Sat Oct 16 01:30:00 2021 @author: Dan Kotlyar

"""

import numpy as np
import matplotlib.pyplot as plt

from pyIsoDep.functions.checkerrors import _inlist, _isarray, _isstr,\
    _isnumber
from pyIsoDep.functions.header import TIME_UNITS_DICT,\
    TIME_UNITS_LIST, FONT_SIZE


class Results:
    """Container to store results and process them


    Parameters
    ----------
    Nt : 2-dim adday
        Nuclide densities as a function of time


    Attributes
    ----------
    Same as Parameters.


    Examples
    --------
    >>> dep = Results(depscenario)

    """

    def __init__(self, results):
        """reset by copying all the attributes from the results container"""

        self.__dict__ = results.__dict__.copy()

    def getvalues(self, attribute, isotopes=None):
        """Obtain the values of a specific property

        The method obtains the values across all the time-points.

        Parameters
        ----------
        attribute : str
            name of the property/attribute
        isotopes : list or array
            identifier/s of the isotopes

        Returns
        -------
        vals : 1-dim or 2-dim array
            values for the property against all the time-points

        Raises
        ------
        KeyError
            If the isotopes do not exist.
        NameError
            If attribute does not exist.

        Examples
        --------
        >>> dep.getvalues('totalQt')
        ... array([21.82682687, 22.79949867])

        """

        _isstr(attribute, "Attribute")
        if isotopes is not None:
            _isarray(isotopes, "Isotopes Id")
        values = getattr(self, attribute)
        if values.ndim == 2:
            vals, idxFull, idxPart =\
                np.intersect1d(self.fullId, isotopes, assume_unique=True,
                               return_indices=True)
            return values[idxFull, :]
        else:
            return values

    def plot(self, attribute, timeUnits="seconds",
             isotopes=None, xlabel=None, ylabel=None, norm=1,
             fontsize=FONT_SIZE, markers="--*", markerfill=False,
             markersize=6, pltType="linear", newFig=True):
        """plot time-dependent results

        The use of this method is similar to the ``getvalues`` method.
        It is important to note that not all the values have time-dependency.

        Parameters
        ----------
        attribute : str
            name of the property/attribute
        isotopes : list or array
            identifier/s of the nuclides
        xlabel : str
            x-axis label with a default ``Time``
        ylabel : str
            y-axis label with a default ``Output``
        fontsize : float
            font size value
        markers : str or list of strings
            markers type
        markerfill : bool
            True if the marking filling to be included and False otherwise
        markersize : int or float
            size of the marker with a default of 8.

        Raises
        ------
        KeyError
            If the isotopes do not exist.
        NameError
            If attribute does not exist.


        """

        # obtain the values
        values = self.getvalues(attribute, isotopes)

        # check variable types
        _isnumber(fontsize, "Font size")
        _isnumber(markersize, "Marker size")
        _isnumber(norm, "Normalization factor")
        _inlist(timeUnits, "time units", TIME_UNITS_LIST)
        _isstr(pltType, "Plot type")
        if not isinstance(markers, list):
            _isstr(markers, "Markers")
            markers = [markers] * len(values)
        elif isinstance(markers, list):
            markers = markers * len(values)

        if ylabel is not None:
            _isstr(ylabel, "ylabel")
        else:
            ylabel = "Output"

        if xlabel is None:
            xlabel = "Time, " + timeUnits
        _isstr(xlabel, "xlabel")

        timepoints =\
            np.cumsum(np.append(0, self.timesteps))/TIME_UNITS_DICT[timeUnits]

        if markerfill:
            mfc = "white"  # marker fill color
        else:
            mfc = None

        if newFig:
            plt.figure()
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)

        # plot the results for multiple isotopes
        if values.ndim == 2:
            if values.shape[1] != self.nsteps+1:
                raise ValueError("The attribute {} has no time-dependency".
                                 format(attribute))
            values = values / norm
            for idx in range(values.shape[0]):
                if pltType == "linear":
                    plt.plot(timepoints, values[idx, :], markers[idx],
                             mfc=mfc, ms=markersize)
                elif pltType == "semilogx":
                    plt.semilogx(timepoints, values[idx, :], markers[idx],
                                 mfc=mfc, ms=markersize)
                elif pltType == "loglog":
                    plt.loglog(timepoints, values[idx, :], markers[idx],
                               mfc=mfc, ms=markersize)
            plt.legend(isotopes)
        else:
            if attribute == "flux" or attribute == "power":
                timepoints = 0.5*(timepoints[1::] + timepoints[:-1:])
            elif len(values) != self.nsteps+1:
                raise ValueError("The attribute {} has no time-dependency".
                                 format(attribute))
            values = values / norm
            if pltType == "linear":
                plt.plot(timepoints, values, markers[0], mfc=mfc,
                         ms=markersize)
            elif pltType == "semilogx":
                plt.semilogx(timepoints, values, markers[0], mfc=mfc,
                             ms=markersize)
            elif pltType == "loglog":
                plt.loglog(timepoints, values, markers[0], mfc=mfc,
                           ms=markersize)

        plt.grid()
        plt.rc('font', size=fontsize)      # text sizes
        plt.rc('axes', labelsize=fontsize)  # labels
        plt.rc('xtick', labelsize=fontsize)  # tick labels
        plt.rc('ytick', labelsize=fontsize)  # tick labels

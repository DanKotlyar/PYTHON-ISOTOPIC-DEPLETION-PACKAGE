"""postprocessresults

Post processing tool to get specific values and plot results.
Data can be obtain for individual isotopes.

Created on Sat Oct 16 01:00:00 2021 @author: Dan Kotlyar
Last updated on Thrus Oct 28 09:12:00 2021 @author: Matt Krecicki

"""

import numpy as np
import pandas as pa
import h5py
import matplotlib.pyplot as plt


from pyIsoDep.functions.checkerrors import _inlist, _isarray, _isstr,\
    _isnumber, _ispositive, _isbool, _iszeropositive
from pyIsoDep.functions.header import TIME_UNITS_DICT,\
    TIME_UNITS_LIST, FONT_SIZE, HDF5_GROUPS, DATA_ATTR, ZAI_DICT,\
        TIME_UNITS_CONV_MTX
from pyIsoDep.functions.generatedata import TransmutationData


class Results:
    """Container to store results and process them, also reconstructs results
    from hdf5 file


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
    >>> dep = Results("dep.h5")

    """

    def __init__(self, results, includeXS=True):
        """reset by copying all the attributes from the results container"""
        if type(results) is str:
            self.__recover(results, includeXS)
        else:
            self.__dict__ = results.__dict__.copy()


    def __buildGroup(self, f, key, attrs): 
        """function reconstructs a single groups results data from hdf5 file"""
        for i in attrs:
            try:
                data = f[key][i][()]
                if type(data) is bytes:
                    data = str(data, "utf-8")
                setattr(self, i, data)
            except:
                print("{} not found in results".format(i))


    def __buildCrossSectionLibary(self, f, xsKeys):
        """function rebuilds cross section libaries from hdf5 file"""
        xslibs = {}
        for i in list(f["xsData"].keys()):
            xslib = TransmutationData(libraryFlag=False)
            for j in xsKeys:
                data = f["xsData"][i][j][()]
                setattr(xslib, j, data)
            xslibs[float(i)] = xslib
        self._xsDataSets = xslibs
    
    
    def __recover(self, file, includeXS=True):
        """function recovers all results data container from hdf5 file"""
        _isstr(file, "results hdf5 output file name")
        _isbool(includeXS, "flag to include xs libaries")
        keys =  list(HDF5_GROUPS.keys())
        keys.remove("xsData")
        with h5py.File(file, "r+") as f:
            for i in keys:
                try:
                    self.__buildGroup(f, i, HDF5_GROUPS[i])
                except:
                    print("{} not found in results".format(i))
            if includeXS:
                try:
                    self.__buildCrossSectionLibary(f, HDF5_GROUPS["xsData"])
                except:
                    print("XS libaries not found in results")


    def __exportGroup(self, group, ATTR, obj=None):
        """function exports a groups results data to hdf5 file"""
        for i in ATTR:
            try:
                _isstr(i, "Attribute")
                if obj is not None:
                    data = getattr(obj, i)
                else:
                    data = getattr(self, i)
                if type(data) in [np.ndarray, list]:                    
                    if type(data) is list: data = np.asarray(data)     
                    group.create_dataset(i, data=data, dtype=str(data.dtype))
                elif type(data) is str:
                    group.create_dataset(i, data=data.encode("ascii", "ignore"))
                else:
                    group.create_dataset(i, data=data, dtype=type(data))
            except:
                pass
            

    def __exportXsSet(self, name, xslib, group):
        """function exports cross section data set to hdf5 file"""
        subgroup = group.create_group(str(name))
        self.__exportGroup(subgroup, list(DATA_ATTR.keys()), obj=xslib)
        

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
    
    
    def export(self, filename, includeXS=True):
        """function exports results to hdf5 file
        
        The method allows the entire simulation results, cross section
        libaries, and meta-data to be written to an hdf5 file. 
        

        Parameters
        ----------
        filename : str
            hdf5 file output name.

        Returns
        -------
        None.
        
        Examples
        --------
        >>> res.export('results.h5', includeXS=True)
        
        """
        
        _isstr(filename, "results hdf5 output file name")
        _isbool(includeXS, "flag to include cross section libaries")
        with h5py.File(filename, "w") as f:
            keys =  list(HDF5_GROUPS.keys())
            keys.remove("xsData")
            for i in keys:
                try:
                    self.__exportGroup(f.create_group(i), HDF5_GROUPS[i])
                except:
                    pass
            if includeXS:
                xs = f.create_group("xsData")
                for i in list(self._xsDataSets.keys()):               
                    self.__exportXsSet(i, self._xsDataSets[i], xs)


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


    def __id2zai(self, Id):
        """converts Id into string ZAI name"""
        zai = ZAI_DICT[int(np.floor(Id/10000))] + "-" +\
            str(int(np.floor(Id/10)) - 1000*int(np.floor(Id/10000)))
        if Id % 2 != 0: zai = zai + "m"
        
        return zai


    def __convertTimepoint(self, timepoint, timeUnit):
        """function covnerts time point based on desired units"""
        return timepoint * TIME_UNITS_CONV_MTX\
            [TIME_UNITS_LIST.index(self.timeunits)]\
                [TIME_UNITS_LIST.index(timeUnit)]


    def rank(self, parameter="Qt", timepoint=None, timeUnit=None):
        """function ranks parameter of interest from most important to least 
        important. 
        

        Parameters
        ----------
        parameter : str, optional
            output of interest to be ranked. The optional are "Qt",
            "toxicityIngestion", "toxicityInhalation", and At. The default is
            "Qt".
        timepoint : float, optional
            specific time point of interest. 

        Returns
        -------
        df : dict
            dictionary containing ranking of parameter

        """
        if timepoint is not None:
            _iszeropositive(timepoint, "time point of interest")
        if timepoint is not None and timeUnit is not None:
            _inlist(timeUnit, "time units", TIME_UNITS_LIST)
            timepoint = self.__convertTimepoint(timepoint, str(timeUnit))
            
        _inlist(parameter, "rank parameter of interest",\
                ["Qt", "toxicityIngestion", "toxicityInhalation", "At",
                 "reactivity"])

        Ids, data = getattr(self, "fullId"), getattr(self, parameter)
        intgrl, zai = [], []
        
        if timepoint is None:
            for i in range(len(data[:,0])): intgrl.append(np.sum(data[i,:]))
        else:
            intgrl = data[:,np.argmin(abs(self.timepoints - timepoint))].tolist()
        
        #calculate rank
        for i in Ids: zai.append(self.__id2zai(i))
        df = pa.DataFrame(data={"Id": Ids, "ZAI": zai, parameter: intgrl})
        df = df.sort_values(parameter, ascending=False)
        df["cumlative sum"] =\
            np.cumsum(df[parameter].values/np.sum(df[parameter].values))
        
        return df
        
        
"""loaddecaydata

Interface into pre-generated decay data.

Provides a very simple method to access decay data that includes:
    - List of isotopes
    - Atomic weights
    - Decay constants 1/sec
    - Q-values, W/Bq
    - Capture to metastable branching ratios
    - Fission yields matrices (thermal and fast)
    - Decay matrix
    - Ingestion and Inhalation coefficients, Sv/Bq

Data is stored in a hierarchical data format (HDF) file.
Package is designated to be used for the purpose of general depletion and decay
analysis.

Created on Sat Oct 10 05:30:00 2021 @author: Dan Kotlyar
Last updated on Sat Oct 10 05:30:00 2021 @author: Dan Kotlyar

"""

import pathlib

import h5py


class DecayData:
    """Interface into the provided data file

    Parameters
    ----------
    h5path : Union[str, pathlib.Path, h5py.File, h5py.Group]
        Potentially supports passing an opened File object
        or Group directly, but best to just provide the file name

    Methods
    -------
    read : read the HDF5 data set for a specific material
    description : the description for a specific property
    materials: get a list of all the materials in the database

    Raises
    ------
    OSError
        If the ``h5path`` is not valid.

    Examples
    --------
    >>> table = DecayData("decay_data.h5")

    """

    def __init__(self, h5path):
        if isinstance(h5path, (str, pathlib.Path)):
            self._h5 = h5py.File(h5path, "r")
        elif not isinstance(h5path, (h5py.File, h5py.Group)):
            raise TypeError(h5path)
        else:
            self._h5 = h5path

    def properties(self):
        """Obtain all existing properties / data fields

        Returns
        -------
        list
            All the data fields that exist

        Raises
        ------
        KeyError
            If the data file has no properties.

        Examples
        --------
        >>> table = DecayData("decay_data.h5)
        >>> table.properties()
        ['IDlist', 'AW', ...]

        """

        allPtys = self._h5.keys()
        if not allPtys:
            raise KeyError("No properties in the current database")
        return list(allPtys)

    def getvalues(self, pty):
        """Obtain a complete data set for a specific property

        Parameters
        ----------
        pty : string
            name of the property in the databse, e.g. "AW"

        Returns
        -------
        Array (2-dim or 1-dim)
            vals

        Raises
        ------
        KeyError
            If the material ``pty`` does not exist.
        TypeError
            If the material ``pty`` is not a string

        Examples
        --------
        >>> table = DecayData("decay_data.h5)
        >>> ID = table.read("IDlist")

        """

        # Obtain the set for the specific material
        if not isinstance(pty, str):
            raise TypeError("property must be a string and not {}".format(pty))
        ptyset = self._h5.get(pty)
        if ptyset is None:
            raise KeyError("Property {} does not exist".format(pty))
        return ptyset[()]

    def description(self, pty):
        """Obtain the description for a specific property

        Parameters
        ----------
        pty : string
            name of the property in the databse, e.g. "AW"

        Returns
        -------
        str
            Description of the property

        Raises
        ------
        KeyError
            If the property does not exist.
            If the property does not exist.
        TypeError
            If ``pty`` is not str type.


        Examples
        --------
        >>> table = DecayData("decay_data.h5)
        >>> ID = table.description("IDlist")

        """

        # Obtain the set for the specific material
        if not isinstance(pty, str):
            raise TypeError("property must be a string and not {}".format(pty))
        ptyset = self._h5.get(pty)
        if ptyset is None:
            raise KeyError("Property {} does not exist".format(pty))
        description = ptyset.attrs.get("description")
        if description is None:
            raise KeyError("No description for property {}".format(pty))
        return description

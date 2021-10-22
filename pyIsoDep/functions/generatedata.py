"""generatedata

Function that collects the transmutation data and
prepares the transmutation matrix required for depletion or decay
calculations.

The following reactions are considered:
    - Radiative capture
    - Capture to ground state
    - Capture to metastable
    - n, 2n
    - n, 3n
    - Fission
    - n, alpha
    - n, proton
    - n, deuteron
    - n, tritium

Created on Sat Oct 09 11:30:00 2021 @author: Dan Kotlyar
Last updated on Sat Oct 16 00:30:00 2021 @author: Dan Kotlyar

"""

import numpy as np
import numbers

from pyIsoDep import setDataPath
from pyIsoDep.functions.loaddecaydata import DecayData
from pyIsoDep.functions.checkerrors import _exp2dshape, _is1darray,\
    _isequallength, _isnonNegativeArray, _is2darray, _inrange
from pyIsoDep.functions.header import H5_PATH, BARN_2_CM2, JOULE_2MEV,\
    IDX_XS, DATA_ATTR


class TransmutationData:
    """Container to store all the information required to perform depletion

    The dacay data can be read from a pre-generated hdf5 file. The latter
    also includes fission yields and branching ratios.
    Cross sections can be directly read via the ``read`` method.
    If fission yields or branching ratios are provided then they will be used
    to overwrite existing values within the library.

    Parameters
    ----------
    libraryFlag : bool
        A flag to indicate whether a pre-generated library is used.
    h5path : str
        Full directory path to the hdf5 data library file
    wgtFY : float
        fission yield weighting factor between the thermal and fast fission
        yields, i.e. <fy> = wgt*fyThermal + (1-wgt)*fyFast. Provided, the
        pre-generated library is used, otherwise it is a redundant parameter.

    Attributes
    ----------
    fullId : 1-dim array
        Nuclides ZAID number (e.g., 541350)
    nIsotopes : int
        Number of isotopes
    AW : 1-dim array
        Atomic weight
    Q : 1-dim array
        Decay heat coefficient for each isotope in W/Bq
    BR : 1-dim array
        branchin ratio that lead to a isomeric state
    lmbda : 1-dim array
        decay constants in 1/sec
    decaymtx : 2-dim array
        Decay matrix
    ingestion : 1-dim array
        Ingestion coefficients in Sv/Bq
    inhalation : 1-dim array
        Inhalattion coefficients in Sv/Bq
    fymtx : 2-dim array
        fission yields matrix for all the fathers-daughters isotopes

    Returns
    -------
    2-dim array
        Transmutation matrix

    Examples
    --------
    >>> xs = TransmutationData(libraryFlag=True)

    """

    def __init__(self, libraryFlag=True, h5path=None, wgtFY=0.0):
        """reset values with a complete list of all the nuclides"""

        self.libraryFlag = libraryFlag
        # load the decay data
        if libraryFlag:
            if h5path is None:
                h5DefaultPath = setDataPath(H5_PATH)
                datalib = DecayData(h5DefaultPath)
            else:
                datalib = DecayData(h5path)

            _inrange(wgtFY, "Fission yield weight", [0.0, 1.0])
            self.fullId = np.array(datalib.getvalues("IDlist"), dtype=int)
            self.nIsotopes = len(self.fullId)
            self.AW = datalib.getvalues("AW")
            self.Q = datalib.getvalues("Q")
            self.BR = datalib.getvalues("BR")
            self.lmbda = datalib.getvalues("lambda")
            self.decaymtx = datalib.getvalues("decayMatrix")
            self.ingestion = datalib.getvalues("ingestion")
            self.inhalation = datalib.getvalues("inhalation")
            self.fymtx = wgtFY*datalib.getvalues("thermalFY") +\
                (1-wgtFY)*datalib.getvalues("fastFY")
        else:
            self.fullId = None
            self.nIsotopes = None
            self.AW = None
            self.Q = None
            self.BR = None
            self.lmbda = None
            self.decaymtx = None
            self.ingestion = None
            self.inhalation = None
            self.fymtx = None

    def ReadData(self, ID, sig_f, sig_c, sig_c2m=None, sig_n2n=None,
                 sig_n3n=None, sig_alpha=None, sig_p=None, sig_d=None,
                 sig_t=None, fymtx=None, EfissMeV=None, BR=None,
                 decaymtx=None, flagBarns=True):
        """Read data and build the transmutation matrix

        Parameters
        ----------
        ID : 1-dim array
            Partial array of nuclides for which cross sections are provided
        sig_f : 1-dim array
            Fission cross sections in barns
        sig_c : 1-dim array
            Radiative capture cross sections in barns
        sig_c2m : 1-dim array
            Radiative capture cross sections leading to metastable in barns
        sig_n2n : 1-dim array
            n, 2n in barns
        sig_n3n : 1-dim array
            n, 3n in barns
        sig_alpha : 1-dim array
            (n, alpha) in barns
        sig_p : 1-dim array
            (n, proton) in barns
        sig_d : 1-dim array
            (n, deuterium) in barns
        sig_t : 1-dim array
            (n, tritium) in barns
        EfissMeV : 1-dim array
            fission energy in MeV for all the isotopes
        BR : 1-dim array
            Branching ratios that lead to an isomeric state
        fymtx : 2-dim array
            fission yields matrix for all the fathers-daughters isotopes
        decaymtx : 2-dim array
            decay matrix
        flagBarns : bool
            if true cross sections are provided in barns otherwise in cm**2


        Returns
        -------
        2-dim array
            Transmutation matrix

        Raises
        ------
        KeyError
            If any of the arrays is not an array.
        TypeError
            If any of the arrays is not an array.
        ValueError
            If any of the arrays have neagtive values or are not of equal size

        Examples
        --------
        >>> xs = TransmutationData(libraryFlag=True)
        >>> xs.ReadData(ID, sig_f, sig_c)

        """

        # store the 1-gr cross sections in a matrix
        # ---------------------------------------------------------------------
        xsData =\
            self._storexs(ID, sig_f, sig_c, sig_c2m, sig_n2n, sig_n3n,
                          sig_alpha, sig_p, sig_d, sig_t, fymtx, EfissMeV, BR,
                          decaymtx, flagBarns)

        # For each parent define the products for all the possible reactions
        # e.g. 922350 absorbs a neutron and leads to 922360
        # ---------------------------------------------------------------------
        # the starting index corresponds to the radiative capture reaction idx
        idxC = IDX_XS["c"]
        # empty matrix to store products ID's
        prodcutsIDs = np.zeros((self.nIsotopes, len(IDX_XS)-idxC))
        # change ID for isomers to find correct daughters
        parents = np.around(self.fullId, -1)
        prodcutsIDs[:, IDX_XS["c"]-idxC] = parents + 10  # n,c - ground
        prodcutsIDs[:, IDX_XS["c2m"]-idxC] = parents + 11  # n,c - meta
        prodcutsIDs[:, IDX_XS["n2n"]-idxC] = parents - 10  # n,2n
        prodcutsIDs[:, IDX_XS["n3n"]-idxC] = parents - 20  # n,3n
        prodcutsIDs[:, IDX_XS["alpha"]-idxC] = parents - 20030  # alpha
        prodcutsIDs[:, IDX_XS["p"]-idxC] = parents - 10000  # proton
        prodcutsIDs[:, IDX_XS["d"]-idxC] = parents - 10010  # deuteron
        prodcutsIDs[:, IDX_XS["t"]-idxC] = parents - 10020  # tritium
        prodcutsIDs[prodcutsIDs < 100] = 0

        # Create the 1-g transmutation matrix (without fission components)
        # ---------------------------------------------------------------------
        trmtx = np.zeros((self.nIsotopes, self.nIsotopes))
        # fill diagonal with negative absorption values
        np.fill_diagonal(trmtx, -xsData[:, IDX_XS["abs"]])

        for idxrow in range(self.nIsotopes):
            for idxcol, product in enumerate(prodcutsIDs[idxrow, :]):
                trmtx[self.fullId == product, idxrow] =\
                    xsData[idxrow, idxcol+idxC]

        # Create the product between fission yields and fission cross sections
        # ---------------------------------------------------------------------
        fissmtrx = np.tile(xsData[:, IDX_XS["f"]], (self.nIsotopes, 1))
        fissmtrx *= self.fymtx

        # Final transmutation matrix
        # ---------------------------------------------------------------------
        trmtx += fissmtrx
        self.transmutationmtx = trmtx

    def Condense(self, ID, printWarnings=False):
        """Condense all the data to contain only a given list of isotopes

        Parameters
        ----------
        ID : 1-dim array
            List of nuclides
        printWarnings : bool
            Flag to indicate whether warnings should be printed

        Raises
        ------
        ValueError
            If any of the arrays have neagtive values or are not of equal size


        """

        # Check that ID is a non-negative 1-dim array
        # ---------------------------------------------------------------------
        ID = np.array(ID)
        _is1darray(ID, "ID array")
        _isnonNegativeArray(ID, "ID array")

        # intersect between the full and given list of nuclides
        # ---------------------------------------------------------------------
        vals, idxFull, idxPart =\
            np.intersect1d(self.fullId, ID, assume_unique=True,
                           return_indices=True)

        # filter the data according to the provided ID array
        # ---------------------------------------------------------------------
        attrFlag = False
        for attr in DATA_ATTR.keys():
            if hasattr(self, attr):
                attrFlag = True
                val = getattr(self, attr)
                if isinstance(val, numbers.Real):
                    setattr(self, attr, len(idxPart))
                elif val.ndim == 1:
                    setattr(self, attr, val[idxFull])
                elif attr == "xsData":
                    setattr(self, attr, val[idxFull, :])
                else:  # 2-dim array
                    setattr(self, attr, val[:, idxFull][idxFull, :])
            else:
                if printWarnings:
                    print("No attribute <{}> in the container".format(attr))
        if not attrFlag:
            raise ValueError("No attributes at all within the class")

    def _storexs(self, ID, sig_f, sig_c, sigc2m, sig_n2n, sig_n3n,
                 sig_alpha, sig_p, sig_d, sig_t, fymtx, EfissMeV, BR,
                 decaymtx, flagBarns):
        """store all the cross section in a specific format"""

        # get all the cross sections in a single matrix (includes IDs)
        # get the fission energy, fission yields matrix, and decay matrix
        xsDataPart, EfissMeVPart, fymtxPart, decaymtxPart = _checkxs(
            ID, sig_f, sig_c, sigc2m, sig_n2n, sig_n3n, sig_alpha, sig_p,
            sig_d, sig_t, fymtx, EfissMeV, BR, decaymtx, flagBarns)

        # External data library is not provided
        if not self.libraryFlag:
            nIsotopes = len(ID)
            self.fullId = ID
            self.nIsotopes = nIsotopes
            self.xsData = xsDataPart
            self.EfissMeV = EfissMeVPart
            self.EfissJoule = EfissMeVPart / JOULE_2MEV
            self.fymtx = fymtxPart
            self.decaymtx = decaymtxPart
            return xsDataPart

        # intersect between the full and partial list of nuclides
        vals, idxFull, idxPart =\
            np.intersect1d(self.fullId, ID, assume_unique=True,
                           return_indices=True)

        # Build a matrix to store cross sections for the full list of nuclides
        xsData = np.zeros((self.nIsotopes, len(IDX_XS)))
        xsData[idxFull, :] = xsDataPart[idxPart, :]  # populate the data

        # Energy per fission for all the nuclides
        EfissMeV = np.zeros(self.nIsotopes)
        EfissMeV[idxFull] = EfissMeVPart[idxPart]
        EfissJoule = EfissMeV / JOULE_2MEV  # convert MeV to Joules

        # Multiply metastable capture cross section with the branching ratio
        if BR is None:
            xsData[:, IDX_XS["c2m"]] = xsData[:, IDX_XS["c"]] * self.BR
            xsData[:, IDX_XS["c"]] = xsData[:, IDX_XS["c"]] * (1 - self.BR)

        # Overwrite a pre-generated fission matrix
        if fymtxPart is not None:
            fymtx = np.zeros((self.nIsotopes, self.nIsotopes))
            for idx in range(len(idxFull)):
                fymtx[idxFull, idxFull[idx]] = fymtxPart[idxPart, idxPart[idx]]
            self.fymtx = fymtx

        # Overwrite a pre-generated decay matrix
        if decaymtxPart is not None:
            decaymtx = np.zeros((self.nIsotopes, self.nIsotopes))
            for idx in range(len(idxFull)):
                decaymtx[idxFull, idxFull[idx]] =\
                    decaymtxPart[idxPart, idxPart[idx]]
            self.decaymtx = decaymtx

        # Store the matrix with the cross sections
        self.xsData = xsData
        self.EfissMeV = EfissMeV
        self.EfissJoule = EfissJoule

        return xsData


# -----------------------------------------------------------------------------
# Supplementary functions to check errors
# -----------------------------------------------------------------------------
def _checkxs(ID, sig_f, sig_c, sigc2m, sig_n2n, sig_n3n, sig_alpha, sig_p,
             sig_d, sig_t, fymtx, EfissMeV, BR, decaymtx, flagBarns):
    """check that all the cross sections are properly provided"""

    _is1darray(ID, "Nuclides Ids")
    numN = len(np.unique(ID))  # number of unique isotopes
    _isequallength(ID, numN, "Nuclides Ids")
    ID = np.array(ID, dtype=int)
    if sig_c is None:
        sig_c = np.zeros(numN)
    else:
        sig_c = np.array(sig_c)

    if sigc2m is None:
        sigc2m = np.zeros(numN)
    else:
        sigc2m = np.array(sigc2m)

    if sig_n2n is None:
        sig_n2n = np.zeros(numN)
    else:
        sig_n2n = np.array(sig_n2n)

    if sig_n3n is None:
        sig_n3n = np.zeros(numN)
    else:
        sig_n3n = np.array(sig_n3n)

    if sig_f is None:
        sig_f = np.zeros(numN)
    else:
        sig_f = np.array(sig_f)

    if sig_alpha is None:
        sig_alpha = np.zeros(numN)
    else:
        sig_alpha = np.array(sig_alpha)

    if sig_p is None:
        sig_p = np.zeros(numN)
    else:
        sig_p = np.array(sig_p)

    if sig_d is None:
        sig_d = np.zeros(numN)
    else:
        sig_d = np.array(sig_d)

    if sig_t is None:
        sig_t = np.zeros(numN)
    else:
        sig_t = np.array(sig_t)

    if EfissMeV is None:
        EfissMeV = _FissionEnergy(ID)
    else:
        EfissMeV = np.array(EfissMeV, dtype=float)

    if BR is not None:
        BR = np.array(BR)

    # check if all variables are 1-dim arrays
    _is1darray(sig_c, "Capture XS")
    _is1darray(sigc2m, "Capture to metastable XS")
    _is1darray(sig_n2n, "(n, 2n) XS")
    _is1darray(sig_n3n, "(n, 3n) XS")
    _is1darray(sig_f, "(n, fission) XS")
    _is1darray(sig_alpha, "(n, alpha) XS")
    _is1darray(sig_p, "(n, proton) XS")
    _is1darray(sig_d, "(n, deutron) XS")
    _is1darray(sig_t, "(n, tritium) XS")
    _is1darray(EfissMeV, "Fission energy in MeV")
    if BR is not None:
        _is1darray(BR, "Branching ratios")
    if fymtx is not None:
        fymtx = np.array(fymtx)
        _is2darray(fymtx, "Fission yields matrix")
    if decaymtx is not None:
        decaymtx = np.array(decaymtx)
        _is2darray(decaymtx, "Decay matrix")

    # check that all variables are with the same size
    _isequallength(sig_c, numN, "Capture XS")
    _isequallength(sigc2m, numN, "Capture to metastable XS")
    _isequallength(sig_n2n, numN, "(n, 2n) XS")
    _isequallength(sig_n3n, numN, "(n, 3n) XS")
    _isequallength(sig_f, numN, "(n, fission) XS")
    _isequallength(sig_alpha, numN, "(n, alpha) XS")
    _isequallength(sig_p, numN, "(n, proton) XS")
    _isequallength(sig_d, numN, "(n, deutron) XS")
    _isequallength(sig_t, numN, "(n, tritium) XS")
    _isequallength(EfissMeV, numN, "Fission energy in MeV")
    if BR is not None:
        _isequallength(BR, numN, "Branching ratios")
    if fymtx is not None:
        _exp2dshape(fymtx, (numN, numN), "Fission yields matrix")
    if decaymtx is not None:
        _exp2dshape(decaymtx, (numN, numN), "Decay matrix")

    # check if all variables do not contain negative values
    _isnonNegativeArray(sig_c, "Capture XS")
    _isnonNegativeArray(sigc2m, "Capture to metastable XS")
    _isnonNegativeArray(sig_n2n, "(n, 2n) XS")
    _isnonNegativeArray(sig_n3n, "(n, 3n) XS")
    _isnonNegativeArray(sig_f, "(n, fission) XS")
    _isnonNegativeArray(sig_alpha, "(n, alpha) XS")
    _isnonNegativeArray(sig_p, "(n, proton) XS")
    _isnonNegativeArray(sig_d, "(n, deutron) XS")
    _isnonNegativeArray(sig_t, "(n, tritium) XS")
    _isnonNegativeArray(EfissMeV, "Fission energy in MeV")
    if BR is not None:
        _isnonNegativeArray(BR, "Branching ratios")
    if fymtx is not None:
        _isnonNegativeArray(fymtx, "Fission yields matrix")

    # builds xs data matrix to store all the cross sections
    xsData = np.zeros((numN, len(IDX_XS)))
    # Absorption cross section
    sig_abs = sig_c + sigc2m + sig_n2n + sig_n3n + sig_f + sig_alpha + sig_p +\
        sig_d + sig_t

    if flagBarns:
        convUnits = BARN_2_CM2
    else:
        convUnits = 1.0

    xsData[:, IDX_XS["id"]] = ID
    xsData[:, IDX_XS["abs"]] = sig_abs * convUnits
    xsData[:, IDX_XS["c"]] = sig_c * convUnits
    xsData[:, IDX_XS["c2m"]] = sigc2m * convUnits
    xsData[:, IDX_XS["n2n"]] = sig_n2n * convUnits
    xsData[:, IDX_XS["n3n"]] = sig_n3n * convUnits
    xsData[:, IDX_XS["f"]] = sig_f * convUnits
    xsData[:, IDX_XS["alpha"]] = sig_alpha * convUnits
    xsData[:, IDX_XS["p"]] = sig_p * convUnits
    xsData[:, IDX_XS["d"]] = sig_d * convUnits
    xsData[:, IDX_XS["t"]] = sig_t * convUnits

    if BR is not None:
        xsData[:, IDX_XS["c2m"]] = sig_c * BR
        xsData[:, IDX_XS["c"]] = sig_c * (1 - BR)

    return xsData, EfissMeV, fymtx, decaymtx


# Obtain the energy per fission for the defined isotopes
# -------------------------------------------------------------------------
def _FissionEnergy(ID):
    """Get the fission energy for the actinides (expression from ORIGEN)"""

    idxNonFiss = ID < 900000  # indices for all the non-actinides
    Z = np.floor(ID / 10000)  # number of protons
    A = np.floor((ID - Z * 10000) / 10)  # num of protons+neutrons

    energyFissMeV = 1.29927e-3 * (Z**2 * A**0.5) + 33.12
    energyFissMeV[idxNonFiss] = 0.0

    return energyFissMeV

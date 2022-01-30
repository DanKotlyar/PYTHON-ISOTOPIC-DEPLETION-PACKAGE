"""maindepletionsolver

Main depletion function that executes depletion and solves the Bateman
equations.
Depletion can be executed with transmutation and decay components or
be carried out for decay only.
The number of coupled ODEs depend on the user inputs.
The default library uses a 1743 x 1743 matrix.


Created on Sat Oct 16 01:30:00 2021 @author: Dan Kotlyar
Last updated on Thrus Oct 21 08:45:00 2021 @author: Matt Krecicki

"""

import numpy as np
import time
from pyIsoDep.functions.batemansolvers import CramSolver, expmSolver,\
    odeintSolver, adaptiveOdeintSolver
from pyIsoDep.functions.checkerrors import _inlist,\
    _isequallength, _anynegative, _is1darray, _ispositive, _isarray
from pyIsoDep.functions.header import NAVO, BQ_2_CURIE,\
    TIME_UNITS_DICT, TIME_UNITS_LIST, DEPLETION_METHODS, DATA_ATTR,\
    BARN_2_CM2, IDX_XS, DECAY_MUST_ATTR, DECAY_EXPECTED_ATTR,\
    TRANSMUATION_ATTR, DECAY_HEAT_ATTR, RADIOTOXICITY_ATTR, ACTIVITY_ATTR,\
    MASS_ATTR

TO_PCM = 1E+5  # convert reactivity to pcm


class MainDepletion:
    """A class to perform depletion or decay analysis

    The depletion is performed by assuming constant 1-group cross section
    values provided by the user. The user can either provide directly the flux
    for each depletion interval or the power.
    If power is provided then the flux is calculated from knowing the isotopic
    concentration and the fission cross sections.

    Parameters
    ----------
    argv : TransmutationData Objects
        Arguments for all the inputted class object/container that contains
        the data, such as decay constants, IDs, etc.
    timeframes : array
        Time frames for all the stored objects (can be used to interpolate xs)

    Attributes
    ----------
    _xsDataSets : dictionary
        Keys represent the timeframes and values the TransmutationData objects
    _timeframes : array
        Time frames for all the stored objects (can be used to interpolate xs)

    Raises
    ------
    TypeError
        If any of ``argv`` contains no basic attributes.

    Examples
    --------
    >>> dep = MainDepletion(0.0, data)

    """

    def __init__(self, timeframes=0.0, *argv):
        """reset values with a complete list of all the nuclides"""

        self._xsintrp = None

        # empty dictionary to store all the cross section and decay data sets
        xsDataSets = {}

        # Check for potential errors
        # ---------------------------------------------------------------------
        timeframes = np.array(timeframes, ndmin=1)  # covert to 1-dim array
        nxssets = len(argv)  # number of xs sets / time frames
        _isequallength(timeframes, nxssets, "Time Frames")
        _anynegative(timeframes, "Time Frames")

        # Check that each data set contains the required attributes
        # ---------------------------------------------------------------------
        timeIdx = 0  # time frame index
        for timeframe, data in zip(timeframes, argv):
            for attr in DATA_ATTR.keys():
                if not hasattr(data, attr):
                    print("No attribute <{}> in data".format(attr))
                elif attr in DECAY_EXPECTED_ATTR and timeIdx == 0:
                    setattr(self, attr, getattr(data, attr))
                elif attr == "nu" and timeIdx == 0:
                    setattr(self, attr, getattr(data, attr))
            for attr in DECAY_MUST_ATTR:
                if not hasattr(self, attr):
                    raise ValueError("No attribute <{}> in data".format(attr))
            timeIdx += 1
            xsDataSets[timeframe] = data

        # Store the data on muted dictionary
        # ---------------------------------------------------------------------
        self._xsDataSets = xsDataSets
        self._timeframes = timeframes
        self._solveTime = None

        # Verify that same IDs and chains are used for all argv (i.e.data sets)
        # ---------------------------------------------------------------------
        if nxssets > 1:
            refId = self.fullId
            for timeframe in timeframes[1::]:
                Id0 = self._xsDataSets[timeframe].fullId
                if (refId != Id0).any():
                    raise ValueError(
                        "fullId for frame {} is not identical to the one in "
                        "frame {}.".format(timeframes[0], timeframe))

    def SetDepScenario(self, power=None, flux=None, timeUnits="seconds",
                       timesteps=None, timepoints=None):
        """Define a depletion or decay scenario

        Parameters
        ----------
        power : array
            Absolute power for each depletion step in Watts
        flux : array
            Absolute flux for each depletion step in n/cm**2/s
        timesteps : array
            Depletion time-steps/intervals in sec, min, hr, days
        timepoints : array
            Depletion time-points in sec, min, hr, days
        timeUnits : string
            Time units {"seconds", "minutes", "hours", "days"}

        Attributes
        ----------
        power : array
            Absolute power for each depletion step in Watts
        flux : array
            Absolute flux for each depletion step in n/cm**2/s
        timesteps : array
            Depletion time-steps/intervals in sec, min, hr, days
        timepoints : array
            Depletion time-points in sec, min, hr, days
        timeUnits : string
            Time units {"seconds", "minutes", "hours", "days"}
        flagPower : bool
            To indicate if power was provided

        Raises
        ------
        ValueError
            If any of the arrays have neagtive values or are not of equal size
            If timesteps and timepoints are not provided
            If timeUnits is not defined in TIME_UNITS_LIST

        Examples
        --------
        >>> dep = MainDepletion(0.0, data)
        >>> dep.SetDepScenario(power=[1, 2], flux=None, timeUnits="hours",
        >>>                    timesteps=[0.5, 0.1], timepoints=None)

        """

        # check for possible errors
        # ---------------------------------------------------------------------
        _inlist(timeUnits, "time units", TIME_UNITS_LIST)
        if timesteps is not None:
            timepoints = np.append(0, timesteps)
            timepoints = np.cumsum(timepoints)
        elif timepoints is not None:
            timesteps = timepoints[1::] - timepoints[0:-1:]
        else:
            raise ValueError("timesteps or timepoints must be provided")

        timesteps = np.array(timesteps)
        _anynegative(timesteps, "timesteps")
        nsteps = len(timesteps)  # number of depletion steps

        if power is not None:
            _isarray(power, "power")
            power = np.array(power)
            self.flagPower = True
            _isequallength(power, nsteps, "power")
            _anynegative(power, "power")
            flux = np.zeros(nsteps)
        elif flux is not None:
            _isarray(flux, "flux")
            flux = np.array(flux)
            _isequallength(flux, nsteps, "flux")
            _anynegative(flux, "flux")
            power = np.zeros(nsteps)
            self.flagPower = False

        # Store the depletion parameters
        # ---------------------------------------------------------------------
        self.usertimesteps = timesteps
        self.timesteps = timesteps * TIME_UNITS_DICT[timeUnits]
        self.timepoints = timepoints
        self.timeunits = timeUnits
        self.nsteps = int(len(timesteps))
        self.power = power
        self.flux = flux

    def SetInitialComposition(self, ID, N0, vol=1.0):
        """Set initial composition

        Parameters
        ----------
        ID : array
            Identification of isotopes following ZZAAA0/1 format
        ND0 : array
            Isotopic concentrations #/cm/b
        vol : float
            Volume of the system in cm**3. Volume is needed if total values
            are required.

        Attributes
        ----------
        providedID : array
            Identification of isotopes following ZZAAA0/1 format
        providedN0 : array
            Isotopic concentrations #/cm/b

        Raises
        ------
        ValueError
            If ID contains negative values
            If power and flux are not provided
            If timesteps and timepoints are not provided
            If timeUnits is not defined in TIME_UNITS_LIST

        Examples
        --------
        >>> dep = MainDepletion(h5path=datapath, fyWgt=0.5)
        >>> dep.SetInitialComposition([541350, 922350], [0.0, 0.021])

        """

        # Check potential errors
        # ---------------------------------------------------------------------
        ID = np.array(ID, dtype=int)
        _is1darray(ID, "Isotopic IDs")
        _anynegative(ID, "Isotopic IDs")
        uniqIsot = len(np.unique(ID))
        _isequallength(ID, uniqIsot, "ID contains identical isotopes")
        N0 = np.array(N0, dtype="float64")
        _is1darray(N0, "Isotopic Concentrations")
        _anynegative(N0, "Isotopic Concentrations")
        _ispositive(vol, "Volume")

        # Remap provided isotopes to the full list pre-generated in datafile
        # ---------------------------------------------------------------------
        # intersect between the full and provided list of IDs
        vals, idxFull, idxPart =\
            np.intersect1d(self.fullId, ID, assume_unique=True,
                           return_indices=True)
        # Define an initial vector with concentrations
        self.N0 = np.zeros(self.nIsotopes)
        self.N0[idxFull] = N0[idxPart]
        self.providedN0 = N0
        self.providedID = ID
        self.volume = vol

    def SolveDepletion(self, method="cram", xsinterp=False, rtol=1E-10):
        """Solve the Bateman equations that include transmutation and decay

        Parameters
        ----------
        method : str
            Method used to solve the Bateman equations
        xsinterp : bool
            Flag to indicate whether interpolation in between timesteps is
            allowed to be performed for the transmutation data.
        rtol : float, optional
            Isotopic concentration convergence criteria, relative difference.
            The default is 1E-10

        Attributes
        ----------
        Nt : 2-dim array
            Concentrations for all the isotopes as a function of time

        Raises
        ------
        ValueError
            If ``method`` is not defined.
            If any of the attribures in ``TRANSMUATION_ATTR`` do not exist.

        Examples
        --------
        >>> dep.SetInitialComposition([541350, 922350], [0.0, 0.021])
        >>> dep.SolveDepletion("cram")
        """

        # Check potential errors
        # ---------------------------------------------------------------------
        _inlist(method, "Method to solve Bateman eqs", DEPLETION_METHODS)
        if method == "cram":
            singleDepletion = CramSolver()
        elif method == "expm":
            singleDepletion = expmSolver()
        elif method == "odeint":
            singleDepletion = odeintSolver(rtol=rtol)

        if self.power is None and self.flux is None:
            raise ValueError("Either power or flux must be defined when "
                             "SetDepScenario is set.")

        # All data sets must contain the required attributes
        for key, data in self._xsDataSets.items():
            for attr in TRANSMUATION_ATTR:
                if not hasattr(data, attr):
                    raise ValueError("No attribute <{}> in data for time={}"
                                     .format(attr, key))

        # Nt will store the concentrations as a function of time
        self.Nt = np.zeros((self.nIsotopes, self.nsteps + 1))
        self.Nt[:, 0] = self.N0  # initial concentrations
        # Nt will store all the weighted cross sections as a function of time
        self.XS = np.zeros((self.nIsotopes, len(IDX_XS), self.nsteps + 1))

        tic = time.perf_counter()  # start timer

        if method == "adaptive":  # if adaptive time mesh required
            adptDepletion = adaptiveOdeintSolver(self, xsinterp, rtol=rtol)
            adptDepletion.solve()
            self = adptDepletion.dep
            toc = time.perf_counter()
            self._solveTime = toc - tic
            return

        for idx, dt in enumerate(self.timesteps):

            # Obtain the interpolated fission energy, xs, and transmutation mtx
            # -----------------------------------------------------------------
            fissE, sigf, transmutationmtx, xsTable =\
                self._getInterpXS(self.timepoints[idx], xsinterp)

            # Store the weighted cross sections:
            # -----------------------------------------------------------------
            self.XS[:, :, idx] = xsTable

            # flux is used directly
            # -----------------------------------------------------------------
            if not self.flagPower:
                # calculate power for this step
                self.power[idx] = (self.flux[idx] * sigf * self.Nt[:, idx] *
                                   fissE * self.volume).sum()

            # power is provided and needs to be converted to flux
            # -----------------------------------------------------------------
            else:
                self.flux[idx] = self.power[idx] / (
                        sigf * self.Nt[:, idx] * fissE * self.volume).sum()

            # define the overall matrix to represent Bateman equations
            # -----------------------------------------------------------------
            mtxA = transmutationmtx*self.flux[idx] + self.decaymtx

            # solve and obtain the concentrations after a single depletion
            # -----------------------------------------------------------------
            self.Nt[:, idx+1] =\
                singleDepletion.solve(mtxA, self.Nt[:, idx], dt)

        toc = time.perf_counter()
        self._solveTime = toc - tic
        self._xsintrp = xsinterp

    def SolveDecay(self, method="cram", rtol=1E-10):
        """Solve the Bateman equations with only the decay chains

        Parameters
        ----------
        method : str
            Method used to solve the decay chains
        rtol : float, optional
            Isotopic concentration convergence criteria, relative difference.
            The default is 1E-10

        Attributes
        ----------
        Nt : 2-dim array
            Concentrations for all the isotopes as a function of time

        Raises
        ------
        ValueError
            If ``method`` is not defined.

        Examples
        --------
        >>> dep.SetInitialComposition([541350, 922350], [0.0, 0.021])
        >>> dep.SolveDepletion("expm")
        """

        # Check potential errors
        # ---------------------------------------------------------------------
        _inlist(method, "Method to solve Bateman eqs", DEPLETION_METHODS)
        if method == "cram":
            singleDepletion = CramSolver()
        elif method == "expm":
            singleDepletion = expmSolver()
        elif method == "odeint":
            singleDepletion = odeintSolver(rtol=rtol)

        # Nt will store the concentrations as a function of time
        self.Nt = np.zeros((self.nIsotopes, self.nsteps + 1))
        self.Nt[:, 0] = self.N0  # initial concentrations
        tic = time.perf_counter()

        if method == "adaptive":  # if adaptive time mesh required
            adptDepletion = adaptiveOdeintSolver(self, False, rtol=rtol)
            adptDepletion.solve()
            self = adptDepletion.dep
            toc = time.perf_counter()
            self._solveTime = toc - tic
            return

        for idx, dt in enumerate(self.timesteps):

            # define the overall matrix to represent Bateman equations
            # -----------------------------------------------------------------
            mtxA = self.decaymtx

            # solve and obtain the concentrations after a single depletion
            # -----------------------------------------------------------------
            self.Nt[:, idx+1] =\
                singleDepletion.solve(mtxA, self.Nt[:, idx], dt)
        toc = time.perf_counter()
        self._solveTime = toc - tic

    def _getInterpXS(self, currtime, interpFlag):
        """Obtains the transmutation data required to solve depletion"""

        timeframes = self._timeframes
        # Find the index of closest data with a time below the current time
        if (currtime <= timeframes).all():
            idx0 = 0
            idx1 = 0
        # Find the index of closest data with a time above the current time
        elif (currtime >= timeframes).all():
            idx0 = len(timeframes) - 1
            idx1 = len(timeframes) - 1
        # Current time is in-between exisiting time frames
        else:
            idx0 = np.where(timeframes <= currtime)[0][-1]
            idx1 = np.where(timeframes >= currtime)[0][0]

        # Cross sections are not interpolated
        if not interpFlag or idx0 == idx1:
            data = self._xsDataSets[timeframes[idx0]]
            fissE = data.EfissJoule  # fission energy in joules
            sigf = data.xsData[:, IDX_XS["f"]] / BARN_2_CM2  # fission xs barns
            transmutationmtx = data.transmutationmtx
            xsTable = data.xsData
        else:
            wgt =\
                (currtime-timeframes[idx0])/(timeframes[idx1]-timeframes[idx0])
            data0 = self._xsDataSets[timeframes[idx0]]
            data1 = self._xsDataSets[timeframes[idx1]]
            fissE0 = data0.EfissJoule
            fissE1 = data1.EfissJoule
            sigf0 = data0.xsData[:, IDX_XS["f"]] / BARN_2_CM2
            sigf1 = data1.xsData[:, IDX_XS["f"]] / BARN_2_CM2
            transmutationmtx0 = data0.transmutationmtx
            transmutationmtx1 = data1.transmutationmtx
            # Weight cross sections
            sigf = (1-wgt)*sigf0 + wgt*sigf1
            fissE = (1-wgt)*fissE0 + wgt*fissE1
            transmutationmtx =\
                (1-wgt)*transmutationmtx0 + wgt*transmutationmtx1
            # weighted cross sections
            xsTable = (1-wgt)*data0.xsData + wgt*data0.xsData
            xsTable[:, 0] = data0.xsData[:, 0]

        return fissE, sigf, transmutationmtx, xsTable

    def DecayHeat(self):
        """Calculate decay heat in Watts

        Attributes
        ----------
        At : 2-dim array
            Activity in Bq for all the isotopes as a function of time
        Qt : 2-dim array
            Decay heat in Watts for all the isotopes as a function of time
        totalQt : 1-dim array
            Total decay heat in Watts as a function of time
        """

        # check that all the attributes exist
        for attr in DECAY_HEAT_ATTR:
            if not hasattr(self, attr):
                raise ValueError("No attribute <{}> in data".format(attr))

        mtxQ = np.tile(self.Q, (self.nsteps + 1, 1))  # Q-values [Watts/Bq]
        mtxL = np.tile(self.lmbda, (self.nsteps + 1, 1))  # Decay constants
        mtxQ = mtxQ.transpose()
        mtxL = mtxL.transpose()
        # Activity, becquerel
        self.At = self.volume * mtxL * self.Nt / BARN_2_CM2
        # Calculate Decay heat
        self.Qt = self.At * mtxQ
        self.totalQt = self.Qt.sum(axis=0)

    def Radiotoxicity(self):
        """Calculate radiotoxicity in Sv

        Attributes
        ----------
        At : 2-dim array
            Activity in Bq for all the isotopes as a function of time
        toxicityIngestion : 2-dim array
            Ingestion radiotoxicity in Sv for all the isotopes against time
        inhalationToxicity : 2-dim array
            Inhalation radiotoxicity in Sv for all the isotopes against time
        totalToxIngestion : 2-dim array
            Total ingestion radiotoxicity in Sv against time
        totalToxInhalation : 2-dim array
            Total inhalation radiotoxicity in Sv against time
        """

        # check that all the attributes exist
        for attr in RADIOTOXICITY_ATTR:
            if not hasattr(self, attr):
                raise ValueError("No attribute <{}> in data".format(attr))

        # Ingestion coefficients [Sv/Bq]
        mtxToxIng = np.tile(self.ingestion, (self.nsteps + 1, 1))
        # Inhalation coefficients [Sv/Bq]
        mtxToxInh = np.tile(self.inhalation, (self.nsteps + 1, 1))
        mtxToxIng = mtxToxIng.transpose()
        mtxToxInh = mtxToxInh.transpose()
        mtxL = np.tile(self.lmbda, (self.nsteps + 1, 1))  # Decay constants
        mtxL = mtxL.transpose()
        # Activity, becquerel
        self.At = self.volume * mtxL * self.Nt / BARN_2_CM2
        # Calculate Radotoxicity
        self.toxicityIngestion = self.At * mtxToxIng
        self.toxicityInhalation = self.At * mtxToxInh
        self.totalToxIngestion = self.toxicityIngestion.sum(axis=0)
        self.totalToxInhalation = self.toxicityInhalation.sum(axis=0)

    def Activity(self):
        """Calculate isotopic and total actitvity in Cuire

        Attributes
        ----------
        At : 2-dim array
            Activity in Bq for all the isotopes as a function of time
        AtCurie : 2-dim array
            Activity in Curie for all the isotopes as a function of time
        totalAtCurie : 1-dim array
            Total acitivity in Curie as a function of time
        """

        # check that all the attributes exist
        for attr in ACTIVITY_ATTR:
            if not hasattr(self, attr):
                raise ValueError("No attribute <{}> in data".format(attr))

        mtxL = np.tile(self.lmbda, (self.nsteps + 1, 1))  # Decay constants
        mtxL = mtxL.transpose()
        # Activity, becquerel
        self.AtCurie = self.volume * mtxL * self.Nt / BARN_2_CM2 / BQ_2_CURIE
        self.totalAtCurie = self.AtCurie.sum(axis=0)

    def Mass(self):
        """Calculate isotopic and total masses in grams

        Attributes
        ----------
        mass : 2-dim array
            Mass in grams for all the isotopes as a function of time
        totalmass : 1-dim array
            Total mass in grams as a function of time
        """
        # check that all the attributes exist
        for attr in MASS_ATTR:
            if not hasattr(self, attr):
                raise ValueError("No attribute <{}> in data".format(attr))

        mtxAW = np.tile(self.AW, (self.nsteps + 1, 1))  # Atomic Weight
        mtxAW = mtxAW.transpose()
        # Mass in grams
        self.massgr = self.volume * mtxAW * self.Nt / NAVO
        self.totalMassgr = self.massgr.sum(axis=0)

    def Reactivity(self, nonLeakageP=1.0):
        """Isotopic reactivity worth using first order perturbation theory.

        Attributes
        ----------
        reactivityworth : 2-dim array
            Reactivity worth in pcm for all the isotopes as a function of time
        reactivity : 1-dim array
            Total reactivity in pcm

        """

        abs_xs = np.zeros((self.nIsotopes, self.nsteps))
        fiss_xs = np.zeros((self.nIsotopes, self.nsteps))
        self.dRho = np.zeros((self.nIsotopes, self.nsteps))
        self.dRhoToRho = np.zeros((self.nIsotopes, self.nsteps))
        self.Rho = np.zeros(self.nsteps)
        self.keff = np.zeros(self.nsteps)

        if self.nu is None:
            raise ValueError("Attribute nu does not exist. Please define in "
                             "Transmutation Data.")

        for i in range(self.nsteps):
            abs_xs[:, i] = self.XS[:, IDX_XS['abs'], i]
            fiss_xs[:, i] = self.XS[:, IDX_XS['f'], i]

            # Multiplication factor at the specific time point
            self.keff[i] = nonLeakageP *\
                (self.Nt[:, i] * self.nu * fiss_xs[:, i]).sum() /\
                (self.Nt[:, i] * abs_xs[:, i]).sum()

            # Total reactivity at the specific time point
            self.Rho[i] = TO_PCM * (1 - 1/self.keff[i])

            # Multiplication factor without the isotope j
            keff_wIsot = np.zeros(self.nIsotopes)
            for j in range(self.nIsotopes):
                keff_wIsot[j] =\
                    nonLeakageP *\
                    ((self.Nt[:, i] * self.nu * fiss_xs[:, i]).sum() -
                     self.Nt[j, i] * self.nu[j] * fiss_xs[j, i]) /\
                    ((self.Nt[:, i] * abs_xs[:, i]).sum() -
                     self.Nt[j, i] * abs_xs[j, i])
            # reactivity without the isotope j
            rho_wIsot = TO_PCM * (1 - 1/keff_wIsot)

            # reactivity worth
            self.dRho[:, i] = self.Rho[i] - rho_wIsot
            self.dRhoToRho[:, i] = self.dRho[:, i] / self.Rho[i]

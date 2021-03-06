"""header

File to define basic parameters shared within the package.
Defines:
    - unit conversions
    - cross section types
    - priniting options
    - required attributes for different containers within the package

@authors: Dan Kotlyar & Matt Krecicki
Created on Sat Oct 16 01:30:00 2021 @author: Dan Kotlyar
Last updated on Sun Feb 26 06:30:00 2021 @author: Dan Kotlyar

List changes or additions:
--------------------------
"Concise description" - MM/DD/YYY - Name initials
Expand the cross section dictionary ``IDX_XS`` - 02/26/2022 - DK

"""

from pyIsoDep.functions.checkerrors import _compare2lists

# -----------------------------------------------------------------------------
#                      CONSTANTS
# -----------------------------------------------------------------------------

NAVO = 0.602214199             # avogadro number
JOULE_2MEV = 6.241507649e+12   # conversion from Joule to MeV
BQ_2_CURIE = 3.7E+10           # Conversion from Bq to Curie
BARN_2_CM2 = 1E-24             # Conversion from barns to cm**2

TIME_UNITS_DICT = {"seconds": 1.0, "minutes": 60.0, "hours": 3600.0,
                   "days": 86400.0}
TIME_UNITS_LIST = list(TIME_UNITS_DICT)

TIME_UNITS_CONV_MTX = [[1.0,       60.0,     3600.0, 86400.0],
                       [1/60.0,    1.0,      60.0,   1440.0],
                       [1/3600.0,  1/60.0,   1.0,    24.0],
                       [1/86400.0, 1/1440.0, 1/24.0, 1.0]]

# -----------------------------------------------------------------------------
#                      DEFAULT OPTIONS
# -----------------------------------------------------------------------------

# current depletion options
DEPLETION_METHODS = ["cram", "expm", "odeint", "adaptive"]
H5_PATH = "bgcore_data.h5"              # Pre-generated librray

# -----------------------------------------------------------------------------
#                      DATA MANAGEMENT
# -----------------------------------------------------------------------------
# Storing indices
IDX_XS = {"id": 0,       # ZAID of the isotope
          "abs": 1,      # absorption
          "f": 2,        # fission (MT=18/19)
          "c": 3,        # radiative capture  (MT=102)
          "c2m": 4,      # capture that leads to metastable
          "n2n": 5,      # n, 2n (MT=16)
          "n3n": 6,      # n, 3n (MT=17)
          "a": 7,        # n, alpha (MT=107)
          "p": 8,        # n, proton  (MT=103)
          "d": 9,        # n, deutron (MT=104)
          "t": 10,       # n, tritium (MT=105)
          "na": 11,      # n, n+alpha (MT=22)
          "np": 12,      # n, n+proton (MT=28)
          "2a": 13,      # n, alpha+alpha (MT=108)
          "t2a": 14,     # n, tritium+alpha+alpha (MT=113)
          "a2n": 15,     # n, alpha+2n (MT=24)
          }

XS_LIST = ["sig_"+xstype for xstype in IDX_XS.keys() if
           (xstype != "id" and xstype != "abs")]
# Check that indices are ordered properly
exp_idxList = [i for i in range(len(IDX_XS))]
act_idxList = list(IDX_XS.values())
_compare2lists(act_idxList, exp_idxList, "Actual idx for cross sections",
               "Expected indices")
# -----------------------------------------------------------------------------
#                      PLOTTING
# -----------------------------------------------------------------------------
FONT_SIZE = 16

# -----------------------------------------------------------------------------
#                      DATA ATTRIBUTES
# -----------------------------------------------------------------------------
DATA_ATTR = {
    "fullId": "Full ZAID list",
    "nIsotopes": "Number of isotopes",
    "AW": "Atomic weight",
    "Q": "Decay heat in W/Bq",
    "BR": "Branching ratios leading to isomeric",
    "lmbda": "Decay constant in 1/s",
    "decaymtx": "Decay matrix",
    "ingestion": "Ingestion coefficients in Sv/Bq",
    "inhalation": "Inhalation coefficients in Sv/Bq",
    "fymtx": "Fission yields matrix",
    "EfissMeV": "Energy per fission in MeV",
    "EfissJoule": "Energy per fission in Joules",
    "xsData": "Cross section data matrix",
    "transmutationmtx": "Transmutation matrix without decay",
    "nu": "Number of neutrons emitted per fission"
    }

# Attributes that must exist for transmutation and decay calculations
TRANSMUATION_ATTR = ["fullId",  "nIsotopes", "decaymtx", "fymtx",
                     "EfissJoule", "xsData", "transmutationmtx"]

# Attributes that must exist to calculate acitivtiy
ACTIVITY_ATTR = ["Nt", "lmbda", "volume"]

# Attributes that must exist to calculate nuclide density in decay calculations
DECAY_MUST_ATTR = ["fullId",  "nIsotopes", "decaymtx"]

# Attributes that must exist to calculate decay heat
DECAY_HEAT_ATTR = ["Nt", "lmbda", "Q", "volume"]

# Attributes that are expected to perform decay calculations
DECAY_EXPECTED_ATTR = ["fullId", "nIsotopes", "AW", "Q", "BR", "lmbda",
                       "decaymtx", "ingestion", "inhalation"]

# Attributes that must exist to calculate radiotoxicity
RADIOTOXICITY_ATTR = ["Nt", "lmbda", "inhalation", "ingestion", "volume"]

# Attributes that must exist to calculate mass
MASS_ATTR = ["Nt", "AW", "volume"]

# Data required for interpolation
INTRP_ATTR = ["fymtx", "EfissJoule", "xsData", "transmutationmtx"]

# Attributes that should not be weighted
NOT_WEIGHT_ATTR = ['_xsintrp', 'fullId', 'nIsotopes', 'AW', 'Q', 'BR', 'lmbda',
                   'decaymtx', 'nu', '_xsDataSets',
                   '_timeframes', '_solveTime', 'flagPower', 'usertimesteps',
                   'timesteps', 'timepoints', 'timeunits', 'nsteps', 'volume']

# hdf5 output file atrributes list
HDF5_GROUPS = {"metaData": ["nIsotopes", "AW", "fullId", "timepoints",
                            "timesteps", "timeunits", "usertimesteps",
                            "flagPower", "_timeframes", "providedID",
                            "nsteps"],
               "results": ["Nt", "At", "AtCurie", "Qt", "flux", "power",
                           "massgr", "totalAtCurie", "totalMassgr", "totalQt",
                           "totalToxIngestion", "totalToxInhalation",
                           "toxicityIngestion", "toxicityInhalation",
                           "Rho", "dRho", "dRhoToRho"],
               "xsData": ["fullId", "nIsotopes", "AW", "Q", "BR", "lmbda",
                          "decaymtx", "ingestion", "inhalation", "fymtx",
                          "EfissMeV", "EfissJoule", "xsData",
                          "transmutationmtx"],
               "deplData": ["BR", "decaymtx", "ingestion", "inhalation",
                            "lmbda", ],
               "inital": ["volume", "N0", "providedN0"]}


INTERVAL_ATTRIBUTES = ["flux", "power", "keff", "Rho", "dRho", "dRhoToRho"]
RANK_PARAMETERS = ["Qt", "toxicityIngestion", "toxicityInhalation", "At",
                   "dRho", "dRhoToRho"]

RANK_ATTRIBUTES = ["Qt", "toxicityIngestion", "toxicityInhalation", "At",
                   "dRho", "dRhoToRho", "AW", "fullId", "Nt", "AtCurie",
                   "massgr", "Q", "BR", "lmbda", "EfissMeV", "EfissJoule"]

# -----------------------------------------------------------------------------
#                      Isotope name conversion dictionary
# -----------------------------------------------------------------------------

ZAI_DICT = {1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O",
            9: "F", 10: "Ne", 11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P",
            16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca", 21: "Sc", 22: "Ti",
            23: "V", 24: "Cr", 25: "Mn", 26: "Fe", 27: "Co", 28: "Ni",
            29: "Cu", 30: "Zn", 31: "Ga", 32: "Ge", 33: "As", 34: "Se",
            35: "Br", 36: "Kr", 37: "Rb", 38: "Sr", 39: "Y", 40: "Zr",
            41: "Nb", 42: "Mo", 43: "Tc", 44: "Ru", 45: "Rh", 46: "Pd",
            47: "Ag", 48: "Cd", 49: "In", 50: "Sn", 51: "Sb", 52: "Te",
            53: "I", 54: "Xe", 55: "Cs", 56: "Ba",  57: "La", 58: "Ce",
            59: "Pr", 60: "Nd", 61: "Pm", 62: "Sm", 63: "Eu", 64: "Gd",
            65: "Tb", 66: "Dy", 67: "Ho", 68: "Er", 69: "Tm", 70: "Yb",
            71: "Lu", 72: "Hf", 73: "Ta", 74: "W", 75: "Re", 76: "Os",
            77: "Ir", 78: "Pt", 79: "Au", 80: "Hg", 81: "Tl", 82: "Pb",
            83: "Bi", 84: "Po", 85: "At", 86: "Rn", 87: "Fr", 88: "Ra",
            89: "Ac", 90: "Th", 91: "Pa", 92: "U", 93: "Np", 94: "Pu",
            95: "Am", 96: "Cm", 97: "Bk", 98: "Cf", 99: "Es", 100: "Fm",
            101: "Md"}

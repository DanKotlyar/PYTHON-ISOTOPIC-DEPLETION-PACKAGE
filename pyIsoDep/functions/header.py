"""header

File to define basic parameters shared within the package.


Created on Sat Oct 16 01:30:00 2021 @author: Dan Kotlyar
Last updated on Sat Oct 16 08:42:00 2021 @author: Matt Krecicki

"""


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


# -----------------------------------------------------------------------------
#                      DEFAULT OPTIONS
# -----------------------------------------------------------------------------

DEPLETION_METHODS = ["cram", "expm", "odeint", "adaptive"]  # current depletion options
H5_PATH = "bgcore_data.h5"              # Pre-generated librray

# -----------------------------------------------------------------------------
#                      DATA MANAGEMENT
# -----------------------------------------------------------------------------
# Storing indices
IDX_XS = {"id": 0,     # ZAID of the isotope
          "abs": 1,    # absorption
          "f": 2,      # fission
          "c": 3,      # radiative capture
          "c2m": 4,    # capture that leads to metastable
          "n2n": 5,    # n, 2n
          "n3n": 6,    # n, 3n
          "alpha": 7,  # n, alpha
          "p": 8,      # n, proton
          "d": 9,      # n, deutron
          "t": 10}     # n, tritium

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

#hdf5 output file atrributes list

HDF5_GROUPS = {"metaData": ["nIsotopes", "AW", "fullId", "timepoints",
                            "timesteps", "timeunits", "usertimesteps",
                            "flagPower", "_timeframes", "providedID",
                            "nsteps"],
               "results": ["Nt", "At", "AtCurie", "Qt", "flux", "power",
                           "massgr", "totalAtCurie", "totalMassgr", "totalQt",
                           "totalToxIngestion", "totalToxInhalation",
                           "toxicityIngestion", "toxicityInhalation"],
               "xsData": ["fullId", "nIsotopes", "AW", "Q", "BR", "lmbda",
                          "decaymtx", "ingestion", "inhalation", "fymtx",
                          "EfissMeV", "EfissJoule", "xsData",
                          "transmutationmtx"],
               "deplData":["BR", "decaymtx", "ingestion", "inhalation",
                           "lmbda", ],
               "inital": ["volume", "N0", "providedN0"]}




.. _res:


Data Stored on ``Results``
--------------------------


The following complete list of **attributes** are stored on the object:

==================== ==========================================
Attribute							Description
==================== ==========================================
fullId	   						Full ZAID list
-------------------- ------------------------------------------
nIsotopes	   					Number of isotopes
-------------------- ------------------------------------------
AW      		  				Atomic weight
-------------------- ------------------------------------------
Q											Decay heat in W/Bq
-------------------- ------------------------------------------
BR										Branching ratios leading to isomeric
-------------------- ------------------------------------------
lmbda									Decay constant in 1/s
-------------------- ------------------------------------------
decaymtx							Decay matrix
-------------------- ------------------------------------------
ingestion							Ingestion coefficients in Sv/Bq
-------------------- ------------------------------------------
inhalation						Inhalation coefficients in Sv/Bq
-------------------- ------------------------------------------
providedID						User-provided IDs
-------------------- ------------------------------------------
providedN0						User-provided number densities in #/b/cm
-------------------- ------------------------------------------
volume								Volume of the system in cm3
-------------------- ------------------------------------------
usertimesteps					User-defined time-steps in the original units
-------------------- ------------------------------------------
timepoints						Time-points in seconds
-------------------- ------------------------------------------
timesteps							Time-steps in seconds
-------------------- ------------------------------------------
N0										Nuclide densities in #/b/cm at t=t0
-------------------- ------------------------------------------
Nt										Nuclide densities in #/b/cm for all the time-points
-------------------- ------------------------------------------
At										Isotopic activity in Bq for all the time-points
-------------------- ------------------------------------------
AtCuire								Isotopic activity in Curie for all the time-points
-------------------- ------------------------------------------
Qt										Isotopic decay heat in Watts for all the time-points
-------------------- ------------------------------------------
flux									Time-dependent flux in n/cm2/s
-------------------- ------------------------------------------
power									Time-dependent power in Watts
-------------------- ------------------------------------------
massgr								Isotopic mass in grams for all the time-points
-------------------- ------------------------------------------
toxicityIngestion			Isotopic ingestion radiotoxicity in Sv for all the time-points
-------------------- ------------------------------------------
toxicityInhalation		Isotopic inhalation radiotoxicity in Sv for all the time-points
-------------------- ------------------------------------------
totalQt								Time-dependent total decay heat in Watts 
-------------------- ------------------------------------------
totalAtCuries					Time-dependent total activity in Curie 
-------------------- ------------------------------------------
totalMassgr						Time-dependent total mass in grams
-------------------- ------------------------------------------
totalToxIngestion			Time-dependent total ingestion radiotoxicity in Sv
-------------------- ------------------------------------------
totalToxInnhalation		Time-dependent total inhalation radiotoxicity in Sv
==================== ==========================================


.. Note::

	* All the isotopic time-dependent properties (e.g., ``Nt``) are stored as matrices. The row corresponds to a specific isotope and the column to a specific time-point.
	* Time-dependent flux and power are stored as 1-dim vectors.
	* All the total-related attributes (e.g., ``totalQt``) are stored as 1-dim vectors.


The following **methods** are stored on the object:

====================== ==========================================
Method									Description
====================== ==========================================
getvalues				   			Obtain the values for an attribute and specific isotope
---------------------- ------------------------------------------
plot									  Plot the results for a specific attribute and specific isotope
====================== ==========================================

.. _predep:


Data Stored on ``TransmutationData``
------------------------------------ 

The following complete list of **attributes** are stored on the object:

================ ==========================================
Attribute					Description
================ ==========================================
fullId	   				Full ZAID list
---------------- ------------------------------------------
nIsotopes	   			Number of isotopes
---------------- ------------------------------------------
AW      		  		Atomic weight
---------------- ------------------------------------------
Q									Decay heat in W/Bq
---------------- ------------------------------------------
BR								Branching ratios leading to isomeric
---------------- ------------------------------------------
lmbda							Decay constant in 1/s
---------------- ------------------------------------------
decaymtx					Decay matrix
---------------- ------------------------------------------
ingestion					Ingestion coefficients in Sv/Bq
---------------- ------------------------------------------
inhalation				Inhalation coefficients in Sv/Bq
---------------- ------------------------------------------
fymtx				  		Fission yields matrix
---------------- ------------------------------------------
EfissMeV					Energy per fission in MeV
---------------- ------------------------------------------
EfissJoule				Energy per fission in Joules
---------------- ------------------------------------------
xsData						Cross section data matrix
---------------- ------------------------------------------
transmutationmtx	Transmutation matrix without decay
---------------- ------------------------------------------
libraryFlag				True if a pre-generated library is used
================ ==========================================


.. Note::

	* The cross section columns provided in the ``xsData`` follows the structure given in the table below.
	* Not all the attributes must be defined. For example in cases with only decay analysis, the cross section data might not apprear.
	* Some data may be omitted if external libraries are used or the user does not provide a complete set of data (i.e., no ingestion coefficients). The analysis can still be conducted, but in this case the user will not be able to evaluate the radiotoxicity post-irradiation.

Description of the ``xsData`` attribute:

================ ==========================================
Column						Parameter description
================ ==========================================
0	   							ZAID of the isotope (e.g., 541351)
---------------- ------------------------------------------
1	   							Absorption cross section in cm2
---------------- ------------------------------------------
2      		  			Fission cross section in cm2
---------------- ------------------------------------------
3									Capture to ground cross section in cm2
---------------- ------------------------------------------
4									Capture to isometic cross section in cm2
---------------- ------------------------------------------
5									(n, 2n) cross section in cm2
---------------- ------------------------------------------
6									(n, 3n) cross section in cm2
---------------- ------------------------------------------
7									(n, alpha) cross section in cm2
---------------- ------------------------------------------
8									(n, proton) cross section in cm2
---------------- ------------------------------------------
9				  				(n, deuterium) cross section in cm2
---------------- ------------------------------------------
10								(n, tritium) cross section in cm2
================ ==========================================


The following **methods** are stored on the object:

================ ==========================================
Method						Description
================ ==========================================
ReadData	   			Read cross sections, decay, and fission yields data
---------------- ------------------------------------------
Condense	   			Filter the transmutation chains and leave only the data correponding to a specific list of isotopes.
================ ==========================================
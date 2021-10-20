.. _project-overview:

================
Project Overview
================

The custom depletion solver was developed as an engineering tool
to perform quick depletion and decay analysis.
The package consist of two major sub-modules:

* The data generation module, which allows the user to read
  pre-generated data and complement it with provided cross sections.
  This module generates the transmutation matrix leading to constructing
  the Bateman equations.
	
* The solution of the Bateman equations is carried out using built-in
  solvers. Namely, CRAM, EXPM, and ODEINT.

The following features are included as part of the unique capabilities 
of the Python Isotopic Depletion package:

*	Problem definition can be specified only using memory-based
	data transfer.

*	The package contains pre-generated data libraries that include
	decay data, fission yields, and post-irradiation coefficients.

*	The depletion can be carried out using time-invariant cross section sets
	or using multiple sets pre-generated in advance. The cross sections
	can then be linearly interpolated between depletion steps or held
	contant for a given depletion step.

*	The user can easily define tailored transmutation and/or decay chains.

*	A simple post-processing capability is embedded within the package.


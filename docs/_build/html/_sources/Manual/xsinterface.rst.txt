.. _xsinterface:


Cross sections Interface
------------------------ 
An interface to store all the cross-section sets for a range of operational
conditions (up to three dependencies).
The interface also contains interpolation techniques to allow constructing
a new cross section set for given operational conditions traces.

This interface is used to link between cross sections generated in
branch-off calculations for multiple operational points and the actual
operational regime. The interface allows to generate on-the-fly cross
section sets that represent different time-dependent scenarios.


1. The user is required to build cross sections objects in advance.
2. These objects must be linked to operational conditions.		


**Load module**:

.. code::

	from pyIsoDep.functions.xsinterface import XsInterface


**Execution**:
  
.. code::

	xs = XsInterface(numdepn, numpert, states, xssets, extrpFlag)
	
where,

============= ==========================================
Input					Description
============= ==========================================
numdepn	  		Number of dependencies, int
------------- ------------------------------------------
numpert      	Total number of perturbations, int
------------- ------------------------------------------
states				The description of all the perturbations in terms of the dependencies
------------- ------------------------------------------
xssets      	Cross sections objects
------------- ------------------------------------------
extrpFlag     A flag to indicate if extrapolation is allowed, default is True
============= ==========================================

.. Note::

	* ``states`` is a 2-dimensional matrix. Columns correpond to the state [T1, P1, ND1] and the rows correspond to different perturbations.
	* ``xssets`` is a ``TransmutationData`` object.
	* The interpolation is always linear. Therefore, having more pre-generated points is recommended.
	* ``numdepn`` higher than 3 is not supported in the current version.
	* The ``states`` must have a complete set. For example, [[T1, P1], [T1, P2], [T2, P1], [T2, P2]].
  
  
**Example: 1-dim**

.. code::

	xs = XsInterface(numdepn=1, numpert=3, states=[[500], [600], [700]],
                 	 xssets=[xs500, xs600, xs700], extrpFlag=True)

**Example: 2-dim**

.. code::

	xs = XsInterface(numdepn=2, numpert=6,
                   states=[[500, 2], [500, 3], [500, 4],
                           [600, 2], [600, 3], [600, 4]],
                   xssets=[xs1, xs2, xs3, xs4, xs5, xs6], extrpFlag=True)

**Example: 3-dim**

.. code::

	xs = XsInterface(numdepn=3, numpert=18,
                   states=[[500, 2, 1E-05], [500, 2, 2E-05], [500, 2, 3E-05],
                           [500, 3, 1E-05], [500, 3, 2E-05], [500, 3, 3E-05],
                           [500, 5, 1E-05], [500, 5, 2E-05], [500, 5, 3E-05],
                           [600, 2, 1E-05], [600, 2, 2E-05], [600, 2, 3E-05],
                           [600, 3, 1E-05], [600, 3, 2E-05], [600, 3, 3E-05],
                           [600, 5, 1E-05], [600, 5, 2E-05], [600, 5, 3E-05]],
                   xssets=[xs1, xs2, xs3, xs4, xs5, xs6, xs7, xs8, xs9, xs10,
                           xs11, xs12, xs13, xs14, xs15, xs16, xs17, xs18], extrpFlag=True)


============
setTimeTrace
============

The use can feed an operational trace for each of the dependencies
defined in the problem. If only one dependency exists (e.g. pressure),
only a pressure trace as a function needs to be provided.
If more than one dependnecy exists then separate vectors for each
dependnecy are required.

**Execution**:
  
.. code::

	timepoints, xsTimeSets = xs.setTimeTrace(timepoints, *argv)

The following inputs can be inputted to the method:

============= ==========================================
Input					Description
============= ==========================================
timepoints		The time points at which interpolated data will be created.
------------- ------------------------------------------
argv	   			Multiple time-dependent traces for each dependency
============= ==========================================


.. Note::

	* The number of trace vectors must correspond to the number of dependencies in the problem.
	* All the provided vectors must be of equal size.
	* The first entrence will always represent time, and the following parameters represent the time-dependent traces.

The method stores attributes/data on the object, but also returns:

============= ==========================================
Output				Description
============= ==========================================
timepoints		Time points
------------- ------------------------------------------
xsTimeSets	  Multiple time-dependent ``TransmutationData`` objects needed for depletion
============= ==========================================


**Example: 1-dependency**

.. code::

	timepoints, xsTimeSets = xs.setTimeTrace([0, 3.5, 7.0], [525, 550, 575])
	

**Example: 2-dependencies**

.. code::

	timepoints, xsTimeSets = xs.setTimeTrace([0, 3.5], [500, 550], [3.0, 3.5])
	
**Example: 3-dependencies**

.. code::

	timepoints, xsTimeSets =\
    xs.setTimeTrace([0, 3.5], [525, 550], [2.5, 3.5], [1.5E-05, 2.5E-05])


**Very Important Note**
 
The ``timepoints``, and ``xsTimeSets`` can  be used to define the depletion data in the following manner:

.. code::

	dep = MainDepletion(timepoints, *xsTimeSets)
	
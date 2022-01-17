.. _user-manual:


User's Manual
=============

The User\'s Manual for the ``pyIsoDep`` package describes how 
to use the code and its various transmutation and decay applications. 

Execution Sequence
------------------

The following sequence describes all the stages required to 
perform a transmutation and/or decay simulation.

	1. :ref:`Decay and transmutation data generation <datagen>`: define a single or multiple ``TransmutationData`` containers.
		
		1.1. The user can move directly to step 2, or
		
		1.2. Define pre-generated cross sections sets and time-dependent traces using :ref:`xsinterface`.
		
	2. :ref:`Execute depletion <deplete>`: set depletion or decay history, and solve the Bateman equations.
	3. :ref:`Post-process results <postprocess>`: a dedicated container to store only result attributes, and methods to obtain specific values and plot results.

		
Description of Stored Data 
^^^^^^^^^^^^^^^^^^^^^^^^^^
Following a successive execution of each of the above stages, the containers will be populated with data unique to each container.
The following sub-sections describe the data stored on each container:

	1. :ref:`Pre- decay and transmutation data storage <predep>`
	2. :ref:`Post- decay and transmutation data storage <postdep>`
	3. :ref:`Results container <res>`

Abstract Example 
^^^^^^^^^^^^^^^^


**Data generation**

.. code::

	data = TransmutationData(...)
	data.ReadData(...)

**Pre-generated XS data and traces** (not mandatory)

.. code::

	xs = XsInterface(..., data)
	timepoints, xsTimeSets = xs.setTimeTrace(...)

	
**Calling Depletion**

.. code::

	dep = MainDepletion(..., data or xsTimeSets)
	dep.SetDepScenario(...)
	dep.SetInitialComposition(...)
	dep.SolveDepletion(method="cram")
	
**Post-process results**

.. code::

	res = Results(dep)
	res.getvalues(...)
	res.plot(...)

	


.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Contents:

   Datageneration.rst
   Depletion.rst
   Postprocess.rst
   predepletion.rst
   postdepletion.rst
   resstorage.rst
   xsinterface.rst


	

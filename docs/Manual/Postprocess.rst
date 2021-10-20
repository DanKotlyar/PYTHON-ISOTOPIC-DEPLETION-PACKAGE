.. _postprocess:


Post-process
------------ 

Basic post-processing capabilities to get specific values and plot results.
Time-dependent data can be obtained for individual isotopes. Two methods are available:

=================== ==========================================
Input					Description
=================== ==========================================
:ref:`getvalues`	  Obtain the values of a specific property
------------------- ------------------------------------------
:ref:`plot`	      	Plot time-dependent results
=================== ==========================================

Before the post-processing class is applied, the user must execute a :ref:`depletion scenario <deplete>`    
and generate a :ref:`MainDepletion <postdep>` object.

**Load module**:

.. code::

	pyIsoDep.functions.postprocessresults import Results

**Execution**:
  
.. code::

	res = Results(results)

where,

``results`` is the :ref:`MainDepletion <postdep>` object.


.. Note::

	* The ``Results`` will be very similar to the :ref:`MainDepletion <postdep>` object, but will contain only the attributes (not the methods inherited from the ``MainDepletion`` object).  



.. _getvalues:


getvalues
^^^^^^^^^

Obtain the values of a specific property

.. code::

	res.getvalues(attribute, isotopes)
	
where,

============= ==========================================
Input					Description
============= ==========================================
attribute			The name of the property (provided as a string)
------------- ------------------------------------------
isotopes			A list or array of isotopes provided in ZZAAA0/1 format (e.g., 541351)
============= ==========================================

.. Note::

	* For the full list of attributes that can exist on the ``Results`` container, please see :ref:`Results container <res>`.
	* If the :ref:`post-irradiation supplementary methods <res>` are not executed then some attributes will not exist on the ``Results`` object. 



**Examples:**

* Isotopic density for a specific isotope:

.. code::

	dep = MainDepletion(0.0, data)
	res = Results(dep)
	res.getvalues("Nt", isotopes=[541350])
	
* Decay heat in Watts for multiple isotopes:

.. code::

	res.getvalues("Qt", isotopes=[541350, 942380])


* Material total values:

.. code::

	res.getvalues("totalQt",)
	


.. _plot:


plot
^^^^

Plot time-dependent results.

.. code::

	res.plot(attribute, timeUnits, isotopes, xlabel, ylabel, norm, fontsize, markers, markerfill, markersize=6, pltType, newFig)
	
where,

============= ==========================================
Input					Description
============= ==========================================
attribute			The name of the property (provided as a string)
------------- ------------------------------------------
timeUnits			String for the units to be used on the x-axis = {"seconds", "minutes", "hours", "days"}
------------- ------------------------------------------
isotopes			A list or array of isotopes provided in ZZAAA0/1 format (e.g., 541351)
------------- ------------------------------------------
xlabel				String for the x-axis label
------------- ------------------------------------------
ylabel				String for the y-axis label
------------- ------------------------------------------
norm					Normalization factor. Should be provided as float. It is used to divide the desired results by the provided factor.
------------- ------------------------------------------
markers				String of List of strings to represent markers. Default is "--^"
------------- ------------------------------------------
markerfill		Boolean flag to indicate whetehr marker fill should be removed. Default is ``False``
------------- ------------------------------------------
markersize		Float that represents the marker size.
------------- ------------------------------------------
pltType				String to represent a plot type: {"linear", "loglog", "semilogx"}. Default is "linear"
------------- ------------------------------------------
newFig				Boolean flag to indicate whether a new figure should be created or an existing can be used. Default is ``True``
============= ==========================================

.. Note::

	* For the full list of attributes that can exist on the ``Results`` container, please see :ref:`Results container <res>`.
	* If the :ref:`post-irradiation supplementary methods <res>` are not executed then some attributes will not exist on the ``Results`` object. 



**Examples:**

* Isotopic density for a specific isotope:

.. code::

	dep = MainDepletion(0.0, data)
	res = Results(dep)
	res.plot("Nt", isotopes=[541350])
	
* Decay heat in Watts for multiple isotopes:

.. code::

	res.plot("Qt", isotopes=[531350, 541350], norm=1E+6, ylabel="Total Decay Heat, MW")


* Non-isotopic total values:

.. code::

	res.plot("flux", ylabel="Flux, n/cm2/s", pltType="semilogx")
	



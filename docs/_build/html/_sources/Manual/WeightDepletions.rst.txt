.. _weightdep:


Weighting Different Solutions
-----------------------------
 
The user may want to solve multiple problems using :ref:`deplete`.
Each problem can be set and solved separately.
All these individual problems can be related (e.g., part of a complete problem).
For example, different components in the reactor-core.
The user can combine these depletion objects by applying the weighting method.
The weighting method corrects according to volume.

**Example**

.. code::
	
	from pyIsoDep.functions.weightdepletionobjects import WeightDepObjects
	weightedDep = WeightDepObjects(dep1, dep2, dep3, ...)

where, ``dep1``, ... are the depletion objects.


.. Note::

	* All the ``dep`` objects must have exactly the same structure and order.
	* Specifically, the utilized transmutation chains must be identical. However, the data provided for each container may be different.
	* The time points used for each depletion object must be identocal. However, the traces values can be different.
	* It is recommended that the same post-processing analysis be applied to all the depletion object. 
	* The wighted depletion object ``weightedDep`` can be post-processed using :ref:`Post-process results <postprocess>` capability.

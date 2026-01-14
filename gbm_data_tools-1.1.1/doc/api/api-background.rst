.. _api-background:

Background Estimation Module
============================
The classes in `gbm.background` are standard background estimation classes
available to fit time history background for pre-binned and unbinned data.

Fitter Interface
----------------
These fitting classes are plugins for :class:`~gbm.background.BackgroundFitter`,
which outputs a :class:`~gbm.background.BackgroundRates` object that can be 
integrated in time to produce a :class:`~gbm.background.BackgroundSpectrum`, 
which is a background spectrum.

BackgroundFitter
^^^^^^^^^^^^^^^^
.. autoclass:: gbm.background.BackgroundFitter
   :show-inheritance:
   :members:
   :inherited-members:
----

BackgroundRates
^^^^^^^^^^^^^^^
.. autoclass:: gbm.background.BackgroundRates
   :show-inheritance:
   :members:
   :inherited-members:
----

BackgroundSpectrum
^^^^^^^^^^^^^^^^^^
.. autoclass:: gbm.background.BackgroundSpectrum
   :show-inheritance:
   :members:
   :inherited-members:
----

For Pre-Binned Data
-------------------
To be used for fitting the background of pre-binned data.

Polynomial
^^^^^^^^^^
.. autoclass:: gbm.background.binned.Polynomial
   :show-inheritance:
   :members:
   :inherited-members:
----

For Unbinned Data
-----------------
To be used for fitting the background of unbinned data.

NaivePoisson
^^^^^^^^^^^^
.. autoclass:: gbm.background.unbinned.NaivePoisson
   :show-inheritance:
   :members:
   :inherited-members:
----


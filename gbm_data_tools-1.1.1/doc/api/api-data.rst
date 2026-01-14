.. _api-data:

Data Module
===========
The classes in ``gbm.data`` are interfaces to the data that GBM produces.

Primary GBM Data Classes
------------------------
The primary data classes for GBM data:

* Lightcurves/Spectra (:class:`~gbm.data.Cspec`, :class:`~gbm.data.Ctime`, 
  :class:`~gbm.data.TTE`)
* Detector Responses (:class:`~gbm.data.RSP`)
* Fermi Position History (:class:`~gbm.data.PosHist`)
* Skymaps (:class:`~gbm.data.GbmHealPix`)
* Quicklook Trigger Data (:class:`~gbm.data.Trigdat`, :class:`~gbm.data.Tcat`)

Cspec 
^^^^^

.. autoclass:: gbm.data.Cspec
   :show-inheritance:
   :members:
   :inherited-members:

----

Ctime 
^^^^^

.. autoclass:: gbm.data.Ctime
   :show-inheritance:
   :members:
   :inherited-members:

----

GbmHealPix 
^^^^^^^^^^

.. autoclass:: gbm.data.GbmHealPix
   :show-inheritance:
   :members:
   :inherited-members:

----

PosHist
^^^^^^^

.. autoclass:: gbm.data.PosHist
   :show-inheritance:
   :members:

----


RSP
^^^

.. autoclass:: gbm.data.RSP
   :show-inheritance:
   :members:
   :inherited-members:

----

Tcat
^^^^

.. autoclass:: gbm.data.Tcat
   :show-inheritance:
   :members:

----

Trigdat
^^^^^^^
.. autoclass:: gbm.data.Trigdat
   :show-inheritance:
   :members:

----

TTE
^^^
.. autoclass:: gbm.data.TTE
   :show-inheritance:
   :members:
   :inherited-members:

----


Scat
^^^^
.. autoclass:: gbm.data.Scat
   :show-inheritance:
   :members:
   :inherited-members:

----

Base Data Classes
-----------------
Some base data classes.  These are not intended to be instantiated, but can
be useful to inherit if designing new data classes.

DataFile
^^^^^^^^
.. autoclass:: gbm.data.data.DataFile
   :members:
   :inherited-members:

----

HealPix
^^^^^^^
.. autoclass:: gbm.data.HealPix
   :show-inheritance:
   :members:
   :inherited-members:

----

PHAII
^^^^^
.. autoclass:: gbm.data.phaii.PHAII
   :show-inheritance:
   :members:
   :inherited-members:

----

Data Class Miscellania
----------------------
Miscellaneous data classes: 
Collections of data, single-spectrum PHA and BAK classes, and Trigdat info classes 

DataCollection
^^^^^^^^^^^^^^
.. autoclass:: gbm.data.DataCollection
   :show-inheritance:
   :members:
   :inherited-members:
----

GbmDetectorCollection
^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: gbm.data.GbmDetectorCollection
   :show-inheritance:
   :members:
   :inherited-members:
----

BAK
^^^
.. autoclass:: gbm.data.BAK
   :show-inheritance:
   :members:
   :inherited-members:
----

PHA
^^^
.. autoclass:: gbm.data.PHA
   :show-inheritance:
   :members:
   :inherited-members:
----

scat.Parameter
^^^^^^^^^^^^^^
.. autoclass:: gbm.data.scat.Parameter
   :show-inheritance:
   :members:
   :inherited-members:
----

scat.GbmModelFit
^^^^^^^^^^^^^^^^
.. autoclass:: gbm.data.scat.GbmModelFit
   :show-inheritance:
   :members:
   :inherited-members:
----

scat.DetectorData
^^^^^^^^^^^^^^^^^
.. autoclass:: gbm.data.scat.DetectorData
   :show-inheritance:
   :members:
   :inherited-members:
----


trigdat.MaxRates
^^^^^^^^^^^^^^^^
.. autoclass:: gbm.data.trigdat.MaxRates
   :show-inheritance:
   :members:
   :inherited-members:
----

trigdat.BackRates
^^^^^^^^^^^^^^^^^
.. autoclass:: gbm.data.trigdat.BackRates
   :show-inheritance:
   :members:
   :inherited-members:
----

trigdat.FswLocation 
^^^^^^^^^^^^^^^^^^^
.. autoclass:: gbm.data.trigdat.FswLocation
   :show-inheritance:
   :members:
   :inherited-members:
----

Data Primitives
---------------
Primitive data classes.  The GBM Primary classes act as wrappers around some of
these classes with additional specialized metadata.  These classes may be useful
in designing new higher-level data classes.

* 1D Binned data (:class:`~gbm.data.primitives.Bins`, 
  :class:`~gbm.data.primitives.EnergyBins`, :class:`~gbm.data.primitives.TimeBins`)
* 2D Binned data (:class:`~gbm.data.primitives.TimeEnergyBins`)
* Unbinned event data (:class:`~gbm.data.primitives.EventList`)
* Time range data (:class:`~gbm.data.primitives.TimeRange`, 
  :class:`~gbm.data.primitives.GTI`)

Bins
^^^^
.. autoclass:: gbm.data.primitives.Bins
   :show-inheritance:
   :members:
   :inherited-members:
----

EnergyBins
^^^^^^^^^^
.. autoclass:: gbm.data.primitives.EnergyBins
   :show-inheritance:
   :members:
   :inherited-members:
----

TimeBins
^^^^^^^^
.. autoclass:: gbm.data.primitives.TimeBins
   :show-inheritance:
   :members:
   :inherited-members:
----

TimeEnergyBins
^^^^^^^^^^^^^^
.. autoclass:: gbm.data.primitives.TimeEnergyBins
   :show-inheritance:
   :members:
   :inherited-members:
----

EventList
^^^^^^^^^
.. autoclass:: gbm.data.primitives.EventList
   :show-inheritance:
   :members:
   :inherited-members:
----

GTI
^^^
.. autoclass:: gbm.data.primitives.GTI
   :show-inheritance:
   :members:
   :inherited-members:
----

TimeRange
^^^^^^^^^
.. autoclass:: gbm.data.primitives.TimeRange
   :show-inheritance:
   :members:
   :inherited-members:
----

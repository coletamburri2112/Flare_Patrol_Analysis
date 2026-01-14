.. _api-utils:

Utilities
=================
Several useful utilities used throughout the various Data Tools modules


Coordinate/Geometry
-------------------
The following functions in ``gbm.coords`` are used in various coordinate 
transforms and geometry calculations.

.. automodule:: gbm.coords
    :members:
----

Detector Definitions
--------------------
The class in ``gbm.detector`` contains the GBM Detector parameters.

.. autoclass:: gbm.detectors.Detector
    :show-inheritance:
    :members:
    :undoc-members: 
    :exclude-members: B0, B1, N0, N1, N2, N3, N4, N5, N6, N7, N8, N9, NA, NB,
                      pointing, short_name
----

Filenaming
----------
The class in ``gbm.file`` contains the standardized definitions for the GBM
file naming scheme and other associated utilities

.. automodule:: gbm.file
    :members:

----

Time Conversions
----------------
The class in ``gbm.time`` contains Fermi MET time definition, a class for time
conversions to and from Fermi MET and other time utilities.

.. automodule:: gbm.time
    :show-inheritance:
    :members:
    :inherited-members:

----


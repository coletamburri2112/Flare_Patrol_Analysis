.. _api-finders:

Data Finders
============
The classes in ``gbm.finder`` are for finding data in the GBM trigger and 
continuous data FTP databases, and for searching the online GBM catalogs.

Base Classes
------------
These base classes define an interface to the HEASARC FTP data server for GBM
data and the online Browse catalog.

FtpFinder
^^^^^^^^^
.. autoclass:: gbm.finder.FtpFinder
   :show-inheritance:
   :members:
   :inherited-members:

----

HeasarcBrowse
^^^^^^^^^^^^^
.. autoclass:: gbm.finder.HeasarcBrowse
   :show-inheritance:
   :members:
   :inherited-members:

----


FTP Finders
-----------
Classes to find and download GBM data. 
Inherits from :class:`~gbm.finder.FtpFinder`. 

TriggerFtp
^^^^^^^^^^
.. autoclass:: gbm.finder.TriggerFtp
   :show-inheritance:
   :members:
   :inherited-members:

----

ContinuousFtp
^^^^^^^^^^^^^
.. autoclass:: gbm.finder.ContinuousFtp
   :show-inheritance:
   :members:
   :inherited-members:

----


HEASARC Catalogs
----------------
Classes to download and query the GBM trigger and burst catalogs.


TriggerCatalog
^^^^^^^^^^^^^^
.. autoclass:: gbm.finder.TriggerCatalog
   :show-inheritance:
   :members:
   :inherited-members:

----

BurstCatalog
^^^^^^^^^^^^
.. autoclass:: gbm.finder.BurstCatalog
   :show-inheritance:
   :members:
   :inherited-members:

----
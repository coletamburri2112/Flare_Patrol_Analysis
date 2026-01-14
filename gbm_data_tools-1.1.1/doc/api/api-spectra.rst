.. _api-spectra:

Spectra Module
==============
The classes and methods in ``gbm.spectra`` are for spectral fitting.


Spectral Fitting
----------------
These classes and functions define the likelihood functions and 
maximum-likelihood fitters associated with those functions.


SpectralFitter 
^^^^^^^^^^^^^^
.. autoclass:: gbm.spectra.fitting.SpectralFitter
   :show-inheritance:
   :members:
   :inherited-members:

----

SpectralFitterChisq
^^^^^^^^^^^^^^^^^^^
.. autoclass:: gbm.spectra.fitting.SpectralFitterChisq
   :show-inheritance:
   :members:

----

chisq
^^^^^
.. autofunction:: gbm.spectra.fitting.chisq

SpectralFitterCstat
^^^^^^^^^^^^^^^^^^^
.. autoclass:: gbm.spectra.fitting.SpectralFitterCstat
   :show-inheritance:
   :members:

----

cstat
^^^^^
.. autofunction:: gbm.spectra.fitting.cstat

SpectralFitterPgstat
^^^^^^^^^^^^^^^^^^^^
.. autoclass:: gbm.spectra.fitting.SpectralFitterPgstat
   :show-inheritance:
   :members:

----

pgstat
^^^^^^
.. autofunction:: gbm.spectra.fitting.pgstat

----

SpectralFitterPstat
^^^^^^^^^^^^^^^^^^^
.. autoclass:: gbm.spectra.fitting.SpectralFitterPstat
   :show-inheritance:
   :members:

----

pstat
^^^^^
.. autofunction:: gbm.spectra.fitting.pstat


----

Spectral Function Base Classes
------------------------------
The base classes represent photon model functions with metadata that allow these
functions to be fit to data.  Specifically, :class:`~gbm.spectra.functions.Function`
is the base class for all functions and :class:`~gbm.spectra.functions.SuperFunction`
inherits from :class:`~gbm.spectra.functions.Function` and is the base class
for either a sum or product of functions.

Function
^^^^^^^^
.. autoclass:: gbm.spectra.functions.Function
   :show-inheritance:
   :members:
   :inherited-members:

----

SuperFunction
^^^^^^^^
.. autoclass:: gbm.spectra.functions.SuperFunction
   :show-inheritance:
   :members:
   :inherited-members:

----

Spectral Functions
------------------
The functions listed below inherit from :class:`~gbm.spectra.functions.Function`
and have the optional metadata described in that class documentation.  

PowerLaw
^^^^^^^^
.. autoclass:: gbm.spectra.functions.PowerLaw
   :show-inheritance:
   
----

Comptonized
^^^^^^^^^^^
.. autoclass:: gbm.spectra.functions.Comptonized
   :show-inheritance:
   
----

Band
^^^^
.. autoclass:: gbm.spectra.functions.Band
   :show-inheritance:
   
----

BandOld
^^^^^^^
.. autoclass:: gbm.spectra.functions.BandOld
   :show-inheritance:
   
----

BrokenPowerLaw
^^^^^^^^^^^^^^
.. autoclass:: gbm.spectra.functions.BrokenPowerLaw
   :show-inheritance:
   
----

DoubleBrokenPowerLaw
^^^^^^^^^^^^^^^^^^^^
.. autoclass:: gbm.spectra.functions.DoubleBrokenPowerLaw
   :show-inheritance:
   
----

SmoothlyBrokenPowerLaw
^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: gbm.spectra.functions.SmoothlyBrokenPowerLaw
   :show-inheritance:
   
----


LogNormal
^^^^^^^^^
.. autoclass:: gbm.spectra.functions.LogNormal
   :show-inheritance:
   
----

GaussianLog
^^^^^^^^^^^
.. autoclass:: gbm.spectra.functions.GaussianLog
   :show-inheritance:
   
----

GaussianLogVaryingFWHM
^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: gbm.spectra.functions.GaussianLogVaryingFWHM
   :show-inheritance:
   
----

SunyaevTitarchuk
^^^^^^^^^^^^^^^^
.. autoclass:: gbm.spectra.functions.SunyaevTitarchuk
   :show-inheritance:
   
----

OTTB
^^^^
.. autoclass:: gbm.spectra.functions.OTTB
   :show-inheritance:
   
----

OTTS
^^^^
.. autoclass:: gbm.spectra.functions.OTTS
   :show-inheritance:
   
----

BlackBody
^^^^^^^^^
.. autoclass:: gbm.spectra.functions.BlackBody
   :show-inheritance:
   
----

GaussLine
^^^^^^^^^
.. autoclass:: gbm.spectra.functions.GaussLine
   :show-inheritance:
   
----

YangSoongPulsar
^^^^^^^^^^^^^^^
.. autoclass:: gbm.spectra.functions.YangSoongPulsar
   :show-inheritance:
   
----

TanakaPulsar
^^^^^^^^^^^^
.. autoclass:: gbm.spectra.functions.TanakaPulsar
   :show-inheritance:
   
----

LowEnergyCutoff
^^^^^^^^^^^^^^^
.. autoclass:: gbm.spectra.functions.LowEnergyCutoff
   :show-inheritance:
   
----

HighEnergyCutoff
^^^^^^^^^^^^^^^^
.. autoclass:: gbm.spectra.functions.HighEnergyCutoff
   :show-inheritance:
   
----

PowerLawMult
^^^^^^^^^^^^
.. autoclass:: gbm.spectra.functions.PowerLawMult
   :show-inheritance:
   
----

GaussLineMult
^^^^^^^^^^^^^
.. autoclass:: gbm.spectra.functions.GaussLineMult
   :show-inheritance:
   
----

LorentzLineMult
^^^^^^^^^^^^^^^
.. autoclass:: gbm.spectra.functions.LorentzLineMult
   :show-inheritance:
   
----

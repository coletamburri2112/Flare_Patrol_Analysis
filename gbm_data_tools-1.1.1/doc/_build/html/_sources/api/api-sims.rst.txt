.. _api-sims:

Simulate Module
===============
The classes and methods in ``gbm.simulate`` are for simulating spectra and
event data.


PHA Simulation
--------------
Simulate count spectra by folding a photon model through a detector response.
Also simulates background from a background model.

PhaSimulator
^^^^^^^^^^^^
.. autoclass:: gbm.simulate.PhaSimulator
   :show-inheritance:
   :members:
   :inherited-members:

----

TTE Simulation
--------------
Simulate event data generated from a spectral and temporal source model, 
as well as event data from a spectral and temporal background model.

TteSourceSimulator
^^^^^^^^^^^^^^^^^^
.. autoclass:: gbm.simulate.TteSourceSimulator
   :show-inheritance:
   :members:
   :inherited-members:

----

TteBackgroundSimulator
^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: gbm.simulate.TteBackgroundSimulator
   :show-inheritance:
   :members:
   :inherited-members:

----


Simulation Generators
---------------------
A variety of specialized random generators

PoissonBackgroundGenerator
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: gbm.simulate.generators.PoissonBackgroundGenerator
   :show-inheritance:
   :members:
   :inherited-members:

----

GaussianBackgroundGenerator
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: gbm.simulate.generators.GaussianBackgroundGenerator
   :show-inheritance:
   :members:
   :inherited-members:

----

SourceSpectrumGenerator
^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: gbm.simulate.generators.SourceSpectrumGenerator
   :show-inheritance:
   :members:
   :inherited-members:

----

VariablePoissonBackground
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: gbm.simulate.generators.VariablePoissonBackground
   :show-inheritance:
   :members:
   :inherited-members:

----

VariableGaussianBackground
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: gbm.simulate.generators.VariableGaussianBackground
   :show-inheritance:
   :members:
   :inherited-members:

----

VariableSourceSpectrumGenerator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: gbm.simulate.generators.VariableSourceSpectrumGenerator
   :show-inheritance:
   :members:
   :inherited-members:

----

EventSpectrumGenerator
^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: gbm.simulate.generators.EventSpectrumGenerator
   :show-inheritance:
   :members:
   :inherited-members:

----

Time-Profile Functions
----------------------
These functions are provided to define the time profile of background and 
sources.

.. automodule:: gbm.simulate.profiles
    :show-inheritance:
    :members:
    :inherited-members:

----

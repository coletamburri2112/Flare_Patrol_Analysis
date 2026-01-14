.. _api-plot:

Plotting Module
===============
The classes and methods in ``gbm.plot`` are for plotting and visualizing the
GBM science and auxiliary data.

Base Classes
-------------
These base classes define specialized plot element and a collection of
various plot elements (forming a plot), respectively.  

GbmPlotElement
^^^^^^^^^^^^^^
.. autoclass:: gbm.plot.gbmplot.GbmPlotElement
   :show-inheritance:
   :members:
   :inherited-members:

----

GbmPlot
^^^^^^^
.. autoclass:: gbm.plot.gbmplot.GbmPlot
   :show-inheritance:
   :members:
   :inherited-members:

----


GBM Plot Classes
----------------
These classes inherit from :class:`~gbm.plot.gbmplot.GbmPlot` and represent
standard visualizations of GBM science data and auxiliary data.


Lightcurve
^^^^^^^^^^
.. autoclass:: gbm.plot.Lightcurve
   :show-inheritance:
   :members:
   :inherited-members:

----

Spectrum
^^^^^^^^
.. autoclass:: gbm.plot.Spectrum
   :show-inheritance:
   :members:
   :inherited-members:

----

ModelFit
^^^^^^^^
.. autoclass:: gbm.plot.ModelFit
   :show-inheritance:
   :members:
   :inherited-members:

----

ResponseMatrix
^^^^^^^^^^^^^^
.. autoclass:: gbm.plot.ResponseMatrix
   :show-inheritance:
   :members:
   :inherited-members:

----

PhotonEffectiveArea
^^^^^^^^^^^^^^^^^^^
.. autoclass:: gbm.plot.PhotonEffectiveArea
   :show-inheritance:
   :members:
   :inherited-members:

ChannelEffectiveArea
^^^^^^^^^^^^^^^^^^^^
.. autoclass:: gbm.plot.ChannelEffectiveArea
   :show-inheritance:
   :members:
   :inherited-members:

----

SkyPlot
^^^^^^^
.. autoclass:: gbm.plot.SkyPlot
   :show-inheritance:
   :members:
   :inherited-members:

----

FermiSkyPlot
^^^^^^^^^^^^
.. autoclass:: gbm.plot.FermiSkyPlot
   :show-inheritance:
   :members:
   :inherited-members:

----

EarthPlot
^^^^^^^^^
.. autoclass:: gbm.plot.EarthPlot
   :show-inheritance:
   :members:
   :inherited-members:

----

GBM Plot Element Classes
------------------------
These classes inherit from :class:`~gbm.plot.gbmplot.GbmPlotElemant` and 
represent individual standard plot elements.


Histo
^^^^^
.. autoclass:: gbm.plot.gbmplot.Histo
   :show-inheritance:
   :members:
   :inherited-members:

----

HistoErrorbars
^^^^^^^^^^^^^^
.. autoclass:: gbm.plot.gbmplot.HistoErrorbars
   :show-inheritance:
   :members:
   :inherited-members:

----

HistoFilled
^^^^^^^^^^^
.. autoclass:: gbm.plot.gbmplot.HistoFilled
   :show-inheritance:
   :members:
   :inherited-members:

----

LightcurveBackground
^^^^^^^^^^^^^^^^^^^^
.. autoclass:: gbm.plot.gbmplot.LightcurveBackground
   :show-inheritance:
   :members:
   :inherited-members:

----

SpectrumBackground
^^^^^^^^^^^^^^^^^^
.. autoclass:: gbm.plot.gbmplot.SpectrumBackground
   :show-inheritance:
   :members:
   :inherited-members:

----

HeatMap
^^^^^^^
.. autoclass:: gbm.plot.gbmplot.HeatMap
   :show-inheritance:
   :members:
   :inherited-members:

----

EffectiveArea
^^^^^^^^^^^^^
.. autoclass:: gbm.plot.gbmplot.EffectiveArea
   :show-inheritance:
   :members:
   :inherited-members:

----

SAA
^^^
.. autoclass:: gbm.plot.gbmplot.SAA
   :show-inheritance:
   :members:
   :inherited-members:

----

McIlwainL
^^^^^^^^^
.. autoclass:: gbm.plot.gbmplot.McIlwainL
   :show-inheritance:
   :members:
   :inherited-members:

----

EarthLine
^^^^^^^^^
.. autoclass:: gbm.plot.gbmplot.EarthLine
   :show-inheritance:
   :members:
   :inherited-members:

----

EarthPoints
^^^^^^^^^^^
.. autoclass:: gbm.plot.gbmplot.EarthPoints
   :show-inheritance:
   :members:
   :inherited-members:

----

FermiIcon
^^^^^^^^^
.. autoclass:: gbm.plot.gbmplot.FermiIcon
   :show-inheritance:
   :members:
   :inherited-members:

----

SkyPoints
^^^^^^^^^
.. autoclass:: gbm.plot.gbmplot.SkyPoints
   :show-inheritance:
   :members:
   :inherited-members:

----

SkyLine
^^^^^^^
.. autoclass:: gbm.plot.gbmplot.SkyLine
   :show-inheritance:
   :members:
   :inherited-members:

----

SkyCircle
^^^^^^^^^
.. autoclass:: gbm.plot.gbmplot.SkyCircle
   :show-inheritance:
   :members:
   :inherited-members:

----

SkyPolygon
^^^^^^^^^^
.. autoclass:: gbm.plot.gbmplot.SkyPolygon
   :show-inheritance:
   :members:
   :inherited-members:

----

SkyAnnulus
^^^^^^^^^^
.. autoclass:: gbm.plot.gbmplot.SkyAnnulus
   :show-inheritance:
   :members:
   :inherited-members:

----

Sun
^^^
.. autoclass:: gbm.plot.gbmplot.Sun
   :show-inheritance:
   :members:
   :inherited-members:

----

DetectorPointing
^^^^^^^^^^^^^^^^
.. autoclass:: gbm.plot.gbmplot.DetectorPointing
   :show-inheritance:
   :members:
   :inherited-members:

----

GalacticPlane
^^^^^^^^^^^^^
.. autoclass:: gbm.plot.gbmplot.GalacticPlane
   :show-inheritance:
   :members:
   :inherited-members:

----

SkyHeatmap
^^^^^^^^^^
.. autoclass:: gbm.plot.gbmplot.SkyHeatmap
   :show-inheritance:
   :members:
   :inherited-members:

----

ModelData
^^^^^^^^^
.. autoclass:: gbm.plot.gbmplot.ModelData
   :show-inheritance:
   :members:
   :inherited-members:

----

ModelSamples
^^^^^^^^^^^^
.. autoclass:: gbm.plot.gbmplot.ModelSamples
   :show-inheritance:
   :members:
   :inherited-members:

----

Collection
^^^^^^^^^^
.. autoclass:: gbm.plot.gbmplot.Collection
   :show-inheritance:
   :members:
   :inherited-members:

----

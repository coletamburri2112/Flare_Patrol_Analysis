.. _api-binning:

Data Binning Module
===================
The functions in ``gbm.binning`` are standard functions available to rebin
pre-binned data and to bin unbinned data.


For Pre-Binned Data
-------------------
The following functions in ``gbm.binning.binned`` can be used to rebin 
pre-binned data.

User-defined algorithms can be used as long as the function has the following
inputs:

1. ``old_counts``:  The current array of counts in each bin  
2. ``old_exposure``: The current array of exposure for each bin  
3. ``old_edges``: The current array of bin edges  

Additionally, you can define any algorithm-specific positional arguments after
the required inputs.

The required outputs are:

1. ``new_counts``:  The new array of counts in the rebinned data  
2. ``new_exposure``: The new array of exposure for the rebinned data  
3. ``new_edges``: The new array of bin edges for the binned data
 
Essentially, the function should read in the old histogram data and provide 
updated values for the new histogram. As long as you design your function like 
this, you can use it to bin PHAII data.

----

combine_by_factor
^^^^^^^^^^^^^^^^^
.. autofunction:: gbm.binning.binned.combine_by_factor
----

combine_into_one
^^^^^^^^^^^^^^^^
.. autofunction:: gbm.binning.binned.combine_into_one
----

rebin_by_edge
^^^^^^^^^^^^^
.. autofunction:: gbm.binning.binned.rebin_by_edge
----

rebin_by_edge_index
^^^^^^^^^^^^^^^^^^^
.. autofunction:: gbm.binning.binned.rebin_by_edge_index
----

rebin_by_snr
^^^^^^^^^^^^
.. autofunction:: gbm.binning.binned.rebin_by_snr
----

rebin_by_time
^^^^^^^^^^^^^
.. autofunction:: gbm.binning.binned.rebin_by_time
----

For Unbinned Data
-------------------
The following functions in ``gbm.binning.unbinned`` can be used to bin 
unbinned data

User-defined algorithms can be used as long as the function has the following
required input:

1. ``times``: The array of event times  

Additionally, you can define any algorithm-specific positional or keyword
arguments after the required inputs.

The required output is:

1. ``bin_edges``: The bin edges
 
The algorithm should read in the event times and return the bin edges.  
A function designed like this can be used to bin event data.

----

bin_by_edges
^^^^^^^^^^^^
.. autofunction:: gbm.binning.unbinned.bin_by_edges
----

bin_by_snr
^^^^^^^^^^
.. autofunction:: gbm.binning.unbinned.bin_by_snr
----

bin_by_time
^^^^^^^^^^^
.. autofunction:: gbm.binning.unbinned.bin_by_time
----

combine_by_factor
^^^^^^^^^^^^^^^^^
.. autofunction:: gbm.binning.unbinned.combine_by_factor
----

combine_into_one
^^^^^^^^^^^^^^^^
.. autofunction:: gbm.binning.unbinned.combine_into_one
----

time_to_spill
^^^^^^^^^^^^^
.. autofunction:: gbm.binning.unbinned.time_to_spill
----

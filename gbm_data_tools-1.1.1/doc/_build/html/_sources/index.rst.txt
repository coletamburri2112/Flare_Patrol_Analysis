Welcome to the Fermi GBM Data Tools documentation!
==================================================

.. figure:: images/gdt_logo.png
   
   Hello, I'm Fermi.  Pleased to meet you!

The Fermi GBM Data Tools is an Application Programming Interface (API) for 
GBM data.  The fundamental purpose of the Data Tools is to allow general users
to incorporate GBM analysis into their scripts and workflows without having to
sweat very many details.  To this end, the Data Tools have a fairly high-level
API layer allowing a user to read, reduce, and visualize GBM data with only a 
few lines of code.  For expert users, and users who want fine control over 
various aspects of their analysis, the Data Tools exposes a lower-level API 
layer, which can also be used to generalize the GBM Data Tools to data from 
other like instruments.


.. rubric:: Architecture
The Data Tools are designed with generalization in mind.  Underlying the 
science data interfaces are primitive data types that allow similar data to be
used with the Data Tools, even if the exact data file definitions are different
from the GBM file definitions. The inheritance structure of the Data Tools allows
generalization of many of the Data Tools functionality to data from other
instruments, once the interface to those data files are defined.

The Data Tools are designed with various aspects to be plugin-like.  For example,
The Data Tools provides a ``binning`` module that contains a number of binning
algorithms for pre-binned and unbinned data.  A user is not restricted to using
only the pre-packaged binning algorithms provided, but can write their own
algorithms by following the instructions on the required inputs and outputs 
expected.  The user-defined algorithms can then be used seamlessly with the
Data Tools.  The same architecture can be leveraged to allow for user-defined
background estimation algorithms, spectral models/functions, and spectral
fitting algorithms.

.. rubric:: Future
We are working to provide the GBM Data Tools via the NASA GitHub to allow
for pull and merge requests from the community.  Planned features in an upcoming
release include a response generator and the GBM localization algorithms used
for on-board GRB triggers.

Funding through the NASA Astrophysics Data Analysis Program (ADAP) will allow 
development of a generalized set of tools to expand these data tools to other 
legacy missions with similar data, including BATSE, HETE-2, Suzaku, and others.
Planned work on these Gamma-ray Data Tools, of which the GBM Data Tools will 
form a core component, will commence in 2021.

.. rubric:: Citing
If you use the GBM Data Tools in your research and publications, we would
definitely appreciate an appropriate acknowledgment and citation! We suggest the
following BibTex:

::

 @misc{GbmDataTools,
       author = {Adam Goldstein and William H. Cleveland and Daniel Kocevski},
       title = {Fermi GBM Data Tools: v1.1.1},
       year = 2022,
       url = {https://fermi.gsfc.nasa.gov/ssc/data/analysis/gbm}
 }  


.. rubric:: Additional Resources
The Fermi Science Support Center is a fantastic resource for all things Fermi.
Specifically, for GBM, a lot of useful information about the data products can 
be found `here <https://fermi.gsfc.nasa.gov/ssc/data/access/gbm/>`_.  For 
questions, bug reports, and comments, please visit the 
`Fermi Help Desk <https://fermi.gsfc.nasa.gov/ssc/help/>`_.

.. rubric:: Acknowledgments
The Fermi GBM Data Tools were partially funded by the Fermi Guest Investigator
program (NNH18ZDA001N). Special appreciation for the volunteer testers.


Table of Contents
=================
.. toctree::
   :maxdepth: 2
   
   install
   notebooks/index
   api/index
   changelog/index
 

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

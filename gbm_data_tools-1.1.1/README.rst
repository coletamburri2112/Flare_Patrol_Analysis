The Fermi GBM Data Tools
========================

.. image:: doc/images/gdt_logo.png
   :scale: 10%

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
       title = {Fermi GBM Data Tools: v1.1.0},
       year = 2021,
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


Installation
============

..  Note:: Requires: Python >=3.5
            
           Tested on:
           
           * macOS El Capitan (10.11.6) - Catalina (10.15.3)
           
           * Ubuntu 16.04 - 18.04
           
           * Windows 10 Subsystem for Linux


How to Install
--------------

The Data Tools can currently be installed from the tar file available at
the `Fermi Science Support Center
<https://fermi.gsfc.nasa.gov/ssc/data/analysis/gbm/>`_.  Download the tar file
to some place on the machine you want it installed on.

One of the easiest and suggested methods of install is to create a virtual 
environment, which allows you to install the Data Tools requirements without
worrying about possible conflicts with other versions of packages you may 
already have installed.  It's not mandatory to install in a virtual environment,
but it could make your life much easier.

.. Note::
    One of the requirements for the Data Tools is Basemap, which is a part of 
    the Matplotlib toolkit.  Basemap requires the installation of the GEOS C
    library.  If you do not have the library already compiled in a standard
    place, the Data Tools installation will attempt to compile it
    for you.  There may be cases where this fails.  If you encounter such a
    case, you should follow the instructions on how to compile the GEOS library
    for basemap `here <https://matplotlib.org/basemap/users/installing.html#installation>`_.

To install via ``pip``, navigate to the directory where you want the virtual
environment to live and::
    
    $ python3 -m venv gbm
    $ source gbm/bin/activate
    
This will create a new virtual environment called "gbm" and activates the 
environment.

To install via ``conda``:: 
    
    $ conda create --name gbm
    $ source activate gbm

This will similarly create a new virtual environment and activates it.

It's also good to ensure you have ``pip`` updated::
    
    $ pip3 install --upgrade pip
    
The complete version of the Data Tools requires Matplotlib Basemap, although
the Data Tools can be installed without Basemap with the loss of only some 
plotting functionality (EarthPlot and some SkyPlot functions). If you wish to 
install the "light" version with Basemap, simply::

    $ pip3 install <path_to_tar>/gbm_data_tools-1.1.0.tar.gz

If you want the full install, including Basemap, you will need the GEOS C
library. If it is not already installed in a standard location, set your 
``GEOS_DIR`` environment variable to a path where you would like it installed 
(and where you have write privileges).  For example, if you use bash, you can 
do something like this: ::
    
    $ export GEOS_DIR=~/.geos_dir

If you already have the GEOS C library installed in a *non-standard* location, 
set ``GEOS_DIR`` to that path instead. Then you can install the data tools::
    
    $ pip3 install <path_to_tar>/gbm_data_tools-1.1.0.tar.gz[basemap]

Finally, if you want to run the Jupyter notebooks provided with the Data Tools 
(trust us, you do)::

    $ pip3 install ipython
    $ pip3 install notebook

If you don't wish to install via a virtual environment, you are welcome to 
install with your preferred method, but you may encounter more difficulties, 
especially if you have existing package installs that conflict with what is 
required for the Data Tools.  Check the requirements.txt for required packages 
and versions.

----

How to Uninstall
----------------

To uninstall::

    $ pip3 uninstall gbm-data-tools

There are also a number of files (documentation, notebooks, etc.) for the tools
that are copied into your ``$HOME/.gbm_data_tools`` directory.  You can delete 
thesefiles if you wish.

Documentation 
-------------
On successful installation, you can launch the local HTML documentation by
calling::

    $ gbm-docs

from the command line.

----

.. _notebook_launch:

Launching the Notebooks 
-----------------------
If you have installed Jupyter as suggested above, you can run the notebooks
provided with the Data Tools.  After successful installation, the notebooks
can be launched by calling::

    $ gbm-demos

from the command line.

----

Quickstart
----------
To load the Data Tools within your python environment, simply::
    
    import gbm

The Data Tools has several different modules. For example, the data module
containing the interfaces to the GBM data files can be loaded by::

    from gbm import data
    

Known Issues
------------
* **On install, Python complains the tar.gz file is not a valid file.** 
  This is a known issue with certain browsers, where they download *and*
  decompress the file without removing the .gz extension.  The easiest fix is
  to rename the file without the .gz extension and retry the install.


* **When running a notebook in Linux, you observe a similar error**::
    
    %matplotlib notebook                                                               
    Warning: Cannot change to a different GUI toolkit: notebook. Using osx instead.
  This is due to some backend plotting issue with Jupyter notebook on Linux.
  Remove the ``%matplotlib inline`` in the notebook cell and re-evaluate the
  cell *twice* to see the plot.


* **The FTP Data Finders fail with a connection error.**  This may be due to the
  underlying OpenSSL library that the Python ``ftplib`` uses.  You may need to
  update your OpenSSL library to get this to work.  Note that this only appears
  to be an issue with the FTPS protocol, not the normal FTP protocol.  
  A potential solution is to try the following:
        
        * $ pip3 install pyopenssl
        * $ pip3 install requests[security]
  
* **The virtual environment is using your system ipython (or other package) 
  install.**  This can sometimes happen if you didn't install ipython (or other
  package) in the virtual environment.  Try installing ipython (or other package) 
  and restart your virtual environment.

* **You observe the following error**::
    
    ImportError: No module named '_tkinter'
  This is a situation where Matplotlib is using the ``tkinter`` backend for
  plotting.  You would see this error if you don't have ``tkinter`` installed. 
  You don't need to install ``tkinter`` if you don't want to; instead, you can
  create a file named `matplotlibrc` in your working directory that contains the
  following::
    backend : Agg
.. LWA Software Library documentation master file, created by
   sphinx-quickstart on Wed Dec  1 16:06:23 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

LWA Software Library Documentation
==================================
This package contains a collection of tools for reading, format shifting, 
and analyzing LWA data.  It also contains the majority of the lwa_user
library developed by NRL.

Contents
========

Installation
------------
.. toctree::
   :maxdepth: 2

   install
   rnotes

   lwa_user

General Utilities
-----------------
.. toctree::
   :maxdepth: 2

   astro
   kurtosis
   dedispersion
   catalog
   transform
   mathutils
   robust
   paths
   busy
   color
   progress

LWA-Specific Utilities
----------------------
.. toctree::
   :maxdepth: 2

   observations
   swarm

   station
   sdf
   idf
   readers
   writers

   spectra

   correlator

   imaging

   beamformer

   fakedata

Miscellaneous
-------------
.. toctree::
   :maxdepth: 2

   config
   misc

   genindex
   modindex
   search

Included Scripts
----------------
.. toctree::
   :maxdepth: 2

   scripts

Tutorials
=========
* `General LSL Routines <https://nbviewer.jupyter.org/github/lwa-project/lsl/blob/master/doc/notebooks/General%20Tutorial.ipynb>`_
* `Building SDFs Programmatically <https://nbviewer.jupyter.org/github/lwa-project/lsl/blob/master/doc/notebooks/SDF%20Generation.ipynb>`_
* `Observation Metadata <https://nbviewer.jupyter.org/github/lwa-project/lsl/blob/master/doc/notebooks/MCS%20Metadata%20Tutorial.ipynb>`_
* `Working with Data <https://nbviewer.jupyter.org/github/lwa-project/lsl/blob/master/doc/notebooks/LDP%20Tutorial.ipynb>`_

Advanced
========

.. toctree::
   :maxdepth: 2

   data-formats
   c-extensions

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

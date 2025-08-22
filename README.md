LSL - The LWA Software Library
==============================

[![GHA](https://github.com/lwa-project/lsl/actions/workflows/main.yml/badge.svg)](https://github.com/lwa-project/lsl/actions/workflows/main.yml)  [![Coverage Status](https://codecov.io/gh/lwa-project/lsl/branch/master/graph/badge.svg?token=CNCKDKBX9R)](https://codecov.io/gh/lwa-project/lsl)
  [![Documentation Status](https://readthedocs.org/projects/lsl/badge/?version=latest)](https://lsl.readthedocs.io/en/latest/?badge=latest)

[![Paper](https://img.shields.io/badge/arXiv-1209.1576-blue.svg)](https://arxiv.org/abs/1209.1576)    [![PyPi](https://img.shields.io/pypi/v/lsl.svg)](https://pypi.org/project/lsl/)

DESCRIPTION
-----------
This package contains a collection of tools for reading, format shifting, and analyzing LWA data.  It also contains the majority of the lwa_user library developed by NRL.  The lwa_user modules includes (and their respective names in the LSL package are):
 * lwa_user.astro           (lsl.astro)
 * lwa_user.catalog         (lsl.catalog)
 * lwa_user.mathutil        (lsl.misc.mathutil)
 * lwa_user.nec_util        (lsl.sim.nec_util)
 * lwa_user.skymap          (lsl.skymap)
 * lwa_user.transform       (lsl.transform)

REQUIREMENTS
------------
 * python >= 3.8
 * C++ compiler with at least C++11 support
 * fftw3 >= 3.2 (single precision version)
 * gsl >= 2.4
 * gdbm >= 1.8
 * numpy >= 1.7
 * scipy >= 0.19
 * astropy >= 5.2
 * pyephem >= 3.7.5.3
 * aipy >= 3.0.1
 * pytz >= 2012c
 * h5py
 * matplotlib >= 1.1 (required for some of the scripts)
 * BeautifulSoup (required for some of the scripts)
 * casacore (required for measurement set support)

INSTALLING
----------
The LSL package is installed as a regular Python package using distutils.  Unzip and untar the source distribution. Setup the Python interpreter you wish to use for running the package applications and switch to the root of the source distribution tree.

Install LSL by running:

	pip install . [--root=<prefix>|--user]

If the '--root' option is not provided, then the installation tree root directory is the same as for the Python interpreter used to run `setup.py`.  For instance, if the Python interpreter is in '/usr/local/bin/python', then '<prefix>' will be set to '/usr/local'.  Otherwise, the explicit \<prefix> value is taken from the command line option.  The package will install files in the following locations:
 * <prefix>/bin
 * <prefix>/lib/python3.8/site-packages
 * <prefix>/share/doc
 * <prefix>/share/install

If an alternate '\<prefix>' value is provided, you should set the PATH environment to include directory '\<prefix>/bin' and the PYTHONPATH environment to include directory '\<prefix>/lib/python3.8/site-packages'.

If the '--user' option is provided, then then installation tree root directory will be in the current user's home directory.

UNIT TESTS
----------
Unit tests for the package may be found in the 'lsl/tests' sub-directory in the package source distribution tree.  To run the complete suite of package unit tests:

    cd tests
    python -m unittest discover

TELEMETRY
---------
By default LSL installs with basic telemetry enabled in order to help inform the LWA how the software is used and to help inform future 
development.  The data collected as part of this consist seven things:
 * a timestamp for when the report is generated,
 * a unique installation identifier,
 * the LSL version being used, 
 * the Python version being used,
 * the execution time of the Python process that imports LSL,
 * which LSL modules are imported,
 * which LSL functions are used and their average execution times, and
 * which LSL scripts are used.
These data are sent to the LWA using a HTTP POST request where they are aggregated.

Users can opt out of telemetry collection using the lsl.misc.telemetry module via:

    python -m lsl.misc.telemetry --disable

This command will set a disk-based flag that disables the reporting process.  This script can also be used to re-enable telemetry and check the unique installation identifier being used.

DOCUMENTATION
-------------
See docs/README.md for documentation information.

RELEASE NOTES
-------------
See the CHANGELOG for a detailed list of changes and notes.

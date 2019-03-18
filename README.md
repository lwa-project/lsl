LSL - The LWA Software Library
==============================

### [![Paper](https://img.shields.io/badge/arXiv-1209.1576-blue.svg)](https://arxiv.org/abs/1209.1576)

DESCRIPTION
-----------
This package contains a collection of tools for reading, format shifting, and analyzing LWA data.  It also contains the majority of the lwa_user library developed by NRL.  The lwa_user modules includes (and their respective names in the LSL package are):
 * lwa_user.astro			(lsl.astro)
 * lwa_user.catalog		(lsl.catalog)
 * lwa_user.mathutil		(lsl.misc.mathutil)
 * lwa_user.nec_util		(lsl.sim.nec_util)
 * lwa_user.skymap			(lsl.skymap)
 * lwa_user.transform		(lsl.transform)

REQUIREMENTS
------------
 * python >= 2.6 and python < 3.0
 * atlas >= 3.6
 * fftw3 >= 3.2
 * gdbm >= 1.8
 * numpy >= 1.2
 * scipy >= 0.7
 * pyfits >= 3.1
 * pyephem >= 3.7.5
 * aipy >= 1.0
 * pytz >= 2011k
 * matplotlib >= 1.1 (required for some of the scripts)
 * BeautifulSoup (required for some of the scripts)
 * casacore (required for measurement set support)

BUILDING
--------
The LSL package is installed as a regular Python package using distutils.  Unzip and untar the source distribution. Setup the Python interpreter you wish to use for running the package applications and switch to the root of the source distribution tree.

A setup configuration script is required to run the build process.  A sample config file may be found in the root directory of the source distribution as `setup.cfg`.  Do not modify any of the lines in the `[build_ext]` section of the config file.

To build the LSL package, run:

    python setup.py build

TESTING
-------
To test the as-build LSL package, run:

	python setup.py test

INSTALLATION
------------
Install LSL by running:

	pip install . [--root=<prefix>|--user]

If the '--root' option is not provided, then the installation tree root directory is the same as for the Python interpreter used to run `setup.py`.  For instance, if the Python interpreter is in '/usr/local/bin/python', then '<prefix>' will be set to '/usr/local'.  Otherwise, the explicit <prefix> value is taken from the command line option.  The package will install files in the following locations:
 * <prefix>/bin
 * <prefix>/lib/python2.6/site-packages
 * <prefix>/share/doc
 * <prefix>/share/install

If an alternate '<prefix>' value is provided, you should set the PATH environment to include directory '<prefix>/bin' and the PYTHONPATH environment to include directory '<prefix>/lib/python2.6/site-packages'.

If the '--user' option is provided, then then installation tree root directory will be in the current user's home directory.

UNIT TESTS
----------
Unit tests for the package may be found in the 'lsl/tests' sub-directory in the package source distribution tree.  To run the complete suite of package unit tests:

    python test_lsl.py
    
The tests may also be run from an interactive Python session after LSL has been successfully installed via: 

	import lsl
	lsl.test()

DOCUMENTATION
-------------
See doc/README for documentation information.

RELEASE NOTES
-------------
See the CHANGELOG for a detailed list of changes and notes.

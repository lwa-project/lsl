Requirements
============
 * python >= 3.8
 * fftw3 >= 3.2 (single precision version)
 * gsl >= 2.4
 * gdbm >= 1.8
 * numpy >= 1.7
 * scipy >= 0.19
 * astropy >= 5.2
 * ephem >= 3.7.5.3
 * aipy >= 3.0.1
 * pytz >= 2012c
 * matplotlib >= 1.1 [1]_
 * BeautifulSoup [1]_
 * casacore [2]_

.. [1] Required for some of the included scripts
.. [2] Required for measurement set support

Installing
==========
The LSL package is installed as a regular Python package using distutils.  
Unzip and untar the source distribution. Setup the python interpreter you 
wish to use for running the package applications and switch to the root of 
the source distribution tree.

Install LSL by running::
	
	pip install . [--root=<prefix>|--user]

If the '--root' option is not provided, then the installation 
tree root directory is the same as for the python interpreter used to run 
setup.py.  For instance, if the python interpreter is in 
'/usr/local/bin/python', then <prefix> will be set to '/usr/local'.
Otherwise, the explicit <prefix> value is taken from the command line
option.  The package will install files in the following locations:
 * <prefix>/bin
 * <prefix>/lib/python3.6/site-packages
 * <prefix>/share/doc
 * <prefix>/share/install

If an alternate <prefix> value is provided, you should set the PATH
environment to include directory '<prefix>/bin' and the PYTHONPATH
environment to include directory '<prefix>/lib/python3.6/site-packages'.

If the '--user' option is provided, then then installation tree root 
directory will be in the current user's home directory.

Unit Tests
==========
Unit tests for the package may be found in the 'lsl/tests' sub-directory in the package source distribution tree.  To run the complete suite of package unit tests::

    cd tests
    python -m unittest discover

Telemetry
=========
By default LSL installs with basic telemetry enabled in order to help
inform the LWA how the software is used and to help inform future 
development.  The data collected as part of this consist seven things:
 * a timestamp for when the report is generated,
 * a unique installation identifier,
 * the LSL version being used, 
 * the execution time of the Python process that imports LSL,
 * which LSL modules are imported,
 * which LSL functions are used and their average execution times, and
 * which LSL scripts are used.

These data are sent to the LWA using a HTTP POST request where they
are aggregated.

Users can opt out of telemetry collection using the lsl.misc.telemetry module
via::

    python -m lsl.misc.telemetry --disable

This command will set a disk-based flag that disables the reporting process.
This script can also be used to re-enable telemetry and check the unique
installation identifier being used.

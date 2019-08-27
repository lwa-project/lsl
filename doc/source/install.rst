Requirements
============
 * python >= 2.7
 * fftw3 >= 3.2
 * gdbm >= 1.8
 * numpy >= 1.7
 * scipy >= 0.19
 * astropy < 2.0 (2.7) or >= 3.0 (3.0+)
 * ephem >= 3.7.5.3
 * aipy >= 3.0.1
 * pytz >= 2012c
 * matplotlib >= 1.1 [1]_
 * BeautifulSoup [1]_
 * casacore [2]_

.. [1] Required for some of the included scripts
.. [2] Required for measurement set support

Building
========
The LSL package is installed as a regular Python package using distutils.  
Unzip and untar the source distribution. Setup the python interpreter you 
wish to use for running the package applications and switch to the root of 
the source distribution tree.

A setup configuration script is required to run the build process.  A sample 
config file may be found in the root directory of the source distribution as 
'setup.cfg'.  Do not modify any of the lines in the [build_ext] section of 
the config file.

To build the LSL package, run::

	python setup.py build

Testing
=======
To test the as-build LSL package, run::

	python setup.py test

Installing
==========
Install LSL by running::
	
	pip install . [--root=<prefix>|--user]

If the '--root' option is not provided, then the installation 
tree root directory is the same as for the python interpreter used to run 
setup.py.  For instance, if the python interpreter is in 
'/usr/local/bin/python', then <prefix> will be set to '/usr/local'.
Otherwise, the explicit <prefix> value is taken from the command line
option.  The package will install files in the following locations:
 * <prefix>/bin
 * <prefix>/lib/python2.6/site-packages
 * <prefix>/share/doc
 * <prefix>/share/install
If an alternate <prefix> value is provided, you should set the PATH
environment to include directory '<prefix>/bin' and the PYTHONPATH
environment to include directory '<prefix>/lib/python2.6/site-packages'.

If the '--user' option is provided, then then installation tree root 
directory will be in the current user's home directory.

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

Users can opt out of telemetry collection using the provided lslTelemetry.py
script via:

    python lslTelemetry.py --disable

This command will set a disk-based flag that disables the reporting process.
This script can also be used to re-enable telemetry and check the unique
installation identifier being used.
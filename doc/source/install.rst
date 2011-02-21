Requirements
============
 * python >= 2.6 and python < 3.0
 * swig >= 1.3
 * numpy >= 1.2
 * scipy >= 0.7
 * pyfits >= 2.1
 * ephem >= 3.7.3
 * aipy >= 0.9.1
 * matplotlib >= 0.98.3 [#]_

.. [#] Required for some of the included scripts

In addition, the libnova C library is required.  This is an astronomical
calculation and ephemeris library.  The lsl.libnova python module is built 
as a SWIG wrapper around the libnova library, and this provides the basis 
of the astronomical code contained in the package.  If you are using binary 
packages, be sure to install the libnova-dev version, as the library header 
include files are required for the build.  Some issues have arisen with the 
system version of swig for some platforms: this default version is too old 
to correctly generate wrapper code for python 2.6.

.. warning::
	Use libnova version 0.13.0 or later; version 0.12 is known to contain errors.

The libnova C library can be found `here <http://libnova.sourceforge.net/>`_.

Building
========
The LSL package is installed as a regular Python package using distutils.  
Unzip and untar the source distribution. Setup the python interpreter you 
wish to use for running the package applications and switch to the root of 
the source distribution tree.

A setup configuration script is required to run the build process.  A sample 
config file may be found in the root directory of the source distribution as 
'setup.cfg'.  This file must be edited to specify the location of the 
libnova library files if that location is not in '/usr/local'.  Edit the 
following lines::

	***********************************************************************
	[DEFAULT]
	
	# The libnova_prefix value should point to the platform installation of
	# the libnova libraries and header files
		
	libnova_prefix = /usr/local
	***********************************************************************
	and change to to the location where libnova is installed on your machine:
	***********************************************************************
	[DEFAULT]
		
	# The libnova_prefix value should point to the platform installation of
	# the libnova libraries and header files
		
	libnova_prefix = /Software/libnova/0.13.0
	***********************************************************************

The setup process will expect to find the libnova library C header files in 
'<libnova_prefix>/include/libnova' and the share library object file in 
'<libnova_prefix>/lib'.  Do not modify any of the lines in the [build_ext] 
section of the config file.

To build the LSL package, run::

	python setup.py build

Testing
=======
To test the as-build LSL package, run::

	python setup.py test

Installing
==========
Install LSL by running::
	
	python setup.py install [--prefix=<prefix>|--user]

If the '--prefix' option is not provided, then the installation 
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

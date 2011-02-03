#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Given a GSM output model in the HEALpix RING format, reduce the resolution
of the map to a more sane value and save the data along with the RA and Dec
coordinates into a numpy NPZ file."""

import sys
import ephem
import numpy
import healpy
import getopt


def sides2res(nSide):
	"""Convert a value for NSide to a mean HEALpix pixel spacing in degrees."""

	return 58.6323 / nSide


def sides2area(nSide):
	"""Convert a vale for NSide to a HEALpix pixel area in steradians."""

	nPix = 12*nSide*nSide
	return 4*numpy.pi / nPix


def usage(exitCode=None):
	print """gsm2npz.py - Read in a GSM RING-style HEALpix file created by gsm
and convert it to a NPZ file with RA and dec. pairs at a lower resolution.

Usage: gsm2npz.py [OPTIONS] gsm_file [gsm_file [...]]

Options:
-h, --help              Display this help information
-n, --nside             Number of sides in output NPZ file (default = 32)
-f, --full-sky          Save all pixels, not just those visible from LWA-1
                        (default = save only visible)
                        
Mean Resolution as a Function of NSide
   1: 58.632 degrees
   2: 29.316 degrees
   4: 14.658 degrees
   8:  7.329 degrees
  16:  3.665 degrees
  32:  1.832 degrees
  64:  0.916 degrees
 128:  0.458 degrees
 256:  0.229 degrees
 512:  0.115 degrees
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	config['nSide'] = 32
	config['fullSky'] = False
	config['args'] = []

	try:
		opts, args = getopt.getopt(args, "hfn:", ["help", "full-sky", "nside="])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)

	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-f', '--full-sky'):
			config['fullSky'] = True
		elif opt in ('-n', '--nside'):
			config['nSide'] = int(value)
		else:
			assert False

	# Add in arguments
	config['args'] = args

	# Return configuration
	return config


def main(args):
	config = parseOptions(args)

	gsmFilenames = config['args']
	for gsmFilename in gsmFilenames:
		print "Working on GSM file '%s'" % gsmFilename
		# Load in the file using the numpy.loadtxt routine and figure out how many 
		# sides are present in the file
		hMap = numpy.loadtxt(gsmFilename)
		nSide = healpy.pixelfunc.npix2nside(len(hMap))
		print "-> Native number of sides: %i (~%.4f degrees)" % (nSide, sides2res(nSide))

		# Downgrade the map to a lower resolution
		map = healpy.pixelfunc.ud_grade(hMap, config['nSide'])
		nSide = healpy.pixelfunc.npix2nside(len(map))
		print "-> Lower-resolution number of sides: %i (~%.4f degrees) " % (nSide, sides2res(nSide))


		# Initialize the array for holding the equatorial coordinates and run with it
		eCoords = numpy.zeros((len(map),2), dtype=numpy.float_)
		i = 0
		for row in eCoords:
			t = healpy.pixelfunc.pix2ang(nSide, i)
			u = ephem.Equatorial(ephem.Galactic(t[1], numpy.pi/2-t[0], epoch='2000'))
			eCoords[i,0] = float(u.ra)
			eCoords[i,1] = float(u.dec)

			i = i + 1

		# Save the output
		outFilename = gsmFilename.replace('.dat', '.npz')
		if config['fullSky']:
			numpy.savez(outFilename, nSide=nSide, pixRes=sides2res(nSide), pixArea=sides2area(nSide), eCoords=eCoords, T=map)
			print "-> Saved %i points to '%s'" % (len(map), outFilename)
		else:
			valid = numpy.where(eCoords[:,1] >= -60*numpy.pi/180.0)[0]
			numpy.savez(outFilename, nSide=nSide, pixRes=sides2res(nSide), pixArea=sides2area(nSide), eCoords=eCoords[valid,:], T=map[valid])
			print "-> Saved %i points to '%s'" % (len(valid), outFilename)


if __name__ == "__main__":
	main(sys.argv[1:])

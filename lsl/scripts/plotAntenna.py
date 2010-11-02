#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Example script to plot the relative response of an isolated LWA antenna 
as a function of azimuth and elevation using an NEC model."""

import os
import sys
import shutil
import tempfile
import numpy as n
from scipy.special import sph_harm

from lsl.common.paths import data as dataPath
from lsl.sim import nec_util

import matplotlib.pyplot as plt

def main(args):
	if len(args) > 0:
		freq = float(args[0])
	else:
		freq = 50.0

	# Create a temporary directory to hold the model that we are
	# going to work on
	modelPath = tempfile.mkdtemp(prefix='plotAntenna-', suffix='.tmp')
	shutil.copy(os.path.join(dataPath, 'lwa-antenna-model.nec'), modelPath)

	# Run the NEC model and extract the results
	nec_file = os.path.join(modelPath, 'lwa-antenna-model.nec')
	pattern_dB = nec_util.NECPattern(nec_file, freq).antenna_pat_dB 
	pattern = 10**(pattern_dB/10)
	pattern /= pattern.max()

	# Plot!
	fig = plt.figure()
	ax = fig.add_subplot(1, 1, 1)
	p = ax.imshow(pattern.T, origin='lower')
	
	cb = fig.colorbar(p, ax=ax)
	cb.ax.set_ylabel('Relative Response')

	ax.set_title('LWA Antenna Model @ %.1f MHz' % freq)
	ax.set_xlabel('Azimuth [deg.]')
	ax.set_ylabel('Elevation [deg.]')
	plt.show()


if __name__ == "__main__":
	main(sys.argv[1:])

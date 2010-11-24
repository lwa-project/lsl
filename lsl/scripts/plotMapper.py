#!/usr/bin/env python

import sys
import numpy
import pyfits

from lsl.correlator import uvUtils

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

hdu = pyfits.open(sys.argv[1])
mapper = hdu['NOSTA_MAPPER']

nosta = mapper.data.field('NOSTA')
noact = mapper.data.field('NOACT')

fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)
ax2 = plt.axes([0.75, 0.75, 0.1, 0.1])

xyz = uvUtils.getXYZ(numpy.arange(1,259))
c = ax1.scatter(xyz[:,0], xyz[:,1], c=xyz[:,2], s=40.0, alpha=0.50)
ax1.plot(xyz[noact-1,0], xyz[noact-1,1], marker='x', linestyle=' ', alpha=1.0, color='k')
for m,a in zip(nosta, noact):
	ax1.text(xyz[a-1,0], xyz[a-1,1], '%i' % m)
ax1.set_xlabel('$\Delta$X [E-W; m]')
ax1.set_xlim([-80, 80])
ax1.set_ylabel('$\Delta$Y [N-S; m]')
ax1.set_ylim([-80, 80])

ax2.scatter(xyz[:,0], xyz[:,1], c=xyz[:,2], s=40.0, alpha=0.50)
ax2.plot(xyz[noact-1,0], xyz[noact-1,1], marker='x', linestyle=' ', alpha=1.0, color='k')
for m,a in zip(nosta, noact):
	ax2.text(xyz[a-1,0], xyz[a-1,1], '%i' % m)
ax2.set_title('Outlier')
ax2.set_xlim([335, 345])
ax2.set_ylim([10, 20])
ax2.xaxis.set_major_formatter( NullFormatter() )
ax2.yaxis.set_major_formatter( NullFormatter() )

plt.show()




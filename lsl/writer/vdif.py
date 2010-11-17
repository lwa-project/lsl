import math
import numpy
import ephem
import struct

from lsl.common import stations as lwa_common
from lsl.common import dp as dp_common
from lsl.reader import tbw

vdifEpoch = ephem.Date('2000/01/01 00:00:00.00')
unixEpoch = ephem.Date('1970/01/01 00:00:00.00')


def log2(numb):
	"""Calculate the log-base 2 of a number or array."""
  
	return numpy.log(numb) / numpy.log(2.0)


class VDIFFrame(object):
	def __init__(self, stand=0, time=0, frame=0, bits=16, sampleRate=dp_common.fS, data=None):
		self.stand = stand
		self.time = time
		self.frame = frame
		self.bits = bits
		self.sampleRate = sampleRate
		self.data = data - data.min()

		self.__setEpoch()

		
	def __setEpoch(self):
		curEpoch = float(unixEpoch) + self.time
		self.epoch = 0
		self.seconds = long(curEpoch - float(vdifEpoch))

	def isValid(self):
		logChan = log2(len(self.data))
		if int(logChan) != logChan:
			return False
		return True

	def createRawFrame(self):
		samplesPerWord = 32 / self.bits
		raw = numpy.zeros(32 + samplesPerWord*len(self.data), dtype=numpy.uint8)

		# Valid data, standard 32-bit header, and seconds since the epoch
		raw[3] = (0 << 7) | (0 << 6) | ((self.seconds >> 24) & 63)
		raw[2] = (self.seconds >> 16) & 255
		raw[1] = (self.seconds >> 8) & 255
		raw[0] = self.seconds & 255

		# Refernce epoch and frame count
		raw[7] = self.epoch & 63
		raw[6] = (self.frame >> 16) & 255
		raw[5] = (self.frame >> 8) & 255
		raw[4] = self.frame & 255

		# VDIF version number, number of channels, and data frame length
		raw[11] = (3 << 6) | int(math.log10(len(data))/math.log10(2.0)) & 31
		raw[10] = (len(raw) >> 16) & 255
		raw[9] = (len(raw) >> 8) & 255
		raw[8] = len(raw) & 255

		# Data type, bits per sample, thread ID, and station ID
		# NB:  thread ID is fixed at `0'
		if self.data.dtype.kind == 'c':
			raw[15] = (1 << 7) | ((self.bits & 31) << 2) | ((0 >> 8) & 3)
		else:
			raw[15] = (0 << 7) | ((self.bits & 31) << 2) | ((0 >> 8) & 3)
		raw[14] = 0 & 255
		raw[13] = (self.stand >> 8) & 255
		raw[12] = self.stand & 255
		
		# User-defined words 4 through 7 (array entries 16 to 31)
		
		# Data values
		if self.data.dtype.kind == 'c':
			for f in range(32,len(raw),4):
				i = (f - 32) / 4
				word = 0L
				for p in range(32/self.bits/2):
					word |= (self.data[f*32/self.bits/2 + p].real << self.bits*p)
					word |= (self.data[f*32/self.bits/2 + p].imag << self.bit*(p+1))
				raw[f+3] = (word >> 24) & 255
				raw[f+2] = (word >> 16) & 255
				raw[f+1] = (word >> 8) & 255
				raw[f+0] = word & 255
		else:
			for f in range(32,len(raw),4):
				i = (f - 32) / 4
				word = 0L
				for p in range(32/self.bits):
					print f, p, i*32/self.bits+p
					word |= (self.data[i*32/self.bits + p] << self.bits*p)
				raw[f+3] = (word >> 24) & 255
				raw[f+2] = (word >> 16) & 255
				raw[f+1] = (word >> 8) & 255
				raw[f+0] = word & 255
				
		return raw

	def writeRawFrame(self, fh):
		print self.isValid()
		rawFrame = self.createRawFrame()
		rawFrame.tofile(fh)


if __name__ == "__main__":
	fh = open('/home/jayce/TBW Data/multiTBW_Sept_19_8pm.dat', 'rb')
	
	# Make sure that the data is TBW and determine the data length
	test = tbw.readFrame(fh)
	print "TBW Data:  %s" % test.header.isTBW()
	if not test.header.isTBW():
		raise errors.notTBWError()
	print "Data Length: %i bits" % test.getDataBits()
	if test.header.getDataBits() == 12:
		nData = 400
	else:
		nData = 1200
	fh.seek(0)

	# Read in the data and add it to the FITS file created above
	count = {}
	syncCount = 0
	masterCount = 0
	data = numpy.zeros((20,12000000), dtype=numpy.int16)
	for i in range(3000):
		# Read in the next frame and anticipate any problems that could occur
		try:
			cFrame = tbw.readFrame(fh, Verbose=False)
		except errors.eofError:
			break
		except errors.syncError:
			#print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/TBWFrameSize-1)
			syncCount = syncCount + 1
			continue
		except errors.numpyError:
			break

		stand = cFrame.header.parseID()
		if cFrame.header.frameCount % 10000 == 0:
			print "%2i  %14i  %6.3f  %5i  %5i" % (stand, cFrame.data.timeTag, cFrame.getTime(), cFrame.header.frameCount, cFrame.header.secondsCount)
		if stand not in count.keys():
			count[stand] = 0

		data[2*stand, count[stand]*nData:(count[stand]+1)*nData] = cFrame.data.xy[0,:]
		data[2*stand+1, count[stand]*nData:(count[stand]+1)*nData] = cFrame.data.xy[1,:]

		count[stand] = count[stand] + 1
		masterCount = masterCount + 1
	fh.close()
	
	fh = open('test.vdif', 'wb')
	for i in range(1):
		frame = VDIFFrame(stand=0, time=i, frame=1, data=data[0,0:512])
		frame.writeRawFrame(fh)
	fh.close()
	
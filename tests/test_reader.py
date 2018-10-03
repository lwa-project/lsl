# -*- coding: utf-8 -*-

# Python3 compatiability
from __future__ import division
import sys
if sys.version_info > (3,):
    xrange = range
    
"""Unit test for lsl.reader modules"""

import os
import unittest

from lsl.common.paths import DATA_BUILD
from lsl.reader import tbw
from lsl.reader import tbn
from lsl.reader import drx
from lsl.reader import vdif
from lsl.reader import drspec
from lsl.reader import errors


__revision__ = "$Rev$"
__version__  = "0.6"
__author__    = "Jayce Dowell"


tbwFile = os.path.join(DATA_BUILD, 'tests', 'tbw-test.dat')
tbnFile = os.path.join(DATA_BUILD, 'tests', 'tbn-test.dat')
drxFile = os.path.join(DATA_BUILD, 'tests', 'drx-test.dat')
vdifFile = os.path.join(DATA_BUILD, 'tests', 'vdif-test.dat')
drspecFile = os.path.join(DATA_BUILD, 'tests', 'drspec-test.dat')


class reader_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.reader
    modules."""
    
    ### TBW ###
    
    def test_tbw_read(self):
        """Test reading in a frame from a TBW file."""
        
        fh = open(tbwFile, 'rb')
        # First frame is really TBW and stores the correct stand ID
        frame1 = tbw.read_frame(fh)
        self.assertTrue(frame1.header.is_tbw)
        self.assertEqual(frame1.id, 2)
        # Second frame
        frame2 = tbw.read_frame(fh)
        self.assertTrue(frame2.header.is_tbw)
        self.assertEqual(frame2.id, 1)
        fh.close()
        
    def test_tbw_bits(self):
        """Test getting the data bits from a TBW file."""
        
        fh = open(tbwFile, 'rb')
        # File contains 12-bit data, two ways
        self.assertEqual(tbw.get_data_bits(fh), 12)
        frame1 = tbw.read_frame(fh)
        self.assertEqual(frame1.data_bits, 12)
        fh.close()
        
    def test_tbw_errors(self):
        """Test reading errors."""
        
        fh = open(tbwFile, 'rb')
        # Frames 1 through 8
        for i in range(1,9):
            frame = tbw.read_frame(fh)
            
        # Last frame should be an error (errors.EOFError)
        self.assertRaises(errors.EOFError, tbw.read_frame, fh)
        fh.close()
        
        # If we offset in the file by 1 byte, we should be a 
        # sync error (errors.SyncError).
        fh = open(tbwFile, 'rb')
        fh.seek(1)
        self.assertRaises(errors.SyncError, tbw.read_frame, fh)
        fh.close()
        
    def test_tbw_comps(self):
        """Test the TBW frame comparison operators (>, <, etc.) for time tags."""
        
        fh = open(tbwFile, 'rb')
        # Frames 1 through 3
        frames = []
        for i in range(1,4):
            frames.append(tbw.read_frame(fh))
        fh.close()
        
        self.assertTrue(0 < frames[0])
        self.assertFalse(0 > frames[0])
        self.assertTrue(frames[-1] >= frames[0])
        self.assertFalse(frames[-1] <= frames[0])
        self.assertTrue(frames[0] == frames[0])
        self.assertFalse(frames[0] == frames[-1])
        self.assertFalse(frames[0] != frames[0])
        
    def test_tbw_sort(self):
        """Test sorting TBW frames by time tags."""
        
        fh = open(tbwFile, 'rb')
        # Frames 1 through 3
        frames = []
        for i in range(1,4):
            frames.append(tbw.read_frame(fh))
        fh.close()
        
        frames.sort()
        frames = frames[::-1]
        
        for i in xrange(1,len(frames)):
            self.assertTrue( frames[i-1] >= frames[i] )
            
    def test_tbw_math(self):
        """Test mathematical operations on TBW frame data via frames."""
        
        fh = open(tbwFile, 'rb')
        # Frames 1 through 3
        frames = []
        for i in range(1,4):
            frames.append(tbw.read_frame(fh))
        fh.close()
        
        # Multiplication
        frameT = frames[0] * 2.0
        for i in range(800):
            self.assertAlmostEqual(frameT.data.xy[i%2, i//2], 2*frames[0].data.xy[i%2, i//2], 2)
        frameT *= 2.0
        for i in range(800):
            self.assertAlmostEqual(frameT.data.xy[i%2, i//2], 4*frames[0].data.xy[i%2, i//2], 2)
        frameT = frames[0] * frames[1]
        for i in range(800):
            self.assertAlmostEqual(frameT.data.xy[i%2, i//2], frames[0].data.xy[i%2, i//2]*frames[1].data.xy[i%2, i//2], 2)
            
        # Addition
        frameA = frames[0] + 2.0
        for i in range(800):
            self.assertAlmostEqual(frameA.data.xy[i%2, i//2], 2+frames[0].data.xy[i%2, i//2], 2)
        frameA += 2.0
        for i in range(800):
            self.assertAlmostEqual(frameA.data.xy[i%2, i//2], 4+frames[0].data.xy[i%2, i//2], 2)
        frameA = frames[0] + frames[1]
        for i in range(800):
            self.assertAlmostEqual(frameA.data.xy[i%2, i//2], frames[0].data.xy[i%2, i//2]+frames[1].data.xy[i%2, i//2], 2)
            
    ### TBN ###
    
    def test_tbn_read(self):
        """Test reading in a frame from a TBN file."""
        
        fh = open(tbnFile, 'rb')
        # First frame is really TBN and stores the correct IDs
        frame1 = tbn.read_frame(fh)
        stand, pol = frame1.id
        self.assertEqual(stand, 1)
        self.assertEqual(pol, 0)
        # Second frame
        frame2 = tbn.read_frame(fh)
        stand, pol = frame2.id
        self.assertEqual(stand, 1)
        self.assertEqual(pol, 1)
        fh.close()
        
    def test_tbn_errors(self):
        """Test reading in all frames from a truncated TBN file."""
        
        fh = open(tbnFile, 'rb')
        # Frames 1 through 29
        for i in range(1,30):
            frame = tbn.read_frame(fh)
            
        # Last frame should be an error (errors.EOFError)
        self.assertRaises(errors.EOFError, tbn.read_frame, fh)
        fh.close()
        
        # If we offset in the file by 1 byte, we should be a 
        # sync error (errors.SyncError).
        fh = open(tbnFile, 'rb')
        fh.seek(1)
        self.assertRaises(errors.SyncError, tbn.read_frame, fh)
        fh.close()
        
    def test_tbn_block(self):
        """Test finding out how many stands are in a file."""
        
        fh = open(tbnFile, 'rb')
        nx, ny = tbn.get_frames_per_obs(fh)
        self.assertEqual(nx, 10)
        self.assertEqual(ny, 10)
        fh.close()
        
    def test_tbn_rate(self):
        """Test finding out the sample rate of a TBN file."""
        
        fh = open(tbnFile, 'rb')
        rate = tbn.get_sample_rate(fh)
        self.assertEqual(rate, 100000)
        code = tbn.get_sample_rate(fh, FilterCode=True)
        self.assertEqual(code, 7)
        fh.close()
        
    def test_tbn_comps(self):
        """Test the TBN frame comparison operators (>, <, etc.) for time tags."""
        
        fh = open(tbnFile, 'rb')
        # Frames 1 through 29
        frames = []
        for i in range(1,30):
            frames.append(tbn.read_frame(fh))
        fh.close()
        
        self.assertTrue(0 < frames[0])
        self.assertFalse(0 > frames[0])
        self.assertTrue(frames[-1] >= frames[0])
        self.assertFalse(frames[-1] <= frames[0])
        self.assertTrue(frames[0] == frames[0])
        self.assertFalse(frames[0] == frames[-1])
        self.assertFalse(frames[0] != frames[0])
        
    def test_tbn_sort(self):
        """Test sorting TBN frames by time tags."""
        
        fh = open(tbnFile, 'rb')
        # Frames 1 through 29
        frames = []
        for i in range(1,30):
            frames.append(tbn.read_frame(fh))
        fh.close()
        
        frames.sort()
        frames = frames[::-1]
        
        for i in xrange(1,len(frames)):
            self.assertTrue( frames[i-1] >= frames[i] )
            
    def test_tbn_math(self):
        """Test mathematical operations on TBN frame data via frames."""
        
        fh = open(tbnFile, 'rb')
        # Frames 1 through 29
        frames = []
        for i in range(1,30):
            frames.append(tbn.read_frame(fh))
        fh.close()
        
        # Multiplication
        frameT = frames[0] * 2.0
        for i in range(512):
            self.assertAlmostEqual(frameT.data.iq[i], 2*frames[0].data.iq[i], 2)
        frameT *= 2.0
        for i in range(512):
            self.assertAlmostEqual(frameT.data.iq[i], 4*frames[0].data.iq[i], 2)
        frameT = frames[0] * frames[1]
        for i in range(512):
            self.assertAlmostEqual(frameT.data.iq[i], frames[0].data.iq[i]*frames[1].data.iq[i], 2)
            
        # Addition
        frameA = frames[0] + 2.0
        for i in range(512):
            self.assertAlmostEqual(frameA.data.iq[i], 2+frames[0].data.iq[i], 2)
        frameA += 2.0
        for i in range(512):
            self.assertAlmostEqual(frameA.data.iq[i], 4+frames[0].data.iq[i], 2)
        frameA = frames[0] + frames[1]
        for i in range(512):
            self.assertAlmostEqual(frameA.data.iq[i], frames[0].data.iq[i]+frames[1].data.iq[i], 2)
            
    ### TBW/TBN Mix-up ###
    
    def test_tbw_tbn_catch(self):
        """Test that tbw will not read tbn files and vice versa."""
        
        fh = open(tbnFile, 'rb')
        frame1 = tbw.read_frame(fh)
        self.assertFalse(frame1.header.is_tbw)
        fh.close()
        
        fh = open(tbwFile, 'rb')
        frame1 = tbn.read_frame(fh)
        self.assertFalse(frame1.header.is_tbn)
        fh.close()
        
    ### DRX ###
    
    def test_drx_read(self):
        """Test reading in a frame from a DRX file."""
        
        fh = open(drxFile, 'rb')
        # First frame is really DRX and stores the IDs
        frame1 = drx.read_frame(fh)
        beam, tune, pol = frame1.id
        self.assertEqual(beam, 4)
        self.assertEqual(tune, 1)
        self.assertEqual(pol,  1)
        # Second frame
        frame2 = drx.read_frame(fh)
        beam, tune, pol = frame2.id
        self.assertEqual(beam, 4)
        self.assertEqual(tune, 2)
        self.assertEqual(pol,  0)
        fh.close()
        
    def test_drx_errors(self):
        """Test reading in all frames from a truncated DRX file."""
        
        fh = open(drxFile, 'rb')
        # Frames 1 through 32
        for i in range(1,33):
            frame = drx.read_frame(fh)
            
        # Last frame should be an error (errors.EOFError)
        self.assertRaises(errors.EOFError, drx.read_frame, fh)
        fh.close()
        
        # If we offset in the file by 1 byte, we should be a 
        # sync error (errors.SyncError).
        fh = open(drxFile, 'rb')
        fh.seek(1)
        self.assertRaises(errors.SyncError, drx.read_frame, fh)
        fh.close()
        
    def test_drx_beam(self):
        """Test finding out how many beams are present in a DRX file."""
        
        fh = open(drxFile, 'rb')
        nBeam = drx.get_beam_count(fh)
        self.assertEqual(nBeam, 1)
        fh.close()
        
    def test_drx_block(self):
        """Test finding out how many tunings/pols. per beam are in a DRX file."""
        
        fh = open(drxFile, 'rb')
        b1, b2, b3, b4 = drx.get_frames_per_obs(fh)
        self.assertEqual(b1, 0)
        self.assertEqual(b2, 0)
        self.assertEqual(b3, 0)
        self.assertEqual(b4, 4)
        fh.close()
        
    def test_drx_rate(self):
        """Test finding out the DRX sample rate."""
        
        fh = open(drxFile, 'rb')
        cFrame = drx.read_frame(fh)
        fh.seek(0)
        
        # Sample rate
        self.assertEqual(cFrame.sample_rate, drx.get_sample_rate(fh))
        
        # Filter code
        self.assertEqual(cFrame.filter_code, drx.get_sample_rate(fh, FilterCode=True))
        fh.close()
        
    def test_drx_comps(self):
        """Test the DRX frame comparison operators (>, <, etc.) for time tags."""
        
        fh = open(drxFile, 'rb')
        # Frames 1 through 10
        frames = []
        for i in range(1,11):
            frames.append(drx.read_frame(fh))
        fh.close()
        
        self.assertTrue(0 < frames[0])
        self.assertFalse(0 > frames[0])
        self.assertTrue(frames[-1] >= frames[0])
        self.assertFalse(frames[-1] <= frames[0])
        self.assertTrue(frames[0] == frames[0])
        self.assertFalse(frames[0] == frames[-1])
        self.assertFalse(frames[0] != frames[0])
        
    def test_drx_sort(self):
        """Test sorting DRX frames by time tags."""
        
        fh = open(drxFile, 'rb')
        # Frames 1 through 10
        frames = []
        for i in range(1,11):
            frames.append(drx.read_frame(fh))
        fh.close()
        
        frames.sort()
        frames = frames[::-1]
        
        for i in xrange(1,len(frames)):
            self.assertTrue( frames[i-1] >= frames[i] )
            
    def test_drx_math(self):
        """Test mathematical operations on DRX frame data via frames."""
        
        fh = open(drxFile, 'rb')
        # Frames 1 through 10
        frames = []
        for i in range(1,11):
            frames.append(drx.read_frame(fh))
        fh.close()
        
        # Multiplication
        frameT = frames[0] * 2.0
        for i in range(4096):
            self.assertAlmostEqual(frameT.data.iq[i], 2*frames[0].data.iq[i], 2)
        frameT *= 2.0
        for i in range(4096):
            self.assertAlmostEqual(frameT.data.iq[i], 4*frames[0].data.iq[i], 2)
        frameT = frames[0] * frames[1]
        for i in range(4096):
            self.assertAlmostEqual(frameT.data.iq[i], frames[0].data.iq[i]*frames[1].data.iq[i], 2)
            
        # Addition
        frameA = frames[0] + 2.0
        for i in range(4096):
            self.assertAlmostEqual(frameA.data.iq[i], 2+frames[0].data.iq[i], 2)
        frameA += 2.0
        for i in range(4096):
            self.assertAlmostEqual(frameA.data.iq[i], 4+frames[0].data.iq[i], 2)
        frameA = frames[0] + frames[1]
        for i in range(4096):
            self.assertAlmostEqual(frameA.data.iq[i], frames[0].data.iq[i]+frames[1].data.iq[i], 2)
            
    ### DR Spectrometer ###
    
    def test_drspec_read(self):
        """Test reading in a frame from a DR spectrometer file."""
        
        fh = open(drspecFile, 'rb')
        # First frame is really DR spectrometer and stores the IDs
        frame1 = drspec.read_frame(fh)
        beam = frame1.id
        self.assertEqual(beam, 1)
        
        # Second frame
        frame2 = drspec.read_frame(fh)
        beam = frame2.id
        self.assertEqual(beam, 1)
        fh.close()
        
    def test_drspec_errors(self):
        """Test reading in all frames from a truncated DR spectrometer file."""
        
        fh = open(drspecFile, 'rb')
        # Frames 1 through 8
        for i in range(1,8):
            frame = drspec.read_frame(fh)
            
        # Last frame should be an error (errors.EOFError)
        self.assertRaises(errors.EOFError, drspec.read_frame, fh)
        fh.close()
        
        # If we offset in the file by 1 byte, we should be a 
        # sync error (errors.SyncError).
        fh = open(drspecFile, 'rb')
        fh.seek(1)
        self.assertRaises(errors.SyncError, drspec.read_frame, fh)
        fh.close()
        
    def test_drspec_metadata(self):
        """Test finding out the DR spectrometer metadata."""
        
        fh = open(drspecFile, 'rb')
        cFrame = drspec.read_frame(fh)
        fh.seek(0)
        
        # Beam
        self.assertEqual(cFrame.id, 1)
        
        # Sample rate
        self.assertAlmostEqual(cFrame.sample_rate, 19.6e6, 1)
        self.assertAlmostEqual(cFrame.sample_rate, drspec.get_sample_rate(fh), 1)
        
        # Filter code
        self.assertEqual(cFrame.filter_code, 7)
        self.assertEqual(cFrame.filter_code, drspec.get_sample_rate(fh, FilterCode=True))
        
        # FFT windows per integration
        self.assertEqual(cFrame.ffts_per_integration, 6144)
        self.assertEqual(cFrame.ffts_per_integration, drspec.get_ffts_per_integration(fh))
        
        # Transform size
        self.assertEqual(cFrame.transform_size, 1024)
        self.assertEqual(cFrame.transform_size, drspec.get_transform_size(fh))
        
        # Integration time
        self.assertAlmostEqual(cFrame.integration_time, 0.32099265, 8)
        self.assertAlmostEqual(cFrame.integration_time, drspec.get_integration_time(fh), 8)
        
        fh.close()
        
    def test_drspec_comps(self):
        """Test the DR spectrometer frame comparison operators (>, <, etc.) for time tags."""

        fh = open(drspecFile, 'rb')
        # Frames 1 through 7
        frames = []
        for i in range(1,8):
            frames.append(drspec.read_frame(fh))
        fh.close()

        self.assertTrue(0 < frames[0])
        self.assertFalse(0 > frames[0])
        self.assertTrue(frames[-1] >= frames[0])
        self.assertFalse(frames[-1] <= frames[0])
        self.assertTrue(frames[0] == frames[0])
        self.assertFalse(frames[0] == frames[-1])
        self.assertFalse(frames[0] != frames[0])
        
    def test_drspec_sort(self):
        """Test sorting DR spectrometer frames by time tags."""
        
        fh = open(drspecFile, 'rb')
        # Frames 1 through 7
        frames = []
        for i in range(1,8):
            frames.append(drspec.read_frame(fh))
        
        frames.sort()
        frames = frames[::-1]
        
        for i in xrange(1,len(frames)):
            self.assertTrue( frames[i-1] >= frames[i] )
        fh.close()
        
    def test_drspec_math(self):
        """Test mathematical operations on DR spectrometer frame data via frames."""
        
        fh = open(drspecFile, 'rb')
        # Frames 1 through 7
        frames = []
        for i in range(1,8):
            frames.append(drspec.read_frame(fh))
        fh.close()
        
        npts = frames[0].data.XX0.size
        
        # Multiplication
        frameT = frames[0] * 2.0
        for i in xrange(npts):
            self.assertAlmostEqual(frameT.data.XX0[i], 2*frames[0].data.XX0[i], 2)
        frameT *= 2.0
        for i in xrange(npts):
            self.assertAlmostEqual(frameT.data.XX1[i], 4*frames[0].data.XX1[i], 2)
        frameT = frames[0] * frames[1]
        for i in xrange(npts):
            self.assertAlmostEqual(frameT.data.YY0[i], frames[0].data.YY0[i]*frames[1].data.YY0[i], 2)
            
        # Addition
        frameA = frames[0] + 2.0
        for i in xrange(npts):
            self.assertAlmostEqual(frameA.data.XX0[i], 2+frames[0].data.XX0[i], 2)
        frameA += 2.0
        for i in xrange(npts):
            self.assertAlmostEqual(frameA.data.XX1[i], 4+frames[0].data.XX1[i], 2)
        frameA = frames[0] + frames[1]
        for i in xrange(npts):
            self.assertAlmostEqual(frameA.data.YY0[i], frames[0].data.YY0[i]+frames[1].data.YY0[i], 2)
            
    ### VDIF ###
    
    def test_vdif_read(self):
        """Test reading in a frame from a VDIF file."""
        
        fh = open(vdifFile, 'rb')
        # First frame
        frame1 = vdif.read_frame(fh)
        
        # Validate header
        station, thread = frame1.id
        self.assertEqual(station, 19284)
        self.assertEqual(thread, 0)
        self.assertEqual(frame1.header.is_legacy, 0)
        self.assertEqual(frame1.header.is_invalid, 0)
        self.assertEqual(frame1.header.ref_epoch, 30)
        self.assertEqual(frame1.header.seconds_from_epoch, 9106862)
        self.assertEqual(frame1.header.frame_in_second, 59866)
        self.assertEqual(frame1.header.frame_length, 8224//8)
        self.assertEqual(frame1.header.nchan, 4)
        self.assertEqual(frame1.header.bits_per_sample, 2)
        self.assertEqual(frame1.header.is_complex, 0)
        
        # Validate (some) data
        for k,d in enumerate((1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ,-1.0)):
            i = k % frame1.header.nchan
            j = k // frame1.header.nchan
            self.assertAlmostEqual(frame1.data.data[i,j], d, 5)
            
        # Second frame
        frame2 = vdif.read_frame(fh)
        
        # Validate header
        station, thread = frame2.id
        self.assertEqual(station, 19284)
        self.assertEqual(thread, 0)
        self.assertEqual(frame2.header.is_legacy, 0)
        self.assertEqual(frame2.header.is_invalid, 0)
        self.assertEqual(frame2.header.ref_epoch, 30)
        self.assertEqual(frame2.header.seconds_from_epoch, 9106862)
        self.assertEqual(frame2.header.frame_in_second, 59867)
        self.assertEqual(frame2.header.frame_length, 8224//8)
        self.assertEqual(frame2.header.nchan, 4)
        self.assertEqual(frame2.header.bits_per_sample, 2)
        self.assertEqual(frame2.header.is_complex, 0)
        
        # Validate (some) data
        for k,d in enumerate((-1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ,1.0)):
            i = k % frame2.header.nchan
            j = k // frame2.header.nchan
            self.assertAlmostEqual(frame2.data.data[i,j], d, 5)
            
        fh.close()
        
    def test_vdif_errors(self):
        """Test reading in all frames from a truncated VDIF file."""
        
        fh = open(vdifFile, 'rb')
        # Frames 1 through 10
        for i in range(1,11):
            frame = vdif.read_frame(fh)
            
        # Last frame should be an error (errors.EOFError)
        self.assertRaises(errors.EOFError, vdif.read_frame, fh)
        fh.close()
        
    def test_vdif_threads(self):
        """Test finding out how many threads file."""
        
        fh = open(vdifFile, 'rb')
        nt = vdif.get_thread_count(fh)
        self.assertEqual(nt, 1)
        fh.close()
        
    def test_vdif_math(self):
        """Test mathematical operations on VDIF frame data via frames."""
        
        fh = open(vdifFile, 'rb')
        # Frames 1 through 10
        frames = []
        for i in range(1,11):
            frames.append(vdif.read_frame(fh))
        fh.close()
        
        nchan, nSamples = frames[0].data.data.shape
        
        # Multiplication
        frameT = frames[0] * 2.0
        for i in range(nchan):
            for j in range(nSamples):
                self.assertAlmostEqual(frameT.data.data[i,j], 2*frames[0].data.data[i,j], 2)
        frameT *= 2.0
        for i in range(nchan):
            for j in range(nSamples):
                self.assertAlmostEqual(frameT.data.data[i,j], 4*frames[0].data.data[i,j], 2)
        frameT = frames[0] * frames[1]
        for i in range(nchan):
            for j in range(nSamples):
                self.assertAlmostEqual(frameT.data.data[i,j], frames[0].data.data[i,j]*frames[1].data.data[i,j], 2)
            
        # Addition
        frameA = frames[0] + 2.0
        for i in range(nchan):
            for j in range(nSamples):
                self.assertAlmostEqual(frameA.data.data[i,j], 2+frames[0].data.data[i,j], 2)
        frameA += 2.0
        for i in range(nchan):
            for j in range(nSamples):
                self.assertAlmostEqual(frameA.data.data[i,j], 4+frames[0].data.data[i,j], 2)
        frameA = frames[0] + frames[1]
        for i in range(nchan):
            for j in range(nSamples):
                self.assertAlmostEqual(frameA.data.data[i,j], frames[0].data.data[i,j]+frames[1].data.data[i,j], 2)


class reader_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.reader units 
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(reader_tests)) 


if __name__ == '__main__':
    unittest.main()

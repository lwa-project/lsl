"""
Unit tests for the lsl.reader modules.
"""

import os
import numpy as np
import unittest
from datetime import timedelta

from lsl.reader import tbn
from lsl.reader import drx
from lsl.reader import vdif
from lsl.reader import drspec
from lsl.reader import errors
from lsl.reader.base import FrameTimestamp


__version__  = "0.9"
__author__    = "Jayce Dowell"


tbnFile = os.path.join(os.path.dirname(__file__), 'data', 'tbn-test.dat')
drxFile = os.path.join(os.path.dirname(__file__), 'data', 'drx-test.dat')
vdifFile = os.path.join(os.path.dirname(__file__), 'data', 'vdif-test.dat')
drspecFile = os.path.join(os.path.dirname(__file__), 'data', 'drspec-test.dat')


class reader_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.reader
    modules."""
    
    ### FrameTimestamp
    
    def test_timestamp(self):
        """Test creating a new FrameTimestamp"""
        
        t = FrameTimestamp(1587495778, 0.5)
        self.assertAlmostEqual(t.unix, 1587495778.5, 6)
        # https://www.epochconverter.com/
        dt = t.datetime
        self.assertEqual(dt.year, 2020)
        self.assertEqual(dt.month, 4)
        self.assertEqual(dt.day, 21)
        self.assertEqual(dt.hour, 19)
        self.assertEqual(dt.minute, 2)
        self.assertEqual(dt.second, 58)
        self.assertEqual(dt.microsecond, 500000)
        
        t = FrameTimestamp(1587495778.5)
        self.assertAlmostEqual(t.unix, 1587495778.5, 6)
        # https://www.epochconverter.com/
        dt = t.datetime
        self.assertEqual(dt.year, 2020)
        self.assertEqual(dt.month, 4)
        self.assertEqual(dt.day, 21)
        self.assertEqual(dt.hour, 19)
        self.assertEqual(dt.minute, 2)
        self.assertEqual(dt.second, 58)
        self.assertEqual(dt.microsecond, 500000)
        
        t = FrameTimestamp(1587495777.4, 1.1)
        self.assertAlmostEqual(t.unix, 1587495778.5, 6)
        # https://www.epochconverter.com/
        dt = t.datetime
        self.assertEqual(dt.year, 2020)
        self.assertEqual(dt.month, 4)
        self.assertEqual(dt.day, 21)
        self.assertEqual(dt.hour, 19)
        self.assertEqual(dt.minute, 2)
        self.assertEqual(dt.second, 58)
        self.assertEqual(dt.microsecond, 500000)
        
        t = FrameTimestamp.from_dp_timetag(1587495778*196000000 + 196000000//2)
        self.assertAlmostEqual(t.unix, 1587495778.5, 6)
        # https://planetcalc.com/503/
        self.assertAlmostEqual(t.mjd, 58960.79372685185+0.5/86400.0, 9)
        self.assertEqual(t.pulsar_mjd[0], 58960)
        self.assertEqual(t.pulsar_mjd[1], 68578)
        self.assertAlmostEqual(t.pulsar_mjd[2], 0.5, 9)
        self.assertEqual(t.dp_timetag, 1587495778*196000000 + 196000000//2)
        
        t = FrameTimestamp.from_mjd_mpm(58962, 60481519)
        # 200423 16:48:01  58962  60481519 T   1099467 1 SHL RPT POWER-OUTAGE|
        dt = t.datetime
        self.assertEqual(dt.year, 2020)
        self.assertEqual(dt.month, 4)
        self.assertEqual(dt.day, 23)
        self.assertEqual(dt.hour, 16)
        self.assertEqual(dt.minute, 48)
        self.assertEqual(dt.second, 1)
        self.assertEqual(dt.microsecond, 519000)
        self.assertTrue(dt.tzinfo is None)
        
        t = FrameTimestamp.from_mjd_mpm(58962, 60481519)
        # 200423 16:48:01  58962  60481519 T   1099467 1 SHL RPT POWER-OUTAGE|
        dt = t.utc_datetime
        self.assertEqual(dt.year, 2020)
        self.assertEqual(dt.month, 4)
        self.assertEqual(dt.day, 23)
        self.assertEqual(dt.hour, 16)
        self.assertEqual(dt.minute, 48)
        self.assertEqual(dt.second, 1)
        self.assertEqual(dt.microsecond, 519000)
        self.assertFalse(dt.tzinfo is None)
        
    def test_timestamp_string(self):
        """Test string representations of a FrameTimestamp"""
        
        t = FrameTimestamp.from_mjd_mpm(58962, 60481519)
        str(t)
        repr(t)
        
    def test_timestamp_add(self):
        """Test adding to a FrameTimestamp"""
        
        t = FrameTimestamp(1587495778, 0.5)
        t = t + 0.1
        self.assertAlmostEqual(t, FrameTimestamp(1587495778, 0.6), 10)
        self.assertAlmostEqual(t, 1587495778.6, 6)
        
        t += 1
        self.assertAlmostEqual(t, FrameTimestamp(1587495779, 0.6), 10)
        self.assertAlmostEqual(t, 1587495779.6, 6)
        
        t = t + timedelta(seconds=1)
        self.assertAlmostEqual(t, FrameTimestamp(1587495780, 0.6), 10)
        self.assertAlmostEqual(t, 1587495780.6, 6)
        
        t += timedelta(seconds=1, microseconds=400000)
        self.assertAlmostEqual(t, FrameTimestamp(1587495782, 0.0), 10)
        self.assertAlmostEqual(t, 1587495782.0, 6)
        
    def test_timestamp_sub(self):
        """Test subtracting from a FrameTimestamp"""
        
        t0 = FrameTimestamp(1587495778, 0.5)
        t1 = FrameTimestamp(1587495778, 0.1)
        self.assertAlmostEqual(t0-t1, 0.4, 9)
        
        t1 = FrameTimestamp(1587495778, 0.7)
        self.assertAlmostEqual(t0-t1, -0.2, 9)
        
        t1 = FrameTimestamp(1587495700, 0.7)
        self.assertAlmostEqual(t0-t1, 77.8, 9)
        
        t0 = t0 - 0.1
        self.assertAlmostEqual(t0, FrameTimestamp(1587495778, 0.4), 10)
        self.assertAlmostEqual(t0, 1587495778.4, 6)
        
        t0 -= 0.4
        self.assertAlmostEqual(t0, FrameTimestamp(1587495778, 0.0), 10)
        self.assertAlmostEqual(t0, 1587495778.0, 6)
        
        t0 = t0 - timedelta(seconds=1)
        self.assertAlmostEqual(t0, FrameTimestamp(1587495777, 0.0), 10)
        self.assertAlmostEqual(t0, 1587495777.0, 6)
        
        t0 -= timedelta(seconds=1, microseconds=500000)
        self.assertAlmostEqual(t0, FrameTimestamp(1587495775, 0.5), 10)
        self.assertAlmostEqual(t0, 1587495775.5, 6)
        
    def test_timestmp_cmp(self):
        """Test FrameTimestamp comparisons"""
        
        t0 = FrameTimestamp(1587495778, 0.5)
        t1 = FrameTimestamp(1587495779, 0.5)
        self.assertTrue(t0 <= t1)
        
        self.assertTrue(t0 > 0)
        self.assertTrue(t0 > 1587495778)
        self.assertTrue(t0 <= 1587495778.5)
        
        t0 = FrameTimestamp(1587495778, 0.0)
        self.assertEqual(t0, 1587495778)
        
    def test_timestamp_conversion(self):
        """Test FrameTimestamp string and numeric conversions"""
        
        t0 = FrameTimestamp(1587495778, 0.5)
        v = f"{t0:s}"
        self.assertEqual(v, str(t0))
        
        v = f"{t0:f}"
        self.assertEqual(float(v), t0.unix)
        
        v = f"{t0:d}"
        self.assertEqual(int(v), int(t0.unix))
        
    ### TBN ###
    
    def test_tbn_read(self):
        """Test reading in a frame from a TBN file."""
        
        fh = open(tbnFile, 'rb')
        # First frame is really TBN and stores the correct IDs
        frame1 = tbn.read_frame(fh)
        stand, pol = frame1.id
        self.assertEqual(stand, 1)
        self.assertEqual(pol, 0)
        str(frame1)
        repr(frame1)
        # Second frame
        frame2 = tbn.read_frame(fh)
        stand, pol = frame2.id
        self.assertEqual(stand, 1)
        self.assertEqual(pol, 1)
        str(frame2)
        repr(frame2)
        fh.close()
        
    def test_tbn_read_ci8(self):
        """Test reading in a frame from a TBN file, ci8 style."""
        
        fh = open(tbnFile, 'rb')
        frame1 = tbn.read_frame(fh)
        frame2 = tbn.read_frame(fh)
        fh.close()
        
        fh = open(tbnFile, 'rb')
        frame3 = tbn.read_frame_ci8(fh)
        frame4 = tbn.read_frame_ci8(fh)
        fh.close()
        
        # Compare
        data1 = frame3.payload.data['re'] + 1j*frame3.payload.data['im']
        data2 = frame4.payload.data['re'] + 1j*frame4.payload.data['im']
        for i in range(512):
            self.assertAlmostEqual(frame1.payload.data[i], data1[i], 1e-6)
            self.assertAlmostEqual(frame2.payload.data[i], data2[i], 1e-6)
            
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
        code = tbn.get_sample_rate(fh, filter_code=True)
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
        
        for i in range(1,len(frames)):
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
        np.testing.assert_allclose(frameT.payload.data, 2*frames[0].payload.data, atol=1e-6)
        frameT *= 2.0
        np.testing.assert_allclose(frameT.payload.data, 4*frames[0].payload.data, atol=1e-6)
        frameT = frames[0] * frames[1]
        np.testing.assert_allclose(frameT.payload.data, frames[0].payload.data*frames[1].payload.data, atol=1e-6)
        
        # Addition
        frameA = frames[0] + 2.0
        np.testing.assert_allclose(frameA.payload.data, 2+frames[0].payload.data, atol=1e-6)
        frameA += 2.0
        np.testing.assert_allclose(frameA.payload.data, 4+frames[0].payload.data, atol=1e-6)
        frameA = frames[0] + frames[1]
        np.testing.assert_allclose(frameA.payload.data, frames[0].payload.data+frames[1].payload.data, atol=1e-6)
        
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
        
    def test_drx_read_ci8(self):
        """Test reading in a frame from a DRX file, ci8 style."""
        
        fh = open(drxFile, 'rb')
        frame1 = drx.read_frame(fh)
        frame2 = drx.read_frame(fh)
        fh.close()
        
        fh = open(drxFile, 'rb')
        frame3 = drx.read_frame_ci8(fh)
        frame4 = drx.read_frame_ci8(fh)
        fh.close()
        
        # Compare
        data1 = frame3.payload.data['re'] + 1j*frame3.payload.data['im']
        data2 = frame4.payload.data['re'] + 1j*frame4.payload.data['im']
        for i in range(512):
            self.assertAlmostEqual(frame1.payload.data[i], data1[i], 1e-6)
            self.assertAlmostEqual(frame2.payload.data[i], data2[i], 1e-6)
            
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
        self.assertEqual(cFrame.filter_code, drx.get_sample_rate(fh, filter_code=True))
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
        
        for i in range(1,len(frames)):
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
        np.testing.assert_allclose(frameT.payload.data, 2*frames[0].payload.data, atol=1e-6)
        frameT *= 2.0
        np.testing.assert_allclose(frameT.payload.data, 4*frames[0].payload.data, atol=1e-6)
        frameT = frames[0] * frames[1]
        np.testing.assert_allclose(frameT.payload.data, frames[0].payload.data*frames[1].payload.data, atol=1e-6)
        
        # Division
        frameT = frames[0] / 2.0
        np.testing.assert_allclose(frameT.payload.data, 0.5*frames[0].payload.data, atol=1e-6)
        frameT /= 2.0
        np.testing.assert_allclose(frameT.payload.data, 0.25*frames[0].payload.data, atol=1e-6)
        frameT = frames[0] / frames[1]
        np.testing.assert_allclose(frameT.payload.data, frames[0].payload.data/frames[1].payload.data, atol=1e-6)
        
        # Addition
        frameA = frames[0] + 2.0
        np.testing.assert_allclose(frameA.payload.data, 2+frames[0].payload.data, atol=1e-6)
        frameA += 2.0
        np.testing.assert_allclose(frameA.payload.data, 4+frames[0].payload.data, atol=1e-6)
        frameA = frames[0] + frames[1]
        np.testing.assert_allclose(frameA.payload.data, frames[0].payload.data+frames[1].payload.data, atol=1e-6)
        
        # Subtraction
        frameA = frames[0] - 2.0
        np.testing.assert_allclose(frameA.payload.data, -2+frames[0].payload.data, atol=1e-6)
        frameA -= 2.0
        np.testing.assert_allclose(frameA.payload.data, -4+frames[0].payload.data, atol=1e-6)
        frameA = frames[0] - frames[1]
        np.testing.assert_allclose(frameA.payload.data, frames[0].payload.data-frames[1].payload.data, atol=1e-6)
            
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
        
        # The special "data" attribute
        data = frame2.payload.data
        self.assertEqual(len(data.shape), 1)
        self.assertEqual(data.shape[0], 1024)
        for attr in ('XX0', 'XX1', 'YY0', 'YY1'):
            d0 = getattr(frame2.payload, attr, None)
            d1 = data[attr]
            for i in range(1024):
                self.assertAlmostEqual(d0[i], d1[i], 6)
                
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
        self.assertEqual(cFrame.filter_code, drspec.get_sample_rate(fh, filter_code=True))
        
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
        
        for i in range(1,len(frames)):
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
        
        npts = frames[0].payload.XX0.size
        
        # Multiplication
        frameT = frames[0] * 2.0
        np.testing.assert_allclose(frameT.payload.XX0, 2*frames[0].payload.XX0, atol=1e-6)
        frameT *= 2.0
        np.testing.assert_allclose(frameT.payload.XX1, 4*frames[0].payload.XX1, atol=1e-6)
        frameT = frames[0] * frames[1]
        np.testing.assert_allclose(frameT.payload.YY0, frames[0].payload.YY0*frames[1].payload.YY0, atol=1e-6)
        
        # Division
        frameT = frames[0] / 2.0
        np.testing.assert_allclose(frameT.payload.XX0, 0.5*frames[0].payload.XX0, atol=1e-6)
        frameT /= 2.0
        np.testing.assert_allclose(frameT.payload.XX1, 0.25*frames[0].payload.XX1, atol=1e-6)
        frameT = frames[0] / frames[1]
        np.testing.assert_allclose(frameT.payload.YY0, frames[0].payload.YY0/frames[1].payload.YY0, atol=1e-6)
        
        # Addition
        frameA = frames[0] + 2.0
        np.testing.assert_allclose(frameA.payload.XX0, 2+frames[0].payload.XX0, atol=1e-6)
        frameA += 2.0
        np.testing.assert_allclose(frameA.payload.XX0, 4+frames[0].payload.XX0, atol=1e-6)
        frameA = frames[0] + frames[1]
        np.testing.assert_allclose(frameA.payload.YY0, frames[0].payload.YY0+frames[1].payload.YY0, atol=1e-6)
        
        # Subtraction
        frameA = frames[0] - 2.0
        np.testing.assert_allclose(frameA.payload.XX0, -2+frames[0].payload.XX0, atol=1e-6)
        frameA -= 2.0
        np.testing.assert_allclose(frameA.payload.XX0, -4+frames[0].payload.XX0, atol=1e-6)
        frameA = frames[0] - frames[1]
        np.testing.assert_allclose(frameA.payload.YY0, frames[0].payload.YY0-frames[1].payload.YY0, atol=1e-6)
        
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
            self.assertAlmostEqual(frame1.payload.data[i,j], d, 5)
            
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
            self.assertAlmostEqual(frame2.payload.data[i,j], d, 5)
            
        fh.close()
        
    def test_vdif_read_i8(self):
        """Test reading in a frame from a VDIF file, i8 style."""
        
        fh = open(vdifFile, 'rb')
        frame1 = vdif.read_frame(fh)
        frame2 = vdif.read_frame(fh)
        fh.close()
        
        fh = open(vdifFile, 'rb')
        frame3 = vdif.read_frame_i8(fh)
        frame4 = vdif.read_frame_i8(fh)
        fh.close()
        
        # Validate the data
        for i in range(frame1.payload.data.shape[0]):
            for j in range(frame1.payload.data.shape[1]):
                self.assertAlmostEqual(round(frame1.payload.data[i,j]), frame3.payload.data[i,j], 5)
                self.assertAlmostEqual(round(frame2.payload.data[i,j]), frame4.payload.data[i,j], 5)
                
    def test_vdif_time(self):
        """Test VDIF time interpretation."""
        
        f = vdif.FrameHeader(seconds_from_epoch=623633800, ref_epoch=0,
                             frame_in_second=0, version=0)
        self.assertEqual(f.time.unix, 1570318600)
        
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
        
    def test_vdif_comps(self):
        """Test the VDIF frame comparison operators (>, <, etc.) for time tags."""

        fh = open(vdifFile, 'rb')
        # Frames 1 through 10
        frames = []
        for i in range(1,11):
            frames.append(vdif.read_frame(fh, sample_rate=16e6))
        fh.close()

        self.assertTrue(0 < frames[0])
        self.assertFalse(0 > frames[0])
        self.assertTrue(frames[-1] >= frames[0])
        self.assertFalse(frames[-1] <= frames[0])
        self.assertTrue(frames[0] == frames[0])
        self.assertFalse(frames[0] == frames[-1])
        self.assertFalse(frames[0] != frames[0])
        
    def test_vdif_sort(self):
        """Test sorting VDIF frames by time tags."""
        
        fh = open(vdifFile, 'rb')
        # Frames 1 through 10
        frames = []
        for i in range(1,11):
            frames.append(vdif.read_frame(fh, sample_rate=16e6))
        
        frames.sort()
        frames = frames[::-1]
        
        for i in range(1,len(frames)):
            self.assertTrue( frames[i-1] >= frames[i] )
        fh.close()
        
    def test_vdif_math(self):
        """Test mathematical operations on VDIF frame data via frames."""
        
        fh = open(vdifFile, 'rb')
        # Frames 1 through 10
        frames = []
        for i in range(1,11):
            frames.append(vdif.read_frame(fh))
        fh.close()
        
        nchan, nSamples = frames[0].payload.data.shape
        
        # Multiplication
        frameT = frames[0] * 2.0
        np.testing.assert_allclose(frameT.payload.data, 2*frames[0].payload.data, atol=1e-6)
        frameT *= 2.0
        np.testing.assert_allclose(frameT.payload.data, 4*frames[0].payload.data, atol=1e-6)
        frameT = frames[0] * frames[1]
        np.testing.assert_allclose(frameT.payload.data, frames[0].payload.data*frames[1].payload.data, atol=1e-6)
        
        # Addition
        frameA = frames[0] + 2.0
        np.testing.assert_allclose(frameA.payload.data, 2+frames[0].payload.data, atol=1e-6)
        frameA += 2.0
        np.testing.assert_allclose(frameA.payload.data, 4+frames[0].payload.data, atol=1e-6)
        frameA = frames[0] + frames[1]
        np.testing.assert_allclose(frameA.payload.data, frames[0].payload.data+frames[1].payload.data, atol=1e-6)


class reader_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.reader units 
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(reader_tests)) 


if __name__ == '__main__':
    unittest.main()

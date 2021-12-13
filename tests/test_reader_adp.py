"""
Unit tests for the lsl.reader modules that are for ADP stations.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import unittest

from lsl.common.paths import DATA_BUILD
from lsl.reader import tbf
from lsl.reader import cor
from lsl.reader import errors


__version__  = "0.2"
__author__    = "Jayce Dowell"


tbfFile = os.path.join(DATA_BUILD, 'tests', 'tbf-test.dat')
corFile = os.path.join(DATA_BUILD, 'tests', 'cor-test.dat')


class reader_adp_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.reader
    modules."""
    
    ### TBF ###
    
    def test_tbf_read(self):
        """Test reading in a frame from a TBF file."""
        
        fh = open(tbfFile, 'rb')
        # First frame is really TBF and stores the first channel
        frame1 = tbf.read_frame(fh)
        self.assertTrue(frame1.header.is_tbf)
        self.assertEqual(frame1.header.first_chan, 2348)
        # Second frame
        frame2 = tbf.read_frame(fh)
        self.assertTrue(frame2.header.is_tbf)
        self.assertEqual(frame2.header.first_chan, 2360)
        fh.close()
        
    def test_tbf_read_ci8(self):
        """Test reading in a frame from a TBF file, ci8 style."""
        
        fh = open(tbfFile, 'rb')
        frame1 = tbf.read_frame(fh)
        frame2 = tbf.read_frame(fh)
        fh.close()
        
        # Compare
        data1 = frame1.payload.data_ci8[:,:,:,0] + 1j*frame1.payload.data_ci8[:,:,:,1]
        data2 = frame2.payload.data_ci8[:,:,:,0] + 1j*frame2.payload.data_ci8[:,:,:,1]
        for i in range(800):
            c = i // 2 // 256
            s = i // 2 % 256
            p = i % 2
            self.assertAlmostEqual(frame1.payload.data[c,s,p], data1[c,s,p], 1e-6)
            self.assertAlmostEqual(frame2.payload.data[c,s,p], data2[c,s,p], 1e-6)
            
    def test_tbf_errors(self):
        """Test TBF reading errors."""
        
        fh = open(tbfFile, 'rb')
        # Frames 1 through 5
        for i in range(1,6):
            frame = tbf.read_frame(fh)
            
        # Last frame should be an error (errors.EOFError)
        self.assertRaises(errors.EOFError, tbf.read_frame, fh)
        fh.close()
        
        # If we offset in the file by 1 byte, we should be a 
        # sync error (errors.SyncError).
        fh = open(tbfFile, 'rb')
        fh.seek(1)
        self.assertRaises(errors.SyncError, tbf.read_frame, fh)
        fh.close()
        
    def test_tbf_comps(self):
        """Test the TBF frame comparison operators (>, <, etc.) for time tags."""
        
        fh = open(tbfFile, 'rb')
        # Frames 1 through 4
        frames = []
        for i in range(1,5):
            frames.append(tbf.read_frame(fh))
        fh.close()
        
        self.assertTrue(0 < frames[0])
        self.assertFalse(0 > frames[0])
        self.assertTrue(frames[-1] >= frames[0])
        self.assertTrue(frames[-1] <= frames[0])
        self.assertTrue(frames[0] == frames[0])
        self.assertTrue(frames[0] == frames[-1])
        self.assertFalse(frames[0] != frames[0])
        
    def test_tbf_sort(self):
        """Test sorting TBF frames by time tags."""
        
        fh = open(tbfFile, 'rb')
        # Frames 1 through 3
        frames = []
        for i in range(1,4):
            frames.append(tbf.read_frame(fh))
        fh.close()
        
        frames.sort()
        frames = frames[::-1]
        
        for i in range(1,len(frames)):
            self.assertTrue( frames[i-1] >= frames[i] )
            
    def test_tbf_math(self):
        """Test mathematical operations on TBF frame data via frames."""
        
        fh = open(tbfFile, 'rb')
        # Frames 1 through 3
        frames = []
        for i in range(1,4):
            frames.append(tbf.read_frame(fh))
        fh.close()
        
        # Multiplication
        frameT = frames[0] * 2.0
        for i in range(12*256*2):
            c = i // 2 // 256
            s = i // 2 % 256
            p = i % 2
            self.assertAlmostEqual(frameT.payload.data[c,s,p], 2*frames[0].payload.data[c,s,p], 2)
        frameT *= 2.0
        for i in range(12*256*2):
            c = i // 2 // 256
            s = i // 2 % 256
            p = i % 2
            self.assertAlmostEqual(frameT.payload.data[c,s,p], 4*frames[0].payload.data[c,s,p], 2)
        frameT = frames[0] * frames[1]
        for i in range(12*256*2):
            c = i // 2 // 256
            s = i // 2 % 256
            p = i % 2
            self.assertAlmostEqual(frameT.payload.data[c,s,p], frames[0].payload.data[c,s,p]*frames[1].payload.data[c,s,p], 2)
            
        # Addition
        frameA = frames[0] + 2.0
        for i in range(800):
            c = i // 2 // 256
            s = i // 2 % 256
            p = i % 2
            self.assertAlmostEqual(frameA.payload.data[c,s,p], 2+frames[0].payload.data[c,s,p], 2)
        frameA += 2.0
        for i in range(800):
            c = i // 2 // 256
            s = i // 2 % 256
            p = i % 2
            self.assertAlmostEqual(frameA.payload.data[c,s,p], 4+frames[0].payload.data[c,s,p], 2)
        frameA = frames[0] + frames[1]
        for i in range(800):
            c = i // 2 // 256
            s = i // 2 % 256
            p = i % 2
            self.assertAlmostEqual(frameA.payload.data[c,s,p], frames[0].payload.data[c,s,p]+frames[1].payload.data[c,s,p], 2)
            
     ### COR ###
    
    def test_cor_read(self):
        """Test reading in a frame from a COR file."""
        
        fh = open(corFile, 'rb')
        # First frame is really COR and check the basic metadata
        frame1 = cor.read_frame(fh)
        self.assertTrue(frame1.header.is_cor)
        self.assertEqual(frame1.header.first_chan, 1584)
        self.assertEqual(frame1.id, (1,1))
        self.assertEqual(frame1.integration_time, 5)
        self.assertEqual(frame1.gain, 1)
        # Second frame
        frame2 = cor.read_frame(fh)
        self.assertTrue(frame2.header.is_cor)
        self.assertEqual(frame2.header.first_chan, 1584)
        self.assertEqual(frame2.id, (1,2))
        self.assertEqual(frame2.integration_time, 5)
        self.assertEqual(frame2.gain, 1)
        fh.close()
        
    def test_cor_errors(self):
        """Test COR reading errors."""
        
        fh = open(corFile, 'rb')
        # Frames 1 through 65
        for i in range(1,66):
            frame = cor.read_frame(fh)
            
        # Last frame should be an error (errors.EOFError)
        self.assertRaises(errors.EOFError, cor.read_frame, fh)
        fh.close()
        
        # If we offset in the file by 1 byte, we should be a 
        # sync error (errors.SyncError).
        fh = open(corFile, 'rb')
        fh.seek(1)
        self.assertRaises(errors.SyncError, cor.read_frame, fh)
        fh.close()
        
    def test_cor_frames(self):
        """Test determing the number of frames per observation in a COR file."""
        
        fh = open(corFile, 'rb')
        self.assertEqual(32896, cor.get_frames_per_obs(fh))
        fh.close()
        
    def test_cor_channels(self):
        """Test determing the number of channels in a COR observation."""
        
        fh = open(corFile, 'rb')
        self.assertEqual(72, cor.get_channel_count(fh))
        fh.close()
        
    def test_cor_baselines(self):
        """Test determing the number of baselines in a COR observation."""
        
        fh = open(corFile, 'rb')
        self.assertEqual(32896, cor.get_baseline_count(fh))
        fh.close()
        
    def test_cor_comps(self):
        """Test the COR frame comparison operators (>, <, etc.) for time tags."""
        
        fh = open(corFile, 'rb')
        # Frames 1 through 29
        frames = []
        for i in range(1,30):
            frames.append(cor.read_frame(fh))
        fh.close()
        
        self.assertTrue(0 < frames[0])
        self.assertFalse(0 > frames[0])
        self.assertTrue(frames[-1] >= frames[0])
        self.assertTrue(frames[-1] <= frames[0])
        self.assertTrue(frames[0] == frames[0])
        self.assertFalse(frames[0] != frames[-1])
        self.assertFalse(frames[0] != frames[0])
        
    def test_cor_sort(self):
        """Test sorting COR frames by time tags."""
        
        fh = open(corFile, 'rb')
        # Frames 1 through 29
        frames = []
        for i in range(1,30):
            frames.append(cor.read_frame(fh))
        fh.close()
        
        frames.sort()
        frames = frames[::-1]
        
        for i in range(1,len(frames)):
            self.assertTrue( frames[i-1] >= frames[i] )
            
    def test_cor_math(self):
        """Test mathematical operations on COR frame data via frames."""
        
        fh = open(corFile, 'rb')
        # Frames 1 through 29
        frames = []
        for i in range(1,30):
            frames.append(cor.read_frame(fh))
        fh.close()
        
        # Multiplication
        frameT = frames[0] * 2.0
        for i in range(72):
            for j in range(2):
                for k in range(2):
                    self.assertAlmostEqual(frameT.payload.data[i,j,k], 2*frames[0].payload.data[i,j,k], 2)
        frameT *= 2.0
        for i in range(72):
            for j in range(2):
                for k in range(2):
                    self.assertAlmostEqual(frameT.payload.data[i,j,k], 4*frames[0].payload.data[i,j,k], 2)
        frameT = frames[0] * frames[1]
        for i in range(72):
            for j in range(2):
                for k in range(2):
                    self.assertAlmostEqual(frameT.payload.data[i,j,k], frames[0].payload.data[i,j,k]*frames[1].payload.data[i,j,k], 2)
            
        # Addition
        frameA = frames[0] + 2.0
        for i in range(72):
            for j in range(2):
                for k in range(2):
                    self.assertAlmostEqual(frameA.payload.data[i,j,k], 2+frames[0].payload.data[i,j,k], 2)
        frameA += 2.0
        for i in range(72):
            for j in range(2):
                for k in range(2):
                    self.assertAlmostEqual(frameA.payload.data[i,j,k], 4+frames[0].payload.data[i,j,k], 2)
        frameA = frames[0] + frames[1]
        for i in range(72):
            for j in range(2):
                for k in range(2):
                    self.assertAlmostEqual(frameA.payload.data[i,j,k], frames[0].payload.data[i,j,k]+frames[1].payload.data[i,j,k], 2)
            
        

class reader_adp_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.reader units 
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(reader_adp_tests)) 


if __name__ == '__main__':
    unittest.main()

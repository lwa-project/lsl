"""
Unit test for the lsl.misc.file_cache module.
"""

import os
import unittest
import tempfile
import shutil

from lsl.misc import file_cache


__version__  = "0.1"
__author__    = "Jayce Dowell"


class cache_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.misc.file_cache
    module."""
    
    def setUp(self):
        self.cache_dir = tempfile.mkdtemp(prefix='test-cache-', suffix='.tmp')
    
    def test_file_cache(self):
        c = file_cache.FileCache(self.cache_dir)
        self.assertFalse('test.txt' in c)
        
        with c.open('test.txt', 'w') as fh:
            fh.write('test contents')
        self.assertTrue('test.txt' in c)
        
        with c.open('test.txt', 'r') as fh:
            contents = fh.read()
        self.assertEqual(contents, 'test contents')
        
        with c.open('test.dat', 'wb') as fh:
            fh.write(b'test contents binary')
        self.assertTrue('test.dat' in c)
        
        with c.open('test.dat', 'rb') as fh:
            contents = fh.read()
        self.assertEqual(contents, b'test contents binary')
        
        c.getsize('test.txt')
        self.assertRaises(OSError, c.getsize, 'notthere.dat')
        
        c.stat('test.txt')
        self.assertRaises(OSError, c.stat, 'notthere.dat')
        
        c.remove('test.txt')
        self.assertFalse('test.txt' in c)
        
    def test_memory_cache(self):
        c = file_cache.MemoryCache()
        self.assertFalse('test.txt' in c)
        
        with c.open('test.txt', 'w') as fh:
            fh.write('test contents')
        self.assertTrue('test.txt' in c)
        
        with c.open('test.txt', 'r') as fh:
            contents = fh.read()
        self.assertEqual(contents, 'test contents')
        
        with c.open('test.dat', 'wb') as fh:
            fh.write(b'test contents binary')
        self.assertTrue('test.dat' in c)
        
        with c.open('test.dat', 'rb') as fh:
            contents = fh.read()
        self.assertEqual(contents, b'test contents binary')
        
        c.getsize('test.txt')
        self.assertRaises(OSError, c.getsize, 'notthere.dat')
        
        c.stat('test.txt')
        self.assertRaises(OSError, c.stat, 'notthere.dat')
        
        c.remove('test.txt')
        self.assertFalse('test.txt' in c)
        
    def test_max_size(self):
        with self.subTest(type='file'):
            c = file_cache.FileCache(self.cache_dir, max_size=16/1024**2)
            
            file_list = []
            for i in range(10):
                file_list.append(f"test{i+1}.txt")
                with c.open(file_list[-1], 'w') as fh:
                    fh.write(f"temp file #{i+1}")
                    
                for j in range(0, len(file_list)-2):
                    self.assertFalse(file_list[j] in c)
                for j in range(len(file_list)-2, len(file_list)):
                    self.assertTrue(file_list[j] in c)
                    
        with self.subTest(type='memory'):
            c = file_cache.MemoryCache(max_size=16/1024**2)
            
            file_list = []
            for i in range(10):
                file_list.append(f"test{i+1}.txt")
                with c.open(file_list[-1], 'w') as fh:
                    fh.write(f"temp file #{i+1}")
                    
                for j in range(0, len(file_list)-2):
                    self.assertFalse(file_list[j] in c)
                for j in range(len(file_list)-2, len(file_list)):
                    self.assertTrue(file_list[j] in c)
                    
            
    def tearDown(self):
        shutil.rmtree(self.cache_dir, ignore_errors=True)
            


class  cache_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.misc.file_cache units 
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(cache_tests)) 


if __name__ == '__main__':
    unittest.main()

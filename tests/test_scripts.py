# -*- coding: utf-8 -*-

"""Unit test for LSL scripts."""

import unittest
import glob
import sys
import re
import os

from lsl.common.paths import MODULE_BUILD

run_scripts_tests = False
try:
    from pylint import epylint as lint
    run_scripts_tests = True
except ImportError:
    pass


__version__  = "0.1"
__revision__ = "$Rev$"
__author__   = "Jayce Dowell"


_LINT_RE = re.compile('(?P<module>.*?)\:(?P<line>\d+)\: \[(?P<type>.*?)\] (?P<info>.*)')


@unittest.skipUnless(run_scripts_tests, "requires the 'pylint' module")
class scripts_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the LSL scripts."""
    
    pass     


def _test_generator(script):
    """
    Function to build a test method for each script that is provided.  
    Returns a function that is suitable as a method inside a unittest.TestCase
    class
    """
    
    def test(self):
        out, err = lint.py_run("%s -E --init-hook='import sys; sys.path=[%s]; sys.path.insert(0, \"%s\")'" % (script, ",".join(['"%s"' % p for p in sys.path]), os.path.dirname(MODULE_BUILD)), return_std=True)
        out_lines = out.read().split('\n')
        err_lines = err.read().split('\n')
        out.close()
        err.close()
        
        for line in out_lines:
            if line.find("Module 'numpy") != -1:
                continue
                
            mtch = _LINT_RE.match(line)
            if mtch is not None:
                line_no, type, info = mtch.group('line'), mtch.group('type'), mtch.group('info')
                self.assertEqual(type.find('undefined'), -1, "%s:%s - %s" % (os.path.basename(script), line_no, info))
                # HACK to deal with some strangeness in lwa_cat_view.py
                if script.find('lwa_cat_view.py') == -1 or (script.find('lwa_cat_view.py') != -1 and info.find('j2000_') == -1):
                    self.assertEqual(type.find('no-member'), -1, "%s:%s - %s" % (os.path.basename(script), line_no, info))
    return test


_SCRIPTS = glob.glob(os.path.join(MODULE_BUILD, '..', 'scripts', '*.py'))
_SCRIPTS.sort()
for script in _SCRIPTS:
    test = _test_generator(script)
    name = 'test_%s' % os.path.splitext(os.path.basename(script))[0]
    doc = """Static analysis of the '%s' script.""" % os.path.basename(script)
    setattr(test, '__doc__', doc)
    setattr(scripts_tests, name, test)


class scripts_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the LSL script
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(scripts_tests))


if __name__ == '__main__':
    unittest.main()
    

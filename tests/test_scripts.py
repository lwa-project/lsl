"""
Unit tests for the various LSL scripts.
"""

import unittest
import glob
import sys
import os
import json
import subprocess

from lsl.common.paths import MODULE_BUILD

run_scripts_tests = False
try:
    from io import StringIO
    from pylint.lint import Run
    from pylint.reporters.json_reporter import JSONReporter
    run_scripts_tests = True
except ImportError:
    pass


__version__  = "0.2"
__author__   = "Jayce Dowell"


@unittest.skipUnless(run_scripts_tests, "requires the 'pylint' module")
class scripts_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the LSL scripts."""
    
    def test_scripts(self):
        """Static analysis of the LSL scripts."""
        
        _SCRIPTS = glob.glob(os.path.join(MODULE_BUILD, '..', 'scripts', '*.py'))
        _SCRIPTS.sort()
        for script in _SCRIPTS:
            name = script.rsplit('scripts'+os.path.sep)[-1]
            with self.subTest(script=name):
                pylint_output = StringIO()
                reporter = JSONReporter(pylint_output)
                pylint_args = [script, "-E", "--extension-pkg-whitelist=numpy", "--init-hook='import sys; sys.path=[%s]; sys.path.insert(0, \"%s\")'" % (",".join(['"%s"' % p for p in sys.path]), os.path.dirname(MODULE_BUILD))]
                Run(pylint_args, reporter=reporter, exit=False)
                results = json.loads(pylint_output.getvalue())
                
                for i,entry in enumerate(results):
                    with self.subTest(error_number=i+1):
                        if entry['symbol'] == 'no-member' and entry['message'].startswith("Instance of 'HDUList'"):
                            continue
                            
                        self.assertTrue(False, f"{entry['path']}:{entry['line']} - {entry['message']}")
                        
    def test_help(self):
        """Help/documentation in the LSL scripts."""
        
        _SCRIPTS = glob.glob(os.path.join(MODULE_BUILD, '..', 'scripts', '*.py'))
        _SCRIPTS.sort()
        for script in _SCRIPTS:
            name = script.rsplit('scripts'+os.path.sep)[-1]
            with self.subTest(script=name):
                try:
                    status = subprocess.check_call([sys.executable, script, '--help'],
                                                   stdout=subprocess.DEVNULL,# stderr=subprocess.DEVNULL,
                                                   cwd=os.path.dirname(MODULE_BUILD))
                    self.assertTrue(status == 0, f"Non-zero exit code when running script with '--help' flag: {status}")
                except subprocess.CalledProcessError as e:
                    self.assertTrue(False, f"Failed to run script with '--help' flag: {str(e)}")


class scripts_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the LSL script
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(scripts_tests))


if __name__ == '__main__':
    unittest.main()

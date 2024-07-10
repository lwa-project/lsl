"""
Unit tests for the various LSL scripts.
"""

import unittest
import glob
import sys
import re
import os

from lsl.common.paths import MODULE_BUILD

run_scripts_tests = False
try:
    from io import StringIO
    from pylint.lint import Run
    from pylint.reporters.text import TextReporter
    run_scripts_tests = True
except ImportError:
    pass


__version__  = "0.1"
__author__   = "Jayce Dowell"


_LINT_RE = re.compile(r'(?P<module>.*?):(?P<line>\d+): \[(?P<type>.*?)\] (?P<info>.*)')


@unittest.skipUnless(run_scripts_tests, "requires the 'pylint' module")
class scripts_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the LSL scripts."""
    
    def test_scripts(self):
        """Static analysis of the LSL scripts."""
        
        _SCRIPTS = glob.glob(os.path.join(MODULE_BUILD, '..', 'scripts', '*.py'))
        _SCRIPTS.sort()
        for script in _SCRIPTS:
            name = script.rsplit('scripts')[-1]
            with self.subTest(script=name):
                pylint_output = StringIO()
                reporter = TextReporter(pylint_output)
                pylint_args = [script, "-E", "--extension-pkg-whitelist=numpy", "--init-hook='import sys; sys.path=[%s]; sys.path.insert(0, \"%s\")'" % (",".join(['"%s"' % p for p in sys.path]), os.path.dirname(MODULE_BUILD))]
                Run(pylint_args, reporter=reporter, exit=False)
                out = pylint_output.getvalue()
                out_lines = out.split('\n')
                
                for line in out_lines:
                    mtch = _LINT_RE.match(line)
                    if mtch is not None:
                        line_no, type, info = mtch.group('line'), mtch.group('type'), mtch.group('info')
                        self.assertEqual(type.find('syntax'), -1, "%s:%s - %s" % (os.path.basename(script), line_no, info))
                        self.assertEqual(type.find('undefined'), -1, "%s:%s - %s" % (os.path.basename(script), line_no, info))
                        # HACK to deal with some strangeness in lwa_cat_view.py
                        if script.find('lwa_cat_view.py') == -1 or (script.find('lwa_cat_view.py') != -1 and info.find('j2000_') == -1):
                            self.assertEqual(type.find('no-member'), -1, "%s:%s - %s" % (os.path.basename(script), line_no, info))


class scripts_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the LSL script
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(scripts_tests))


if __name__ == '__main__':
    unittest.main()

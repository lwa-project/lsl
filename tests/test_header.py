"""
Unit test for the lsl.common.header module.
"""

import argparse
import warnings
import unittest
import sys

from lsl.common import header
from lsl.version import version as LSL_VERSION


__version__ = "0.1"
__author__  = "Jayce Dowell"


class header_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.common.header
    module."""
    
    def test_width_too_small(self):
        """Test that a width that leaves too little room is rejected."""
        
        # 'left': ewidth = width - comment - indent
        self.assertRaises(ValueError, header.CommentHeader, width=5)
        # 'box': ewidth = width - 2*comment - indent - base_indent
        self.assertRaises(ValueError, header.BoxHeader, width=10)
        
    def test_effective_width(self):
        """Test the effective width calculation."""
        
        h = header.CommentHeader(width=40)
        self.assertEqual(h.effective_width, 40 - 1 - 1)
        
        h = header.BoxHeader(width=40)
        self.assertEqual(h.effective_width, 40 - 2*1 - 1 - 1)
        
    def test_narrow_width_warns(self):
        """Test that a narrow-but-usable width logs a warning."""
        
        with self.assertLogs(header.LSL_LOGGER, level='WARNING'):
            header.CommentHeader(width=17)    # ewidth = 15, in [8, 16)
            
    def test_comment_prefix(self):
        """Test that every rendered line starts with the comment string."""
        
        h = header.CommentHeader(width=0)
        with h.section('Parameters'):
            h.key_value('APARM', 'this')
            h.line("a free-form comment")
        for line in str(h).split('\n'):
            if line:
                self.assertTrue(line.startswith('#'))
                
    def test_section_indent(self):
        """Test that a section's contents are indented past its header."""
        
        h = header.CommentHeader(width=0)
        with h.section('Parameters'):
            h.key_value('APARM', 'this')
        lines = [l for l in str(h).split('\n') if l]
        head = lines[0]                   # '# Parameters:'
        item = lines[1]                   # '#  APARM: ...'
        head_indent = len(head) - len(head.lstrip('# '))
        item_indent = len(item) - len(item.lstrip('# '))
        self.assertGreater(item_indent, head_indent)
        
    def test_key_value_json(self):
        """Test that key_value renders values as JSON."""
        
        h = header.CommentHeader(width=0)
        h.key_value('S', 'text')
        h.key_value('I', 3)
        h.key_value('F', 0.5)
        h.key_value('B', True)
        h.key_value('N', None)
        
        text = str(h)
        self.assertIn('S: "text"', text)
        self.assertIn('I: 3', text)
        self.assertIn('F: 0.5', text)
        self.assertIn('B: true', text)
        self.assertIn('N: null', text)
        
    def test_blank_and_rule(self):
        """Test the shapes of blank() and rule() in the 'left' style."""
        
        h = header.CommentHeader(width=0)
        h.key_value('K', 'value')
        h.blank()
        h.rule()
        
        lines = str(h).split('\n')
        self.assertIn('#', lines)         # blank -> bare comment
        self.assertIn('####', lines)      # rule -> comment*4
        
    def test_clear(self):
        """Test that clear() empties the header."""
        
        h = header.CommentHeader(width=0)
        h.key_value('K', 'value')
        self.assertEqual(len(h._contents), 1)
        
        h.clear()
        self.assertEqual(len(h._contents), 0)
        
    def test_empty(self):
        """Test that an empty header renders as nothing."""
        
        self.assertEqual(str(header.CommentHeader(width=0)), '')
        self.assertEqual(str(header.BoxHeader(width=40)), '')
        
    def test_box_borders(self):
        """Test that the box style is fully enclosed and uniform in width."""
        
        h = header.BoxHeader(width=0)
        with h.section('Parameters'):
            h.key_value('APARM', 'this')
        lines = [l for l in str(h).split('\n') if l]
        
        # Top and bottom are solid rules
        self.assertEqual(set(lines[0]), set('#'))
        self.assertEqual(set(lines[-1]), set('#'))
        # All lines share a single width and have side walls
        width = len(lines[0])
        for line in lines:
            self.assertEqual(len(line), width)
            self.assertTrue(line.startswith('#'))
            self.assertTrue(line.endswith('#'))
            
    def test_roundtrip_types(self):
        """Test that each JSON value type survives a round-trip."""
        
        h = header.CommentHeader(width=0)
        with h.section('Parameters'):
            h.key_value('S', 'text')
            h.key_value('I', 3)
            h.key_value('F', 0.5)
            h.key_value('B', True)
            h.key_value('N', None)
            
        d = header.parse_header(str(h))
        self.assertEqual(d, {'Parameters': {'S': 'text', 'I': 3, 'F': 0.5,
                                            'B': True, 'N': None}})
        
    def test_roundtrip_nested(self):
        """Test that nested sections round-trip into nested dictionaries."""
        
        h = header.CommentHeader(width=0)
        with h.section('Comments'):
            with h.section('General'):
                h.key_value('A', 1)
            with h.section('Specific'):
                h.key_value('B', 2)
        h.key_value('Top', 3)
        
        d = header.parse_header(str(h))
        self.assertEqual(d, {'Comments': {'General': {'A': 1},
                                          'Specific': {'B': 2}},
                             'Top': 3})
        
    def test_roundtrip_box(self):
        """Test that the box style round-trips the same as the left style."""
        
        h = header.BoxHeader(width=0)
        with h.section('Parameters'):
            h.key_value('A', 'x')
            h.key_value('N', 3)
            
        d = header.parse_header(str(h))
        self.assertEqual(d, {'Parameters': {'A': 'x', 'N': 3}})
        
    def test_parse_accepts_list(self):
        """Test that parse_header accepts a list of lines as well as a string."""
        
        h = header.CommentHeader(width=0)
        h.key_value('K', 'value')
        
        text = str(h)
        self.assertEqual(header.parse_header(text),
                         header.parse_header(text.split('\n')))
        
    def test_roundtrip_value_with_colon(self):
        """Test that a value containing a colon survives a round-trip."""
        
        h = header.CommentHeader(width=0)
        h.key_value('When', '2026-06-08T21:42:53+00:00')
        
        d = header.parse_header(str(h))
        self.assertEqual(d['When'], '2026-06-08T21:42:53+00:00')
        
    def test_line_into_freeform(self):
        """Test that free-form line() content is parsed under '_freeform'."""
        
        h = header.CommentHeader(width=0)
        with h.section('Notes'):
            h.line("just some prose")
            
        d = header.parse_header(str(h))
        self.assertIn('_freeform', d['Notes'])
        
    def test_timestamp(self):
        """Test that timestamp() adds a parseable key."""
        
        h = header.CommentHeader(width=0)
        h.timestamp('Updated')
        d = header.parse_header(str(h))
        self.assertIn('Updated', d)
        self.assertIsInstance(d['Updated'], str)
        
    def test_lsl_version(self):
        """Test that lsl_version() records the current LSL version."""
        
        h = header.CommentHeader(width=0)
        h.lsl_version()
        d = header.parse_header(str(h))
        self.assertEqual(d['LSL Version'], LSL_VERSION)
        
    def test_fixed_width_wraps(self):
        """Test that a long value wraps under a fixed width."""
        
        long_value = 'word ' * 8
        
        # A pinned, narrow width forces wrapping onto extra lines
        h = header.CommentHeader(width=20)
        with h.section('P'):
            h.key_value('K', long_value)
            h.key_value('I', sys.maxsize)
            h.key_value('F', sys.float_info.max)
            
        lines = [l for l in str(h).split('\n') if l.strip()]
        self.assertGreater(len(lines), 2)         # section + key + continuation(s)
        d = header.parse_header(str(h))
        self.assertEqual(d['P']['I'], sys.maxsize)
        self.assertEqual(d['P']['F'], sys.float_info.max)
        
        # The auto-sized default keeps the value whole and round-trips it
        h0 = header.CommentHeader(width=0)
        with h0.section('P'):
            h0.key_value('K', long_value)
            h0.key_value('I', sys.maxsize)
            h0.key_value('F', sys.float_info.max)
        d = header.parse_header(str(h0))
        self.assertEqual(d['P']['K'], long_value)
        self.assertEqual(d['P']['I'], sys.maxsize)
        self.assertEqual(d['P']['F'], sys.float_info.max)

    def test_namespace_default_title(self):
        """Test that namespace() files the flags under a 'Parameters' section."""

        ns = argparse.Namespace(x=1)
        h = header.CommentHeader(width=0)
        h.namespace(ns)

        d = header.parse_header(str(h))
        self.assertEqual(d, {'Parameters': {'x': 1}})

    def test_namespace_custom_title(self):
        """Test that namespace() honors a custom section title."""

        ns = argparse.Namespace(x=1)
        h = header.CommentHeader(width=0)
        h.namespace(ns, title='Settings')

        d = header.parse_header(str(h))
        self.assertIn('Settings', d)
        self.assertNotIn('Parameters', d)

    def test_namespace_roundtrip_types(self):
        """Test that every flag in a namespace round-trips with its type."""

        ns = argparse.Namespace(alpha=1, beta='two', gamma=True, delta=None,
                                epsilon=0.5)
        h = header.CommentHeader(width=0)
        h.namespace(ns, title='Settings')

        d = header.parse_header(str(h))
        self.assertEqual(d['Settings'], {'alpha': 1, 'beta': 'two',
                                         'gamma': True, 'delta': None,
                                         'epsilon': 0.5})

    def test_namespace_accepts_dict(self):
        """Test that namespace() accepts a dict as well as a Namespace and
        that the two render identically."""

        mapping = {'alpha': 1, 'beta': 'two', 'gamma': True, 'delta': None,
                   'eps': [0, 1, 2]}

        hd = header.CommentHeader(width=0)
        hd.namespace(mapping, title='Settings')
        d = header.parse_header(str(hd))
        self.assertEqual(d['Settings'], mapping)

        # A dict and the equivalent Namespace produce the same output
        hn = header.CommentHeader(width=0)
        hn.namespace(argparse.Namespace(**mapping), title='Settings')
        self.assertEqual(str(hd), str(hn))

    def test_columns_numbering_and_roundtrip(self):
        """Test that columns() numbers labels from one and round-trips."""

        labels = ['LST [hr]', 'Temperature [K]']
        h = header.CommentHeader(width=0)
        h.columns(labels)

        d = header.parse_header(str(h))
        self.assertEqual(d, {'Columns': {'1': 'LST [hr]',
                                         '2': 'Temperature [K]'}})

    def test_columns_unquoted(self):
        """Test that columns() renders labels without JSON quotation marks,
        unlike the key_value-based namespace()."""

        labels = ['LST [hr]', 'Temperature [K]']

        h = header.CommentHeader(width=0)
        h.columns(labels)
        text = str(h)
        self.assertIn('1: LST [hr]', text)
        self.assertNotIn('"LST [hr]"', text)

        # namespace() with the same content does quote the labels
        hn = header.CommentHeader(width=0)
        hn.namespace({'1': 'LST [hr]', '2': 'Temperature [K]'}, title='Columns')
        self.assertIn('"LST [hr]"', str(hn))

    def test_columns_custom_title(self):
        """Test that columns() honors a custom section title."""

        h = header.CommentHeader(width=0)
        h.columns(['A', 'B'], title='Fields')
        d = header.parse_header(str(h))
        self.assertIn('Fields', d)
        self.assertNotIn('Columns', d)

    def test_columns_units(self):
        """Test that column_units are appended in brackets, and an empty unit
        leaves a name bare."""

        h = header.CommentHeader(width=0)
        h.columns(['LST', 'Index', 'Temperature'],
                  column_units=['hr', '', 'K'])

        d = header.parse_header(str(h))
        self.assertEqual(d, {'Columns': {'1': 'LST [hr]',
                                         '2': 'Index',
                                         '3': 'Temperature [K]'}})

    def test_columns_units_length_mismatch(self):
        """Test that a units list of the wrong length is rejected."""

        h = header.CommentHeader(width=0)
        self.assertRaises(RuntimeError, h.columns, ['A', 'B'],
                          column_units=['only one'])

    def test_text_no_comment_prefix(self):
        """Test that TextHeader emits plain lines with no comment prefix."""

        h = header.TextHeader(width=0)
        with h.section('Parameters'):
            h.key_value('APARM', 'this')
        for line in str(h).split('\n'):
            if line:
                self.assertFalse(line.startswith('#'))

    def test_text_section_indent(self):
        """Test that a TextHeader section's contents are indented past it."""

        h = header.TextHeader(width=0)
        with h.section('Parameters'):
            h.key_value('APARM', 'this')
        lines = [l for l in str(h).split('\n') if l]
        head = lines[0]                   # 'Parameters:'
        item = lines[1]                   # ' APARM: ...'
        head_indent = len(head) - len(head.lstrip(' '))
        item_indent = len(item) - len(item.lstrip(' '))
        self.assertGreater(item_indent, head_indent)

    def test_text_roundtrip_nested(self):
        """Test that a TextHeader round-trips when parsed with comment=''."""

        h = header.TextHeader(width=0)
        with h.section('Comments'):
            with h.section('General'):
                h.key_value('A', 1)
            with h.section('Specific'):
                h.key_value('B', 2)
        h.key_value('Top', 3)

        d = header.parse_header(str(h), comment='')
        self.assertEqual(d, {'Comments': {'General': {'A': 1},
                                          'Specific': {'B': 2}},
                             'Top': 3})

    def test_text_namespace_roundtrip(self):
        """Test that namespace() on a TextHeader round-trips its values."""

        ns = argparse.Namespace(alpha=1, beta='two', gamma=True, delta=None)
        h = header.TextHeader(width=0)
        h.namespace(ns, title='Settings')

        d = header.parse_header(str(h), comment='')
        self.assertEqual(d['Settings'], {'alpha': 1, 'beta': 'two',
                                         'gamma': True, 'delta': None})


class header_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.common.header
    units tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(header_tests))


if __name__ == '__main__':
    unittest.main()

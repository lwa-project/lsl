"""
Module to help write simple structured headers to a terminal or file.

..versionadded:: 4.0.1
"""

import os
import sys
import json
import argparse
import warnings

from io import StringIO
from textwrap import fill as tw_fill
from contextlib import contextmanager
from datetime import datetime, timezone

from lsl.version import version as LSL_VERSION
from lsl.logger import LSL_LOGGER

from typing import Any, Optional, TextIO, List, Union, Dict

__all__ = ['CommentHeader', 'BoxHeader', 'TextHeader', 'parse_header']
__version__ = '0.1'


class _BaseHeader:
    """
    Base class to make quick and easy headers for text files.  This class
    supports two header styles: 'box' and 'left'.
    """
    
    def __init__(self, style: str, comment: str='#', width: Optional[int]=0, indent: int=1):
        if style not in ('box', 'left'):
            raise ValueError(f"Invalid header style '{style}', must be one of 'box' or 'left'")
            
        self.style = style
        self.comment = comment
        if width is None:
            width, _ = os.get_terminal_size()
        self.width = width
        self._curr_indent = indent    # Active indent level
        self._base_indent = indent   # Base indent level
        
        if self.width > 0:
            ewidth = self.effective_width
            if ewidth < 8:
                raise ValueError(f"Effective width of {ewidth} is less than 8 characters")
            elif ewidth < 16:
                LSL_LOGGER.warning(f"Effective width of {ewidth} is less than 16 characters")
                
        self._contents = []
        self.clear()
        
    @property
    def effective_width(self) -> int:
        """
        Effective width of the header, i.e., the width of the contents, after
        accounting for the comments, the indentation, and the trailing space.
        Lines longer than this will be wrapped with an additional space used
        on the subsequent lines.
        
        .. note:: This value is only accurate when `width>0`.  If `width<=0`
                  then the size of the header is auto-scaled to the longest
                  line.
        """
        
        ewidth = self.width - len(self.comment) - self._curr_indent
        if self.style == 'box':
            ewidth -= self._base_indent + len(self.comment)
        return ewidth
        
    def clear(self):
        """
        Reset to a clean state.
        """
        
        self._contents.clear()
        
    def rule(self):
        """
        Add a horizontal rule to the header.
        """
        
        self._contents.append({'type': 'rule',
                               'value': self.comment,
                               'indent': 0
                              })
        
    def line(self, contents: str):
        """
        Add a line with the given contents to the header.  This will be indented
        to the current level inside the header.
        """
        
        self._contents.append({'type': 'line',
                               'value': contents,
                               'indent': self._curr_indent
                              })
        
    def blank(self):
        """
        Add a single blank line to the header.
        """
        
        self._contents.append({'type': 'blank',
                               'value': '',
                               'indent': 0
                              })
        
    @contextmanager
    def section(self, name: str):
        """
        Start a new section in the header and increase the indent level as long
        as the section is open.
        
        .. note:: Section names that do not end in ':' have a ':' appended.
        """
        
        if not name.endswith(':'):
            name = name + ':'
        self._contents.append({'type': 'section',
                               'value': name,
                               'indent': self._curr_indent
                              })
        
        self._curr_indent += self._base_indent
        try:
            yield self
        finally:
            self._curr_indent -= self._base_indent
            
    def key_value(self, key: str, value: Any):
        """
        Add a "key: value"-style line to the header.  This will be indented to
        the current level inside the header.
        
        .. note:: Internally this uses `json.dumps()` to serialize the value.
        
        .. note:: For more explicit control over how the line is formatted see
                  `line()`.
        """
        
        self._contents.append({'type': 'key_value',
                               'key': key,
                               'raw_value': value,
                               'value': f"{key}: {json.dumps(value)}",
                               'indent': self._curr_indent
                              })
        
    def namespace(self, namespace: Union[Dict,argparse.Namespace], title: str='Parameters'):
        """
        Helper to make it easy to document the command line flags/arguments set
        for an argparse-powered script.  The contents of the `namespace` will be
        placed under a new section called `title`.
        
        .. note:: This also accepts a dictionary as well.
        """
        
        if isinstance(namespace, argparse.Namespace):
            ns = vars(namespace)
        else:
            ns = namespace
        with self.section(title):
            for k,v in ns.items():
                self.key_value(k, v)
                
    def timestamp(self, prefix: str='Timestamp', timespec: str='seconds',
                  tz: Optional[timezone]=None):
        """
        Add a ISO formatted timestamp line to the header.  If `tz` is provided
        then that timezone is used instead of UTC.
        """
        
        if tz is None:
            tz = timezone.utc
        t = datetime.now(tz=tz)
        ts= t.isoformat(timespec='seconds')
        self.key_value(prefix, ts)
        
    def lsl_version(self):
        """
        Add the full LSL version to the header.
        """
        
        self.key_value('LSL Version', LSL_VERSION)
        
    def write(self, fh: Optional[TextIO]=None):
        """
        Write the current contents of the header to `fh`.  If `fh` is None then
        the output is written to `sys.stdout`.
        """
        
        # Bail out if there's nothing to write
        if not self._contents:
            return
            
        # If width <= 0 then we are in auto-width mode.  Figure out what the
        # correct width is to fit the longest line entry
        width = self.width
        if width <= 0:
            width = max([len(str(entry['value']))+entry['indent'] for entry in self._contents])
            width += len(self.comment) + 2*self._base_indent
            if self.style == 'box':
                width += len(self.comment)
                
            if width > 80:
                LSL_LOGGER.warning(f"Computed header width is {width} characters")
                
        if fh is None:
            fh = sys.stdout
            
        # Add in the lead in part of the box
        if self.style == 'box':
            line = self.comment*(width//len(self.comment)+1)
            fh.write(line[:width]+'\n')
            
            line = self.comment + ' '*(width-2*len(self.comment)) + self.comment
            fh.write(line+'\n')
            
        # Main header contents
        for entry in self._contents:
            if entry['type'] == 'rule':
                line = entry['value']*(1+3)
                if self.style == 'box':
                    line = self.comment*(width//len(self.comment)+1)
                fh.write(line[:width]+'\n')
                
            elif entry['type'] == 'blank':
                line = self.comment
                if self.style == 'box':
                    line += ' '*(width-2*len(self.comment)) + self.comment
                fh.write(line+'\n')
                
            else:
                contents = str(entry['value'])
                
                indent = entry['indent']
                ewidth = width - len(self.comment) - indent
                if self.style == 'box':
                    ewidth -= self._base_indent + len(self.comment)
                    
                line = tw_fill(contents, width=ewidth, subsequent_indent='... ')
                lines = line.split('\n')
                for line in lines:
                    if self.style == 'box':
                        line = line.ljust(ewidth, ' ')
                    line = f"{self.comment}{' '*indent}{line}"
                    if self.style == 'box':
                        line += ' '*self._base_indent + self.comment
                    fh.write(line+'\n')
                    
        # Add in the lead out part of the box
        if self.style == 'box':
            line = self.comment + ' '*(width-2*len(self.comment)) + self.comment
            fh.write(line+'\n')
            
            line = self.comment*(width//len(self.comment)+1)
            fh.write(line[:width]+'\n')
            
    def __str__(self) -> str:
        """
        Return the current contents of the header as a string.
        """
        
        s = StringIO()
        self.write(s)
        return s.getvalue()


class CommentHeader(_BaseHeader):
    """
    Class to make comment-based header blocks like:
    
    # Parameters:
    #  APARM: this
    #  BPARM: that
    """
    
    def __init__(self, comment: str='#', width: Optional[int]=0, indent: int=1):
        super().__init__('left', comment=comment, width=width, indent=indent)


class BoxHeader(_BaseHeader):
    """
    Class to make comment-based header blocks with surrounding boxes like:
    
     
    ################################
    #                              #
    # Parameters:                  #
    #  APARM: this                 #
    #  BPARM: that                 #
    #                              #
    ################################
    """
    
    def __init__(self, comment: str='#', width: Optional[int]=0, indent: int=1):
        super().__init__('box', comment=comment, width=width, indent=indent)


class TextHeader(_BaseHeader):
    """
    Class to make text header blocks like:
    
    Parameters:
     APARM: this
     BPARM: that
    """
    
    def __init__(self, width: Optional[int]=0, indent: int=1):
        super().__init__('left', comment='', width=width, indent=indent)
        self._curr_indent = 0
        
    def rule(self):
        """
        Add a horizontal rule to the header.
        """
        
        _BaseHeader.rule(self)
        self._contents[-1]['value'] = '#'


def parse_header(contents: Union[str,List[str]], comment: str='#') -> Dict[str,Any]:
    """
    Given a header written by CommentHeader/BoxHeader, parse the header's
    contents and return it as a dictionary.
    
    .. note:: For fixed-width headers, i.e. those with `width>0`, the loading
              of wrapped lines may produce unexpected out, ranging from missing
              spaces in text to incorrectly parsed values.
              
    .. note:: Free-form line that do not fit within the expected "key: value"
              structure will be loaded into "_freeform" keys in the returned
              dictionary.
    """
    
    if isinstance(contents, str):
        contents = contents.split('\n')
        
    curr_indent = -1
    base_indent = -1
    
    header = {}
    active = [header,]
    for i,line in enumerate(contents):
        if not line.startswith(comment):
            continue
            
        line = line.replace(comment, '').rstrip()
        if not line:
            continue
            
        # Parse the current (non-empty) line.  If there's a colon then we
        # assume we are either in a section or a key/value pair.  Otherwise
        # this must be a comment-style line.
        LSL_LOGGER.debug(f"Working on '{line}'")
        next_indent = len(line) - len(line.lstrip())
        if line.find(':') != -1:
            next_key, next_value = line.split(':', 1)
            next_key = next_key.strip()
            next_value = next_value.strip()
        else:
            next_key = ''
            next_value = line.strip()
            
        if next_indent > curr_indent:
            ## The indent level is higher than the current value.  We have either
            ## entered a section or we are on a line that's been wrapped by
            ## `textwrap.fill()`.
            if curr_indent == -1:
                if comment == '':
                    if next_indent >= 1 and base_indent == -1:
                        LSL_LOGGER.debug(f"  Setting base indent to {next_indent}")
                        base_indent = next_indent
                else:
                    LSL_LOGGER.debug(f"  Setting base indent to {next_indent}")
                    base_indent = next_indent
            else:
                if comment == '':
                    if next_indent >= 1 and base_indent == -1:
                        LSL_LOGGER.debug(f"  Setting base indent to {next_indent}")
                        base_indent = next_indent
            curr_indent = next_indent
            
            if next_value == '':
                ### It's a section
                LSL_LOGGER.debug(f"  Starting new section '{next_key}'")
                active[-1][next_key] = {}
                active.append(active[-1][next_key])
            else:
                will_continue = False
                try:
                    next_line = contents[i+1]
                    if next_line.startswith(comment + ' '*curr_indent + '... '):
                        will_continue = True
                except IndexError:
                    pass
                    
                if next_key == '':
                    ### It's a continued line
                    LSL_LOGGER.debug(f"  Processing as continuation of last key/value pair")
                    
                    try:
                        last_item = list(active[-1].keys())[-1]
                    except IndexError:
                        last_item = '_freeform'
                        active[-1][last_item] = ''
                    active[-1][last_item] += next_value.replace('... ', '')
                    if not will_continue:
                        try:
                            value = json.loads(active[-1][last_item])
                            active[-1][last_item] = value
                        except json.decoder.JSONDecodeError:
                            pass
                else:
                    ### It's a key/value pair
                    LSL_LOGGER.debug(f"  Processing as key/value pair")
                    
                    if not will_continue:
                        try:
                            next_value = json.loads(next_value)
                        except json.decoder.JSONDecodeError:
                            pass
                    active[-1][next_key] = next_value
                    
        elif next_indent < curr_indent:
            ## The indent level is lower than the current value.  We have
            ## either left a section or have finished a continued line
            ## that's been wrapped by `textwrap.fill()`.
            for _ in range((curr_indent - next_indent) // base_indent):
                ### Leave the section if the indent level has decreased enough
                ### .. note:: This can't distiguish between a section ending
                ###           and a continued line ending!
                LSL_LOGGER.debug(f"  Active is now {active[-1]}")
                if len(active) > 1:
                    active.pop()
            curr_indent = next_indent
            
            if next_value == '':
                ## It's a section
                LSL_LOGGER.debug(f"  Starting new section '{next_key}'")
                active[-1][next_key] = {}
                active.append(active[-1][next_key])
            else:
                will_continue = False
                try:
                    next_line = contents[i+1]
                    if next_line.startswith(comment + ' '*curr_indent + '... '):
                        will_continue = True
                except IndexError:
                    pass
                    
                ## It's a key/value pair
                LSL_LOGGER.debug(f"  Processing as key/value pair")
                
                if not will_continue:
                    try:
                        next_value = json.loads(next_value)
                    except json.decoder.JSONDecodeError:
                        pass
                active[-1][next_key] = next_value
                
        else:
            ## No change in the indent level, continued line or a key/value
            ## pair.
            will_continue = False
            try:
                next_line = contents[i+1]
                if next_line.startswith(comment + ' '*curr_indent + '... '):
                    will_continue = True
            except IndexError:
                pass
                
            if next_key == '':
                ### It's a continued line
                LSL_LOGGER.debug(f"  Processing as continuation of last key/value pair")
                
                try:
                    last_item = list(active[-1].keys())[-1]
                except IndexError:
                    last_item = '_freeform'
                    active[-1][last_item] = ''
                active[-1][last_item] += next_value.replace('... ', '')
                if not will_continue:
                    try:
                        value = json.loads(active[-1][last_item])
                        active[-1][last_item] = value
                    except json.decoder.JSONDecodeError:
                        pass
            else:
                ## It's a key/value pair
                LSL_LOGGER.debug(f"  Processing as key/value pair")
                
                if not will_continue:
                    try:
                        next_value = json.loads(next_value)
                    except json.decoder.JSONDecodeError:
                        pass
                active[-1][next_key] = next_value
                
    return header
                

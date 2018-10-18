"""
Module to help display color on the command line of ANSI-compliant
termainals.

..versionadded:: 1.2.1
"""

import re

__version__ = '0.1'
__revision__ = '$Rev$'
__all__ = ['colorfy',]


# Dictionary of ANSI color codes
_COLORS = {'black':       u'\u001b[30m', 
           'red':         u'\u001b[31m', 
           'green':       u'\u001b[32m', 
           'yellow':      u'\u001b[33m', 
           'blue':        u'\u001b[34m', 
           'magenta':     u'\u001b[35m', 
           'cyan':        u'\u001b[36m', 
           'white':       u'\u001b[37m', 
           'bkg-black':   u'\u001b[40m', 
           'bkg-red':     u'\u001b[41m', 
           'bkg-green':   u'\u001b[42m', 
           'bkg-yellow':  u'\u001b[43m', 
           'bkg-blue':    u'\u001b[44m', 
           'bkg-magenta': u'\u001b[45m', 
           'bkg-cyan':    u'\u001b[46m', 
           'bkg-white':   u'\u001b[47m', 
           'reset':       u'\u001b[0m', 
           'bold':        u'\u001b[1m', 
           'underline':   u'\u001b[4m', 
           'blink':       u'\u001b[5m',
           'reverse':     u'\u001b[7m'}


def _label(text):
    """
    Private function to tag a ANSI tagged string and split it into various
    tag closure levels.  Returns the level-tagged string.
    """
    
    # Setup the control variables for position in the string and nesting depth
    level = 0
    pos = 0
    
    # Go, until we run out of string
    while pos < len(text):
        ## Split the string at the current position so that we can search for
        ## '{{' and '}}'
        pre  = text[:pos]
        tag  = text[pos:pos+2]
        post = text[pos+2:]
        
        ## Search for tag opening/closing - For the opening we need to make sure 
        ## there is a '%' immediately afterwards.  For closing we need to make 
        ## sure we actually have a valid level to close.
        if tag == '{{':
            if post[0] == '%':
                tag = '{{'+str(level)+tag
                level += 1
                pos = pos + 4
        if tag == '}}':
            level -=1
            if level >= 0:
                tag = tag+str(level)+'}}'
                pos = pos + 4
        text = pre+tag+post
        pos += 1
        
    # Done
    return text


def colorfy(text):
    """
    Given a string encoded with color and/or style information, return a string with
    embedded ANSI terminal codes.  For example:
    {{%red text}} {{%underline underline}}
    will return the text 'text underline' with 'text' displayed as red and 'underline'
    underlined.
    
    Valid colors are:
     * black
     * red
     * green
     * yellow
     * blue
     * magenta
     * cyan
     * white
    
    Valid type styles are:
     * bold
     * underline
     * blink
     * reverse
    
    .. note::
        Background colors can be set by appending 'bkg-' to the color name.
    """
    
    # Parse 1 - Convert nested tags into numbered tags
    text = _label(text)
    
    # Parse 2 - Find the maximum nesting depth of the tags
    nestRE = re.compile(r'{{(?P<depth>\d+){{%')
    levels = [int(l, 10) for l in nestRE.findall(text)]
    levels.append(0)
    
    # Parse 3 - Work at each level to convert the tags to ANSI sequences
    for level in xrange(max(levels), -1, -1):
        ## Build the regular expression needed for this depth
        tagRE = re.compile('{{'+str(level)+'{{%(?P<tag>[a-zA-Z-]+) (?P<text>.*?)}}'+str(level)+'}}')
        
        ## Match and replace until there is nothing left to match
        pos = 0
        mtch = tagRE.search(text, pos)
        while mtch is not None and pos < len(text):
            ### Parse the match
            tag, txt = mtch.group('tag'), mtch.group('text')
            tag = tag.lower()
            try:
                tag = _COLORS[tag]
            except KeyError:
                tag = ''
                
            ### Update the search start position so we can move on down the string
            pos = mtch.start() + 1
            
            ### Substitute, taking into account any higher levels that have already been encoded
            repl = tag+txt
            repl = repl.replace(_COLORS['reset'], _COLORS['reset']+tag)
            repl = repl+_COLORS['reset']
            text = tagRE.sub(repl, text, count=1)
            
            ### Get the next match
            mtch = tagRE.search(text, pos)
            
    # Done
    return text

# -*- coding: utf-8 -*-

"""
Module that provides argparse-compatible conversion functions for a variety 
of value formats, including:
 * positive integers, 
 * ephem.hours instances, and
 * lists of integers from a comma-separated list.

.. versionadded:: 1.2.4
"""

import re
import ephem
from argparse import ArgumentTypeError
from astropy import units

__version__ = '0.1'
__revision__ = '$Rev$'
__all__ = ['positive_or_zero_int', 'positive_int', 'positive_or_zero_float', 
           'positive_float', 'frequency', 'frequency_range', 'wavelength', 
           'wavelength_range', 'hours', 'csv_hours_list', 'degrees', 
           'csv_degrees_list', 'csv_int_list', 'csv_baseline_list', 'csv_hostname_list']


def positive_or_zero_int(string):
    """
    Convert a string to a positive (>=0) integer.
    """
    
    try:
        value = int(string, 10)
    except ValueError:
        msg = "%r is a non-integer value" % string
        raise ArgumentTypeError(msg)
    if value < 0:
        msg = "%r < 0" % string
        raise ArgumentTypeError(msg)
    return value


def positive_int(string):
    """
    Convert a string to a positive (>0) integer.
    """
    
    value = positive_or_zero_int(string)
    if value <= 0:
        msg = "%r <= 0" % string
        raise ArgumentTypeError(msg)
    return value


def positive_or_zero_float(string):
    """
    Convert a string to a positive (>=0.0) float.
    """
    
    try:
        value = float(string)
    except ValueError:
        msg = "%r is a non-float value" % string
        raise ArgumentTypeError(msg)
    if value < 0.0:
        msg = "%r < 0.0" % string
        raise ArgumentTypeError(msg)
    return value


def positive_float(string):
    """
    Convert a string to a positive (>0.0) float.
    """
    
    value = positive_or_zero_float(string)
    if value <= 0.0:
        msg = "%r <= 0.0" % string
        raise ArgumentTypeError(msg)
    return value


def _get_units(string):
    """
    Function to search a string, starting at the end, to find units.
    """
    
    units = None
    for i in xrange(len(string), 0, -1):
        try:
            float(string[:i])
            units = string[i:]
            if units == '':
                units = None
            break
        except:
            pass
    return units


def _quantitiy_to_hz(value):
    """
    Convert a string into a frequency.  If no unit is provided, MHz is 
    assumed.
    """
    
    try:
        value = float(value)
        value *= 1e6
    except ValueError:
        try:
            value = units.quantity.Quantity(value)
            value = value.to(units.Hz, equivalencies=units.spectral()).value
        except Exception as e:
            msg = "%r %s" % (value, str(e))
            raise ArgumentTypeError(msg)
    return value


def _quantitiy_to_m(value):
    """
    Convert a string into a wavelength.  If no unit is provided, m is 
    assumed.
    """
    
    try:
        value = float(value)
    except ValueError:
        try:
            value = units.quantity.Quantity(value)
            value = value.to(units.meter, equivalencies=units.spectral()).value
        except Exception as e:
            msg = "%r %s" % (value, str(e))
            raise ArgumentTypeError(msg)
    return value


def _frequency_conversion_base(string):
    """
    Convert a frequency to a float Hz value.  This function accepts a variety 
    of string formats:
     * pure float values are intepreted to be in MHz (45.0 -> 45e6)
     * number/unit pairs are allowed so long as they are in:
        * [prefix]m, AA, or Angstrom for wavelength and 
        * [prefix]Hz for frequency
     * a 'number~number' is interpretted as a range in MHz
     * a 'number/unit~number/unit' is converted to a range in Hz
    
    .. note::
        For ranges, a two-element list is returned where the first value
        is less than the second.
    """
    
    try:
        value = float(string)*1e6
    except ValueError:
        fields = string.split('~', 1)
        try:
            start, stop = fields
            units1 = _get_units(start)
            units2 = _get_units(stop)
            print units1, units2
            if units1 is not None and units2 is None:
                msg = "%r must have units specified for the second value" % string
                raise ArgumentTypeError(msg)
            elif units2 is not None and units1 is None:
                start = start+units2
        except ValueError:
            start, stop = fields[0], None
        value = _quantitiy_to_hz(start)
        if stop is not None:
            value = [value, _quantitiy_to_hz(stop)]
            if value[1] < value[0]:
                value.reverse()
    return value


def frequency(string):
    """
    Convert a frequency to a float Hz value.  This function accepts a variety 
    of string formats:
     * pure float values are intepreted to be in MHz (45.0 -> 45e6)
     * number/unit pairs are allowed so long as they are in:
        * [prefix]m, A, or ang for wavelength and 
        * [prefix]Hz for frequency
    """
    
    value = _frequency_conversion_base(string)
    try:
        len(value)
        msg = "%r does not appear to be a single frequency" % string
        raise ArgumentTypeError(msg)
    except TypeError:
        pass
    return value


def frequency_range(string):
    """
    Convert a frequency to a float Hz value.  This function accepts a variety 
    of string formats:
     * a 'number~number' is interpretted as a range in MHz
     * a 'number/unit~number/unit' is converted to a range in Hz
    
    .. note::
        For ranges, a two-element list is returned where the first value
        is less than the second.
    """
    
    value = _frequency_conversion_base(string)
    try:
        len(value)
    except TypeError:
        msg = "%r does not appear to be a frequency range" % string
        raise ArgumentTypeError(msg)
    return value


def _wavelength_conversion_base(string):
    """
    Convert a wavelength to a float m value.  This function accepts a variety 
    of string formats:
     * pure float values are intepreted to be in m (45.0 -> 45.0)
     * number/unit pairs are allowed so long as they are in:
        * [prefix]m, A, or ang for wavelength and 
        * [prefix]Hz for frequency
     * a 'number~number' is interpretted as a range in m
     * a 'number/unit~number/unit' is converted to a range in m
    
    .. note::
        For ranges, a two-element list is returned where the first value
        is less than the second.
    """
    
    try:
        value = float(string)
    except ValueError:
        fields = string.split('~', 1)
        try:
            start, stop = fields
            units1 = _get_units(start)
            units2 = _get_units(stop)
            if units1 is not None and units2 is None:
                msg = "%r must have units specified for the second value" % string
                raise ArgumentTypeError(msg)
            elif units2 is not None and units1 is None:
                start = start+units2
        except ValueError:
            start, stop = fields[0], None
        value = _quantitiy_to_m(start)
        if stop is not None:
            value = [value, _quantitiy_to_m(stop)]
            if value[1] < value[0]:
                value.reverse()
    return value


def wavelength(string):
    """
    Convert a wavelength to a float m value.  This function accepts a variety 
    of string formats:
     * pure float values are intepreted to be in m (45.0 -> 45.0)
     * number/unit pairs are allowed so long as they are in:
        * [prefix]m, AA, or Angstrom for wavelength and 
        * [prefix]Hz for frequency
    """
    
    value = _wavelength_conversion_base(string)
    try:
        len(value)
        msg = "%r does not appear to be a single wavelength" % string
        raise ArgumentTypeError(msg)
    except TypeError:
        pass
    return value


def wavelength_range(string):
    """
    Convert a wavelength to a float m value.  This function accepts a variety 
    of string formats:
     * a 'number~number' is interpretted as a range in m
     * a 'number/unit~number/unit' is converted to a range in m
    
    .. note::
        For ranges, a two-element list is returned where the first value
        is less than the second.
    """
    
    value = _wavelength_conversion_base(string)
    try:
        len(value)
    except TypeError:
        msg = "%r does not appear to be a wavelength range" % string
        raise ArgumentTypeError(msg)
    return value


def hours(string):
    """
    Convert a 'HH[:MM[:SS[.SSS]]]' string into an ephem.hours instance.
    """
    
    try:
        value = ephem.hours(string)
    except ValueError as e:
        msg = str(e) % string
        raise ArgumentTypeError(msg)
    return value


def csv_hours_list(string):
    """
    Convert a comma-separated list of 'HH[:MM[:SS.[SSS]]]' string into a list 
    of ephem.hours instances.
    """
    
    string = string.rstrip()
    if string[-1] == ',':
        string = string[:-1]
        
    value = []
    for item in string.split(','):
        value.append( hours(item) )
    return value


def degrees(string):
    """
    Convert a 'sDD[:MM[:SS[.SSS]]]' string into an ephem.degrees instance.
    """
    
    try:
        value = ephem.degrees(string)
    except ValueError as e:
        msg = str(e) % string
        raise ArgumentTypeError(msg)
    return value


def csv_degrees_list(string):
    """
    Convert a comma-separated list of 'sDD[:MM[:SS.[SSS]]]' string into a list 
    of ephem.degrees instances.
    """
    
    string = string.rstrip()
    if string[-1] == ',':
        string = string[:-1]
        
    value = []
    for item in string.split(','):
        value.append( degrees(item) )
    return value


def _int_item_or_range(string):
    if string.find('~') != -1:
        start, stop = string.split('~', 1)
        start, stop = int(start, 10), int(stop, 10)
        value = list(range(start, stop+1))
    else:
        value = [int(string, 10),]
    return value


def csv_int_list(string):
    """
    Convert a comma-separated list of integers into a list of integers.  This 
    function also allows for ranges to be specifed using the '~' character.  
    This formatting converts 'start~stop' to 'range(start, stop+1)'.
    """
    
    if string in ('all', '*'):
        value = 'all'
    elif string in ('none', ''):
        value = 'none'
    else:
        value = []
        for item in string.split(','):
            item = item.strip().rstrip()
            if item == '':
                continue
            try:
                subvalue = _int_item_or_range(item)
            except ValueError:
                msg = "%r contains non-integer values" % string
                raise ArgumentTypeError(msg)
            value.extend(subvalue)
    return value


def csv_baseline_list(string):
    """
    Convert a comma-separated list of baslines pairs into a list of baseline
    pairs.  Baseline pairs are defined as 'antenna1-antenna2' where 'antenna1'
    and 'antenna2' are both integers or ranges denoted by the '~' character.
    """
    
    if string in ('all', '*'):
        value = 'all'
    elif string in ('none', ''):
        value = 'none'
    else:
        value = []
        for item in string.split(','):
            item = item.strip().rstrip()
            if item == '':
                continue
            try:
                ant1, ant2 = item.split('-', 1)
            except ValueError:
                msg = "%s contains non-baseline or non-integer values" % string
                raise ArgumentTypeError(msg)
            try:
                ant1 = _int_item_or_range(ant1)
                ant2 = _int_item_or_range(ant2)
            except ValueError:
                msg = "%r contains non-integer values" % string
                raise ArgumentTypeError(msg)
            for i in ant1:
                for j in ant2:
                    value.append( (i,j) )
    return value


_IPV4_RANGE_RE = re.compile('^(?P<byte1>\d{1,3})(~(?P<byte1e>\d{1,3}))?\.(?P<byte2>\d{1,3})(~(?P<byte2e>\d{1,3}))?\.(?P<byte3>\d{1,3})(~(?P<byte3e>\d{1,3}))?\.(?P<byte4>\d{1,3})(~(?P<byte4e>\d{1,3}))?$')

_HOSTNAME_RE = re.compile('^(?P<hostname>[a-zA-Z\-0-9]*?)$')
_HOSTNAME_RANGE_RE=re.compile('^(?P<hostbase>[a-zA-Z\-]*?)(?P<start>[0-9]+)~(?P<stop>[0-9]+)$')
    

def csv_hostname_list(string):
    """
    Convert a comma-separated list of IPv4 addresses/hostnames into a list 
    IPv4 addresses/hostnames.  This function support indexing with the '~' 
    character so long as:
     * the character is in any one of the IPv4 bytes or
     * the character is at the end of a hostname which ends with a number
    """
    
    value = []
    for item in string.split(','):
        item = item.strip().rstrip()
        if item == '':
            continue
        mtch = _IPV4_RANGE_RE.match(item)
        if mtch is not None:
            ## IPv4 address or IPv4 address range
            b1b = int(mtch.group('byte1'), 10)
            b1e = mtch.group('byte1e')
            b1e = int(b1e, 10) if b1e is not None else b1b
            
            b2b = int(mtch.group('byte2'), 10)
            b2e = mtch.group('byte2e')
            b2e = int(b2e, 10) if b2e is not None else b2b
            
            b3b = int(mtch.group('byte3'), 10)
            b3e = mtch.group('byte3e')
            b3e = int(b3e, 10) if b3e is not None else b3b
            
            b4b = int(mtch.group('byte4'), 10)
            b4e = mtch.group('byte4e')
            b4e = int(b4e, 10) if b4e is not None else b4b
            
            items = []
            for b1 in range(b1b, b1e+1):
                for b2 in range(b2b, b2e+1):
                    for b3 in range(b3b, b3e+1):
                        for b4 in range(b4b, b4e+1):
                            items.append( '%i.%i.%i.%i' % (b1, b2, b3, b4) )
        else:
            mtch = _HOSTNAME_RANGE_RE.match(item)
            if mtch is not None:
                ## Hostname range
                hostbase = mtch.group('hostbase')
                try:
                    start = int(mtch.group('start'), 10)
                    stop = int(mtch.group('stop'), 10)
                except ValueError:
                    msg = "%r contains non-integer hostname values" % string
                    raise ArgumentTypeError(msg)
                items = ['%s%i' % (hostbase, i) for i in range(start, stop+1)]
            else:
                mtch = _HOSTNAME_RE.match(item)
                if mtch is not None:
                    ## Single hostname
                    items = [mtch.group('hostname'),]
                else:
                    msg = "%r contains invalid hostname values" % string
                    raise ArgumentTypeError(msg)
        value.extend( items )
    return value

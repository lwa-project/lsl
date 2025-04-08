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
from datetime import datetime
from astropy import units

from lsl.common.mcs import datetime_to_mjdmpm, mjdmpm_to_datetime

from lsl.misc import telemetry
telemetry.track_module()

from typing import List, Tuple, Union


__version__ = '0.1'
__all__ = ['positive_or_zero_int', 'positive_int', 'positive_or_zero_float', 
           'positive_float', 'frequency', 'frequency_range', 'wavelength', 
           'wavelength_range', 'date', 'mjd', 'time', 'mpm', 'hours', 
           'csv_hours_list', 'degrees', 'csv_degrees_list', 'csv_int_list', 
           'csv_baseline_list', 'csv_hostname_list']


def positive_or_zero_int(string: str) -> int:
    """
    Convert a string to a positive (>=0) integer.
    """
    
    try:
        value = int(string, 10)
    except ValueError:
        msg = f"{string} is a non-integer value"
        raise ArgumentTypeError(msg)
    if value < 0:
        msg = f"{string} < 0"
        raise ArgumentTypeError(msg)
    return value


def positive_int(string: str) -> int:
    """
    Convert a string to a positive (>0) integer.
    """
    
    value = positive_or_zero_int(string)
    if value <= 0:
        msg = f"{string} <= 0"
        raise ArgumentTypeError(msg)
    return value


def positive_or_zero_float(string: str) -> float:
    """
    Convert a string to a positive (>=0.0) float.
    """
    
    try:
        value = float(string)
    except ValueError:
        msg = f"{string} is a non-float value"
        raise ArgumentTypeError(msg)
    if value < 0.0:
        msg = f"{string} < 0.0"
        raise ArgumentTypeError(msg)
    return value


def positive_float(string: str) -> float:
    """
    Convert a string to a positive (>0.0) float.
    """
    
    value = positive_or_zero_float(string)
    if value <= 0.0:
        msg = f"{string} <= 0.0"
        raise ArgumentTypeError(msg)
    return value


def _get_units(string: str) -> Union[str,None]:
    """
    Function to search a string, starting at the end, to find units.
    """
    
    units = None
    for i in range(len(string), 0, -1):
        try:
            float(string[:i])
            units = string[i:]
            if units == '':
                units = None
            break
        except:
            pass
    return units


def _quantitiy_to_hz(value: str) -> float:
    """
    Convert a string into a frequency.  If no unit is provided, MHz is 
    assumed.
    """
    
    try:
        pvalue = float(value)
        pvalue *= 1e6
    except ValueError:
        try:
            qvalue = units.quantity.Quantity(value)
            pvalue = qvalue.to(units.Hz, equivalencies=units.spectral()).value
        except Exception as e:
            msg = f"{value} {str(e)}"
            raise ArgumentTypeError(msg)
    return pvalue


def _quantitiy_to_m(value: str) -> float:
    """
    Convert a string into a wavelength.  If no unit is provided, m is 
    assumed.
    """
    
    try:
        pvalue = float(value)
    except ValueError:
        try:
            qvalue = units.quantity.Quantity(value)
            pvalue = qvalue.to(units.meter, equivalencies=units.spectral()).value
        except Exception as e:
            msg = f"{value} {str(e)}"
            raise ArgumentTypeError(msg)
    return pvalue


def _frequency_conversion_base(string: str) -> Union[float,List[float]]:
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
        pvalue = float(string)*1e6
    except ValueError:
        fields = string.split('~', 1)
        try:
            start, stop = fields
            units1 = _get_units(start)
            units2 = _get_units(stop)
            if units1 is not None and units2 is None:
                msg = f"{string} must have units specified for the second value"
                raise ArgumentTypeError(msg)
            elif units2 is not None and units1 is None:
                start = start+units2
        except ValueError:
            start, stop = fields[0], None
        pvalue = _quantitiy_to_hz(start)
        if stop is not None:
            rvalue = [pvalue, _quantitiy_to_hz(stop)]
            if rvalue[1] < rvalue[0]:
                rvalue.reverse()
            pvalue = rvalue     # type: ignore
    return pvalue


def frequency(string: str) -> float:
    """
    Convert a frequency to a float Hz value.  This function accepts a variety 
    of string formats:
     * pure float values are intepreted to be in MHz (45.0 -> 45e6)
     * number/unit pairs are allowed so long as they are in:
        * [prefix]m, A, or ang for wavelength and 
        * [prefix]Hz for frequency
    """
    
    value = _frequency_conversion_base(string)
    if isinstance(value, List):
        msg = f"{string} does not appear to be a single frequency"
        raise ArgumentTypeError(msg)
    return value


def frequency_range(string: str) -> List[float]:
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
    if not isinstance(value, List):
        msg = f"{string} does not appear to be a frequency range"
        raise ArgumentTypeError(msg)
    return value


def _wavelength_conversion_base(string: str) -> Union[float,List[float]]:
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
        pvalue = float(string)
    except ValueError:
        fields = string.split('~', 1)
        try:
            start, stop = fields
            units1 = _get_units(start)
            units2 = _get_units(stop)
            if units1 is not None and units2 is None:
                msg = f"{string} must have units specified for the second value"
                raise ArgumentTypeError(msg)
            elif units2 is not None and units1 is None:
                start = start+units2
        except ValueError:
            start, stop = fields[0], None
        pvalue = _quantitiy_to_m(start)
        if stop is not None:
            rvalue = [pvalue, _quantitiy_to_m(stop)]
            if rvalue[1] < rvalue[0]:
                rvalue.reverse()
            pvalue = rvalue     # type: ignore
    return pvalue


def wavelength(string: str) -> float:
    """
    Convert a wavelength to a float m value.  This function accepts a variety 
    of string formats:
     * pure float values are intepreted to be in m (45.0 -> 45.0)
     * number/unit pairs are allowed so long as they are in:
        * [prefix]m, AA, or Angstrom for wavelength and 
        * [prefix]Hz for frequency
    """
    
    value = _wavelength_conversion_base(string)
    if isinstance(value, List):
        msg = f"{string} does not appear to be a single wavelength"
        raise ArgumentTypeError(msg)
    return value


def wavelength_range(string: str) -> List[float]:
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
    if not isinstance(value, List):
        msg = f"{string} does not appear to be a wavelength range"
        raise ArgumentTypeError(msg)
    return value


def date(string: str) -> str:
    """
    Convert a data as either a YYYY[-/]MM[-/]DD or MJD string into a 
    YYYY/MM/DD string.
    """
    
    try:
        mjd = int(string, 10)
        dt = mjdmpm_to_datetime(mjd, 0)
    except ValueError:
        cstring = string.replace('-', '/')
        try:
            dt = datetime.strptime("%s 00:00:00" % cstring, "%Y/%m/%d %H:%M:%S")
        except ValueError:
            msg = f"{string} cannot be interpretted as an MJD or date string"
            raise ArgumentTypeError(msg)
            
    date = dt.strftime('%Y/%m/%d')
    return date


def mjd(string: str) -> int:
    """
    Convert a data as either a YYYY[-/]MM[-/]DD or MJD string into an integer
    MJD.
    """
    
    try:
        mjd = int(string, 10)
    except ValueError:
        cstring = string.replace('-', '/')
        try:
            dt = datetime.strptime("%s 00:00:00" % cstring, "%Y/%m/%d %H:%M:%S")
            mjd, mpm = datetime_to_mjdmpm(dt)
        except ValueError:
            msg = f"{string} cannot be interpretted as an MJD or date string"
            raise ArgumentTypeError(msg)
            
    return mjd


def time(string: str) -> str:
    """
    Covnert a time as HH:MM:SS[.SSS] or MPM string into a HH:MM:SS.SSSSSS 
    string.
    """
    
    try:
        mpm = int(string, 10)
        if mpm < 0 or mpm > (86400*1000 + 999):
            msg = f"{mpm} is out of range for an MPM value"
            raise ArgumentTypeError(msg)
        s, f = mpm/1000, mpm%1000
        h = s / 3600
        m = (s / 60) % 60
        s = s % 60
        stime = "%i:%02i:%02i.%06i" % (h, m, s, f*1000)
    except ValueError:
        try:
            dt = datetime.strptime("2000/1/1 %s" % string, "%Y/%m/%d %H:%M:%S.%f")
        except ValueError:
            try:
                dt = datetime.strptime("2000/1/1 %s" % string, "%Y/%m/%d %H:%M:%S")
            except ValueError:
                msg = f"{string} cannot be interpretted as a time string"
                raise ArgumentTypeError(msg)
        stime = "%i:%02i:%02i.%06i" % (dt.hour, dt.minute, dt.second, dt.microsecond)
    return stime


def mpm(string: str) -> int:
    """
    Covnert a time as HH:MM:SS[.SSS] or MPM string into an MPM integer.
    """
    
    try:
        mpm = int(string, 10)
        if mpm < 0 or mpm > (86400*1000 + 999):
            msg = f"{mpm} is out of range for an MPM value"
            raise ArgumentTypeError(msg)
    except ValueError:
        try:
            dt = datetime.strptime("2000/1/1 %s" % string, "%Y/%m/%d %H:%M:%S.%f")
        except ValueError:
            try:
                dt = datetime.strptime("2000/1/1 %s" % string, "%Y/%m/%d %H:%M:%S")
            except ValueError:
                msg = f"{string} cannot be interpretted as a time string"
                raise ArgumentTypeError(msg)
        mjd, mpm = datetime_to_mjdmpm(dt)
    return mpm


def hours(string: str) -> ephem.hours:
    """
    Convert a 'HH[:MM[:SS[.SSS]]]' string into an ephem.hours instance.
    """
    
    try:
        value = ephem.hours(string)
    except ValueError as e:
        msg = f"{str(e)}: {string}"
        raise ArgumentTypeError(msg)
    return value


def csv_hours_list(string:str) -> List[ephem.hours]:
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


def degrees(string: str) -> ephem.degrees:
    """
    Convert a 'sDD[:MM[:SS[.SSS]]]' string into an ephem.degrees instance.
    """
    
    try:
        value = ephem.degrees(string)
    except ValueError as e:
        msg = f"{str(e)}: {string}"
        raise ArgumentTypeError(msg)
    return value


def csv_degrees_list(string: str) -> List[ephem.degrees]:
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


def _int_item_or_range(string: str) -> List[int]:
    if string.find('~') != -1:
        start, stop = string.split('~', 1)
        startv, stopv = int(start, 10), int(stop, 10)
        value = list(range(startv, stopv+1))
    else:
        value = [int(string, 10),]
    return value


def csv_int_list(string: str) -> Union[List[int],str]:
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
        lvalue = []
        for item in string.split(','):
            item = item.strip().rstrip()
            if item == '':
                continue
            try:
                subvalue = _int_item_or_range(item)
            except ValueError:
                msg = f"{string} contains non-integer values"
                raise ArgumentTypeError(msg)
            lvalue.extend(subvalue)
        value = lvalue  # type: ignore
    return value


def csv_baseline_list(string: str) -> Union[List[Tuple[int,int]],str]:
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
        lvalue = []
        for item in string.split(','):
            item = item.strip().rstrip()
            if item == '':
                continue
            try:
                ant1, ant2 = item.split('-', 1)
            except ValueError:
                msg = f"{string} contains non-baseline or non-integer values"
                raise ArgumentTypeError(msg)
            try:
                ant1l = _int_item_or_range(ant1)
                ant2l = _int_item_or_range(ant2)
            except ValueError:
                msg = f"{string} contains non-integer values"
                raise ArgumentTypeError(msg)
            for i in ant1l:
                for j in ant2l:
                    lvalue.append( (i,j) )
        value = lvalue      # type: ignore
    return value


_IPV4_RANGE_RE = re.compile(r'^(?P<byte1>\d{1,3})(~(?P<byte1e>\d{1,3}))?\.(?P<byte2>\d{1,3})(~(?P<byte2e>\d{1,3}))?\.(?P<byte3>\d{1,3})(~(?P<byte3e>\d{1,3}))?\.(?P<byte4>\d{1,3})(~(?P<byte4e>\d{1,3}))?$')

_HOSTNAME_RE = re.compile(r'^(?P<hostname>[a-zA-Z0-9-]*?)$')
_HOSTNAME_RANGE_RE=re.compile(r'^(?P<hostbase>[a-zA-Z-]*?)(?P<start>[0-9]+)~(?P<stop>[0-9]+)$')
    

def csv_hostname_list(string: str) -> List[str]:
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
                    msg = f"{string} contains non-integer hostname values"
                    raise ArgumentTypeError(msg)
                items = ['%s%i' % (hostbase, i) for i in range(start, stop+1)]
            else:
                mtch = _HOSTNAME_RE.match(item)
                if mtch is not None:
                    ## Single hostname
                    items = [mtch.group('hostname'),]
                else:
                    msg = f"{string} contains invalid hostname values"
                    raise ArgumentTypeError(msg)
        value.extend( items )
    return value

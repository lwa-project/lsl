"""
Module for managing configuration parameters for LSL/
"""

import os
import copy
import contextlib
from datetime import datetime
from collections import OrderedDict
from textwrap import fill as tw_fill

# Create the .lsl directory and set the config filename
if not os.path.exists(os.path.join(os.path.expanduser('~'), '.lsl')):
    os.mkdir(os.path.join(os.path.expanduser('~'), '.lsl'))
_CONFIG_FILENAME = os.path.join(os.path.expanduser('~'), '.lsl', 'lsl.cfg')


__version__ = "0.1"
__all__ = ['LSL_CONFIG',]


# Default values
## lsl.misc.ionosphere
DEFAULTS_IONO = OrderedDict()
DEFAULTS_IONO['igs_url'] = {'value': 'ftp://gssc.esa.int/gnss/products/ionex/',
                            'help':  'primary URL for accessing the IGS data products'}
DEFAULTS_IONO['igs_mirror'] = {'value': 'ftp://igs.ensg.ign.fr/pub/igs/products/ionosphere/',
                               'help':  'mirror URL for accessing the IGS data products'}
DEFAULTS_IONO['jpl_url'] = {'value': 'ftp://gssc.esa.int/gnss/products/ionex/',
                            'help':  'primary URL for accessing the JPL data products'}
DEFAULTS_IONO['jpl_mirror'] = {'value': 'ftp://igs.ensg.ign.fr/pub/igs/products/ionosphere/',
                               'help':  'mirror URL for accessing the JPL data products'}
DEFAULTS_IONO['uqr_url'] = {'value': 'ftp://gssc.esa.int/gnss/products/ionex/',
                            'help':  'primary URL for accessing the UQR data products'}
DEFAULTS_IONO['uqr_mirror'] = {'value': 'ftp://igs.ensg.ign.fr/pub/igs/products/ionosphere/',
                               'help':  'mirror URL for accessing the UQR data products'}
DEFAULTS_IONO['code_url'] = {'value': 'ftp://gssc.esa.int/gnss/products/ionex/',
                             'help':  'primary URL for accessing the CODE data products'}
DEFAULTS_IONO['code_mirror'] = {'value': 'ftp://igs.ensg.ign.fr/pub/igs/products/ionosphere/',
                                'help':  'mirror URL for accessing the CODE data products'}
DEFAULTS_IONO['ustec_url'] = {'value': 'http://www.ngdc.noaa.gov/stp/iono/ustec/products/',
                              'help':  'primary URL for accessing the USTEC data products'}
DEFAULTS_IONO['ustec_mirror'] = {'value': None,
                                 'help':  'mirror URL for accessing the USTEC data products'}
DEFAULTS_IONO['max_cache_size'] = {'value': -1,
                                   'help':  'maximum cache size in MB; <= 0 disables cache size limiting'}

## lsl.misc.telemetry
DEFAULTS_TELEMETRY = OrderedDict()
DEFAULTS_TELEMETRY['max_entries'] = {'value': 100,
                                     'help':  'maximum number of entries to accumlate before reporting'}
DEFAULTS_TELEMETRY['timeout'] = {'value': 120,
                                 'help':  'upload timeout in seconds'}

## Download parameters
DEFAULTS_DOWNLOAD = OrderedDict()
DEFAULTS_DOWNLOAD['timeout'] = {'value': 120,
                                'help':  'download timeout in seconds'}

## Everything
DEFAULTS_ALL = OrderedDict()
DEFAULTS_ALL['ionosphere'] = DEFAULTS_IONO
DEFAULTS_ALL['telemetry'] = DEFAULTS_TELEMETRY
DEFAULTS_ALL['download'] = DEFAULTS_DOWNLOAD


def _build_repr(name, attrs=[]):
    name = '.'.join(name.split('.')[-2:])
    output = "<%s" % name
    first = True
    for key,value in attrs:
        output += "%s %s=%s" % (('' if first else ','), key, value)
        first = False
    output += ">"
    return output


class LSLConfigParameter(object):
    """
    Class that hold a single configuration parameter.
    """
    
    def __init__(self, name, value, help=None):
        self.name = name
        self.value = value
        self.help = help
        
    def __repr__(self):
        n = self.__class__.__module__+'.'+self.__class__.__name__
        a = [(attr,getattr(self, attr, None)) for attr in ('name', 'value', 'help')]
        return tw_fill(_build_repr(n,a), subsequent_indent='    ')
        
    def __str__(self):
        output = ""
        if self.help is not None:
            output += "# %s\n" % tw_fill(self.help, subsequent_indent='# ')
            
        section, item = self.name.rsplit('.', 1)
        if self.value is None:
            output += "#%s = \n" % item
        else:
            output += "%s = %s\n" % (item, self.value)
        return output


class LSLConfigContainer(object):
    """
    Class for working with all LSL configuration parameters.
    """
    
    def __init__(self, filename=_CONFIG_FILENAME):
        self.filename = filename
        self._changed = False
        
        self._parameters = OrderedDict()
        self._load_config()
        
    def __repr__(self):
        n = self.__class__.__module__+'.'+self.__class__.__name__
        a = [(attr,getattr(self, attr, None)) for attr in ('filename', '_changed', '_parameters')]
        return tw_fill(_build_repr(n,a), subsequent_indent='    ')
        
    def __str__(self):
        output = "# LSL Configuration File\n"
        
        last = None
        for name in self._parameters:
            param = self._parameters[name]
            section, item = name.rsplit('.', 1)
            if section != last:
                output += "\n"
                output += "[ %s ]\n" % section
                last = section
            output += "%s\n" % str(param)
        return output
        
    def __del__(self):
        if self._changed:
            self._save_config()
            
    def _load_config(self):
        """
        Load the configuation from disk.
        """
        
        for section in DEFAULTS_ALL:
            for item in DEFAULTS_ALL[section]:
                name = section+'.'+item
                value = DEFAULTS_ALL[section][item]['value']
                help = DEFAULTS_ALL[section][item]['help']
                param = LSLConfigParameter(name, value, help)
                self._parameters[name] = param
                
        try:
            with open(self.filename, 'r') as fh:
                section = None
                for line in fh:
                    line = line.strip().rstrip()
                    if len(line) < 3:
                        continue
                    if line[0] == '#':
                        continue
                        
                    if line[0] == '[':
                        section = line.replace('[', '').replace(']', '')
                        section = section.strip().rstrip()
                        continue
                    else:
                        name, value = line.split('=', 1)
                        name = name.strip().rstrip()
                        if section is not None:
                            name = section+'.'+name
                        value = value.strip().rstrip()
                        try:
                            value = int(value, 10)
                        except ValueError:
                            try:
                                value = float(value)
                            except ValueError:
                                pass
                        if value == '':
                            continue
                            
                        if name in self._parameters:
                            self._parameters[name].value = value
                        else:
                            param = LSLConfigParameter(name, value)
                            self._parameters[name] = param
                            
        except IOError:
            self._changed = True
            
    def _save_config(self):
        """
        Save the configuation to disk.
        """
        
        with open(self.filename, 'w') as fh:
            fh.write(str(self))
        self._changed = False
        
    def view(self, section):
        """
        Return a configuration sub-container that defaults to looking up 
        values in the specified section.
        """
        
        if section not in DEFAULTS_ALL:
            raise ValueError("Unknown section '%s'" % section)
            
        return LSLConfigSubContainer(self, section)
        
    def get(self, name):
        """
        Return the value of a parameter.
        """
        
        try:
            value = self._parameters[name].value
        except KeyError:
            raise ValueError("Unknown parameter '%s'" % name)
        return value
        
    def set(self, name, value):
        """
        Set the value of a parameter.
        """
        
        try:
            old_value = self._parameters[name].value
            if type(value) != type(old_value):
                raise TypeError("Expected %s but found %s for '%s'" % (type(old_value).__name__,
                                                                       type(value).__name__,
                                                                       name))
            self._parameters[name].value = value
            if value != old_value:
                self._changed = True
                
        except KeyError:
            raise ValueError("Unknown parameter '%s'" % name)
            
    @contextlib.contextmanager
    def set_temp(self, name, value):
        """
        Temporarily set the value of a parameter.  This value will not persist
        across Python sessions.
        """
        
        try:
            old_value = self._parameters[name].value
            if type(value) != type(old_value):
                raise TypeError("Expected %s but found %s for '%s'" % (type(old_value).__name__,
                                                                       type(value).__name__,
                                                                       name))
            self._parameters[name].value = value
            try:
                yield
            finally:
                self._paremeters[name].value = old_value
                
        except KeyError:
            raise ValueError("Unknown parameter '%s'" % name)


class LSLConfigSubContainer(object):
    def __init__(self, container, section):
        self.container = container
        self.section = section
        
    def __repr__(self):
        n = self.__class__.__module__+'.'+self.__class__.__name__
        a = [(attr,getattr(self, attr, None)) for attr in ('container', 'section',)]
        return tw_fill(_build_repr(n,a), subsequent_indent='    ')
        
    def get(self, name):
        """
        Return the value of a parameter.
        """
        
        try:
            value = self.container.get(self.section+'.'+name)
        except ValueError:
            value = self.container.get(name)
        return value
        
    def set(self, name, value):
        """
        Set the value of a parameter.
        """
        
        try:
            self.container.set(self.section+'.'+name, value)
        except KeyError:
            self.container.set(name, value)
        
    def set_temp(self, name, value):
        """
        Temporarily set the value of a parameter.  This value will not persist
        across Python sessions.
        """
        
        try:
            self.container.set_temp(self.section+'.'+name, value)
        except KeyError:
            self.container.set_temp(name, value)


# The LSLConfigContainer that users shoudl use
LSL_CONFIG = LSLConfigContainer()

"""
Module for managing configuration parameters for LSL
"""

import os
import warnings
import contextlib
from datetime import datetime
from collections import OrderedDict
from textwrap import fill as tw_fill

from lsl.version import full_version as lsl_version
from lsl.misc.file_lock import FileLock

_CONFIG_DIR = os.getenv('LSLCONFIGDIR',
                        os.path.join(os.path.expanduser('~'), '.lsl')
                       )

# Create the .lsl directory and set the config filename
try:
    if not os.path.exists(_CONFIG_DIR):
        os.mkdir(_CONFIG_DIR)
    with FileLock(os.path.join(_CONFIG_DIR, 'write.test')):
        with open(os.path.join(_CONFIG_DIR, 'write.test'), 'w') as fh:
            fh.write('test')
        os.unlink(os.path.join(_CONFIG_DIR, 'write.test'))
    _IS_READONLY = False
except OSError:
    _IS_READONLY = True
    warnings.warn('\u001b[33mCannot create or write to on-disk configuration cache\u001b[0m', RuntimeWarning)
    
_CONFIG_FILENAME = os.path.join(_CONFIG_DIR, 'lsl.cfg')


__version__ = "0.3"
__all__ = ['LSL_CONFIG',]


# Default values
## lsl.common.sdf/sdfADP/idf
DEFAULTS_OBS = OrderedDict()
DEFAULTS_OBS['observer_name'] = {'value': None,
                                 'help':  'Observer name for auto-filling Observer classes',
                                 'is_optional': True}
DEFAULTS_OBS['observer_id'] = {'value': None,
                               'help':  'Observer ID number for auto-filling Observer classes',
                               'is_optional': True}
DEFAULTS_OBS['project_name'] = {'value': None,
                                'help':  'Project name for auto-filling Project classes',
                                'is_optional': True}
DEFAULTS_OBS['project_id'] = {'value': None,
                              'help':  'Project ID for auto-filling Project classes',
                              'is_optional': True}

## lsl.reader.ldp
DEFAULTS_LDP = OrderedDict()
DEFAULTS_LDP['tbn_buffer_size'] = {'value': 20,
                                   'help':  'TBN ring buffer size in timestamps',
                                   'is_optional': False}
DEFAULTS_LDP['drx_buffer_size'] = {'value': 20,
                                   'help':  'DRX ring buffer size in timestamps',
                                   'is_optional': False}
DEFAULTS_LDP['drx_autofill_size'] = {'value': 50,
                                     'help':  'maximum DRX gap in timestamps that can be auto-filled with zeros without throwing a timetag error/warning',
                                     'is_optional': False}
DEFAULTS_LDP['tbf_buffer_size'] = {'value': 25,
                                   'help':  'TBF ring buffer size in timestamps',
                                   'is_optional': False}
DEFAULTS_LDP['cor_buffer_size'] = {'value': 5,
                                   'help':  'COR ring buffer size in timestamps',
                                   'is_optional': False}

## lsl.astro
DEFAULTS_ASTRO = OrderedDict()
DEFAULTS_ASTRO['leapsec_url'] = {'value': 'https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat',
                                 'help':'URL for accessing leap second information',
                                 'is_optional': False}

## lsl.misc.ionosphere
DEFAULTS_IONO = OrderedDict()
DEFAULTS_IONO['igs_url'] = {'value': 'ftps://gdc.cddis.eosdis.nasa.gov/gps/products/ionex/',
                            'help':  'primary URL for accessing the IGS data products',
                            'is_optional': False}
DEFAULTS_IONO['igs_mirror'] = {'value': 'ftp://gssc.esa.int/gnss/products/ionex/',
                               'help':  'mirror URL for accessing the IGS data products',
                               'is_optional': True}
DEFAULTS_IONO['jpl_url'] = {'value': 'ftps://gdc.cddis.eosdis.nasa.gov/gps/products/ionex/',
                            'help':  'primary URL for accessing the JPL data products',
                            'is_optional': False}
DEFAULTS_IONO['jpl_mirror'] = {'value': 'ftp://gssc.esa.int/gnss/products/ionex/',
                               'help':  'mirror URL for accessing the JPL data products',
                               'is_optional': True}
DEFAULTS_IONO['emr_url'] = {'value': 'ftps://gdc.cddis.eosdis.nasa.gov/gps/products/ionex/',
                            'help':  'primary URL for accessing the EMR data products',
                            'is_optional': False}
DEFAULTS_IONO['emr_mirror'] = {'value': 'ftp://gssc.esa.int/gnss/products/ionex/',
                               'help':  'mirror URL for accessing the EMR data products',
                               'is_optional': True}
DEFAULTS_IONO['uqr_url'] = {'value': 'ftps://gdc.cddis.eosdis.nasa.gov/gps/products/ionex/',
                            'help':  'primary URL for accessing the UQR data products',
                            'is_optional': False}
DEFAULTS_IONO['uqr_mirror'] = {'value': 'ftp://gssc.esa.int/gnss/products/ionex/',
                               'help':  'mirror URL for accessing the UQR data products',
                               'is_optional': True}
DEFAULTS_IONO['code_url'] = {'value': 'ftps://gdc.cddis.eosdis.nasa.gov/gps/products/ionex/',
                             'help':  'primary URL for accessing the CODE data products',
                             'is_optional': False}
DEFAULTS_IONO['code_mirror'] = {'value': 'ftp://gssc.esa.int/gnss/products/ionex/',
                                'help':  'mirror URL for accessing the CODE data products',
                                'is_optional': True}
DEFAULTS_IONO['ustec_url'] = {'value': 'http://www.ngdc.noaa.gov/stp/iono/ustec/products/',
                              'help':  'primary URL for accessing the USTEC data products',
                              'is_optional': False}
DEFAULTS_IONO['ustec_mirror'] = {'value': None,
                                 'help':  'mirror URL for accessing the USTEC data products',
                                 'is_optional': True}
DEFAULTS_IONO['glotec_url'] = {'value': 'https://services.swpc.noaa.gov/experimental/products/glotec/geojson_2d_urt/',
                              'help':  'primary URL for accessing the GloTEC data products',
                              'is_optional': False}
DEFAULTS_IONO['glotec_mirror'] = {'value': None,
                                 'help':  'mirror URL for accessing the GloTEC data products',
                                 'is_optional': True}
DEFAULTS_IONO['max_cache_size'] = {'value': -1,
                                   'help':  'maximum cache size in MB; <= 0 disables cache size limiting',
                                   'is_optional': False}

## lsl.misc.telemetry
DEFAULTS_TELEMETRY = OrderedDict()
DEFAULTS_TELEMETRY['enabled'] = {'value': True,
                                 'help': 'whether or not LSL telemetry reporting is enabled',
                                 'is_optional': False}
DEFAULTS_TELEMETRY['max_entries'] = {'value': 100,
                                     'help':  'maximum number of entries to accumlate before reporting',
                                     'is_optional': False}
DEFAULTS_TELEMETRY['timeout'] = {'value': 30,
                                 'help':  'upload timeout in seconds',
                                 'is_optional': False}

## Download parameters
DEFAULTS_DOWNLOAD = OrderedDict()
DEFAULTS_DOWNLOAD['block_size'] = {'value': 8192,
                                   'help':  'download block size in bytes',
                                   'is_optional': False}
DEFAULTS_DOWNLOAD['timeout'] = {'value': 120,
                                'help':  'download timeout in seconds',
                                'is_optional': False}
DEFAULTS_DOWNLOAD['refresh_age'] = {'value': 14,
                                    'help':  'data cache refresh age in days',
                                    'is_optional': False}

## Everything
DEFAULTS_ALL = OrderedDict()
DEFAULTS_ALL['observing'] = DEFAULTS_OBS
DEFAULTS_ALL['ldp'] = DEFAULTS_LDP
DEFAULTS_ALL['astro'] = DEFAULTS_ASTRO
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
    
    def __init__(self, name, value, help=None, is_optional=False):
        self.name = name
        self.value = value
        self.help = help
        self.is_optional = is_optional
        
    def __repr__(self):
        n = self.__class__.__module__+'.'+self.__class__.__name__
        a = [(attr,getattr(self, attr, None)) for attr in ('name', 'value', 'help', 'is_optional')]
        return tw_fill(_build_repr(n,a), subsequent_indent='    ')
        
    def __str__(self):
        output = ""
        if self.help is not None:
            output += "# %s\n" % tw_fill(self.help, subsequent_indent='# ')
            
        section, item = self.name.rsplit('.', 1)
        if self.value is None:
            output += f"#{item} = \n"
        else:
            output += f"{item} = {self.value}\n"
        return output


class LSLConfigContainer(object):
    """
    Class for working with all LSL configuration parameters.
    """
    
    def __init__(self, filename=_CONFIG_FILENAME):
        self.filename = filename
        self._loaded = 0.0
        self._version_changed = False
        self._param_changed = False
        self._param_changed_list = {}

        self._parameters = OrderedDict()
        self._load_config()
        
    @property
    def dirname(self):
        """
        Directory where the LSL configuration information is stored.
        """
        
        return os.path.dirname(self.filename)
        
    def __repr__(self):
        n = self.__class__.__module__+'.'+self.__class__.__name__
        a = [(attr,getattr(self, attr, None)) for attr in ('filename', '_param_changed', '_version_changed', '_parameters')]
        return tw_fill(_build_repr(n,a), subsequent_indent='    ')
        
    def __str__(self):
        output = "# LSL Configuration File\n"
        output += f"# LSL Version: {lsl_version}\n"
        
        last = None
        for name in self._parameters:
            param = self._parameters[name]
            section, item = name.rsplit('.', 1)
            if section != last:
                output += "\n"
                output += f"[ {section} ]\n"
                last = section
            output += f"{str(param)}\n"
        return output
        
    def _load_config(self):
        """
        Load the configuation from disk.
        """
        
        for section in DEFAULTS_ALL:
            for item in DEFAULTS_ALL[section]:
                name = section+'.'+item
                value = DEFAULTS_ALL[section][item]['value']
                help = DEFAULTS_ALL[section][item]['help']
                is_optional = DEFAULTS_ALL[section][item]['is_optional']
                param = LSLConfigParameter(name, value, help, is_optional)
                self._parameters[name] = param
                
        try:
            with open(self.filename, 'r') as fh:
                self._loaded = os.path.getmtime(self.filename)
                
                section = None
                for line in fh:
                    line = line.strip().rstrip()
                    if len(line) < 3:
                        continue
                    if line[0] == '#':
                        if line.find('LSL Version:') != -1:
                            existing_version = line.split('LSL Version:', 1)[1]
                            existing_version = existing_version.strip().rstrip()
                            if existing_version != lsl_version:
                                self._version_changed = True
                                self._param_changed = True
                                
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
                        if value in ('True', 'False'):
                            value = True if value == 'True' else False
                        else:
                            try:
                                value = int(value, 10)
                            except ValueError:
                                try:
                                    value = float(value)
                                except ValueError:
                                    pass
                            if value == '':
                                value = None
                                
                        if name in self._parameters:
                            self._parameters[name].value = value
                        else:
                            param = LSLConfigParameter(name, value)
                            self._parameters[name] = param
                            
        except IOError:
            # File doesn't exist, will need to create it
            self._version_changed = True
            self._param_changed = True
            
    def _save_config(self):
        """
        Save the configuation to disk.
        """
        
        if not _IS_READONLY:
            if os.path.exists(self.filename) and len(self._param_changed_list) > 0:
                mtime = os.path.getmtime(self.filename)
                if mtime > self._loaded:
                    # The file was modified on disk after we loaded it
                    # Abandon our parameter changes but still save if version changed
                    warnings.warn('\u001b[33mConfiguration file changed on disk, abandoning changes from this session\u001b[0m', RuntimeWarning)
                    for key,value in self._param_changed_list.items():
                        warnings.warn("\u001b[33m  %s: %s -> %s\u001b[0m" % (key, value[0], value[1]), RuntimeWarning)
                    self._param_changed_list.clear()
                    # Don't reset _param_changed if version changed - we still need to update the file
                    if not self._version_changed:
                        self._param_changed = False
                        
            if self._param_changed:
                with FileLock(self.filename) as lock:
                    assert(lock.locked())
                    
                    with open(self.filename, 'w') as fh:
                        fh.write(str(self))
                    self._version_changed = False
                    self._param_changed = False
                    self._param_changed_list.clear()
                    
    def view(self, section):
        """
        Return a configuration sub-container that defaults to looking up 
        values in the specified section.
        """
        
        if section not in DEFAULTS_ALL:
            raise ValueError(f"Unknown section '{section}'")
            
        return LSLConfigSubContainer(self, section)
        
    def get(self, name):
        """
        Return the value of a parameter.
        """
        
        try:
            value = self._parameters[name].value
        except KeyError:
            raise ValueError(f"Unknown parameter '{name}'")
        return value
        
    def set(self, name, value):
        """
        Set the value of a parameter.
        """
        
        try:
            old_value = self._parameters[name].value
            is_optional = self._parameters[name].is_optional
            
            # Allow type changes only for optional parameters (those with None defaults)
            # For non-optional parameters, enforce strict type checking
            if type(value) != type(old_value) and not is_optional:
                raise TypeError("Expected %s but found %s for '%s'" % (type(old_value).__name__,
                                                                       type(value).__name__,
                                                                       name))
            self._parameters[name].value = value
            if value != old_value:
                if name in self._param_changed_list:
                    # Get the original value from when we first changed this parameter
                    original_value = self._param_changed_list[name][0]
                else:
                    # This is the first change to this parameter
                    original_value = old_value
                    
                # Check if we're setting back to the original value
                if value == original_value:
                    # Remove from changed list since we're back to the original value
                    if name in self._param_changed_list:
                        del self._param_changed_list[name]
                    # Update _param_changed flag based on whether we still have changes
                    self._param_changed = (len(self._param_changed_list) > 0) or self._version_changed
                else:
                    # Track the change from original to new value
                    self._param_changed_list[name] = (original_value, value)
                    self._param_changed = True
                    
        except KeyError:
            raise ValueError(f"Unknown parameter '{name}'")
            
    @contextlib.contextmanager
    def set_temp(self, name, value):
        """
        Temporarily set the value of a parameter.  This value will not persist
        across Python sessions.
        """
        
        try:
            old_value = self._parameters[name].value
            is_optional = self._parameters[name].is_optional
            
            # Allow type changes only for optional parameters (those with None defaults)
            # For non-optional parameters, enforce strict type checking
            if type(value) != type(old_value) and not is_optional:
                raise TypeError("Expected %s but found %s for '%s'" % (type(old_value).__name__,
                                                                       type(value).__name__,
                                                                       name))
            self._parameters[name].value = value
            try:
                yield
            finally:
                self._parameters[name].value = old_value
                
        except KeyError:
            raise ValueError(f"Unknown parameter '{name}'")


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
            return self.container.set_temp(self.section+'.'+name, value)
        except KeyError:
            return self.container.set_temp(name, value)


#: The LSLConfigContainer that users should use
LSL_CONFIG = LSLConfigContainer()

import atexit
atexit.register(LSL_CONFIG._save_config)

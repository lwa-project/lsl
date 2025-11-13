"""
Common utilities shared by the various metabundle modules.
"""

import os
import re
import shutil
import tarfile
import tempfile

from functools import lru_cache

from lsl.common.mcs import MIB


__version__ = '0.1'
__all__ = ['FILENAME_RE', '_open_tarball', '_get_members', 'managed_mkdtemp',
           'get_beamformer_min_delay', 'get_session_metadata',
           'get_asp_configuration', 'get_asp_configuration_summary']


# Regular expression for figuring out filenames
FILENAME_RE = re.compile(r'(?P<projectID>[a-zA-Z0-9]{1,8})_(?P<sessionID>\d+)(_(?P<obsID>\d+)(_(?P<obsOutcome>\d+))?)?.*\..*')


@lru_cache(maxsize=1)
def _open_tarball(tarname):
    return tarfile.open(tarname, mode='r:*')


@lru_cache(maxsize=1)
def _get_members(tarname):
    tf = _open_tarball(tarname)
    return tf.getmembers()


class managed_mkdtemp(object):
    """
    Wrapper class around tempfile.mkdtemp to enable 'with' statements with 
    automatic cleanup.
    """
    
    def __init__(self, suffix='', prefix='tmp', dir=None):
        self._dir = tempfile.mkdtemp(suffix, prefix, dir)
        
    def __enter__(self):
        return self._dir
        
    def __exit__(self, type, value, tb):
        shutil.rmtree(self._dir, ignore_errors=True)
        
    @property
    def name(self):
        return self._dir


def get_beamformer_min_delay(tarname):
    """
    Given an MCS meta-data tarball, extract the minimum beamformer delay in 
    samples and return it.  If no minimum delay can be found in the tarball, 
    None is returned.
    """
    
    with managed_mkdtemp(prefix='metadata-bundle-') as tempDir:
        # Extract the mindelay.txt file.  If mindelay.txt cannot be found, None
        # is returned via the try...except block.
        tf = _open_tarball(tarname)
        try:
            ti = tf.getmember('mindelay.txt')
        except KeyError:
            return None
        tf.extractall(path=tempDir, members=[ti,])
        
        # Parse the SDM file and build the SDM instance
        with open(os.path.join(tempDir, 'mindelay.txt'), 'r') as fh:
            try:
                mindelay = int(fh.read(), 10)
            except ValueError:
                mindelay = None
                
    return mindelay


def get_session_metadata(tarname):
    """
    Given an MCS meta-data tarball, extract the session meta-data file (MCS0030, 
    Section 7) and return a dictionary of observations that contain dictionaries 
    of the OP_TAG (tag), DRSU Barcode (drsu), OBS_OUTCOME (outcome), and the 
    MSG (msg).
    
    .. versionchanged:: 0.6.5
        Update to the new _metadata.txt format
    """
    
    with managed_mkdtemp(prefix='metadata-bundle-') as tempDir:
        path, basename = os.path.split(tarname)
        basename, ext = os.path.splitext(basename)
        
        # Extract the session meta-data file (_metadata.txt)
        tf = _open_tarball(tarname)
        try:
            ti = tf.getmember('%s_metadata.txt' % basename)
        except KeyError:
            for ti in _get_members(tarname):
                if ti.name[-13:] == '_metadata.txt':
                    break
        tf.extractall(path=tempDir, members=[ti,])
        
        # Read in the SMF
        filename = os.path.join(tempDir, ti.name)
        with open(filename, 'r') as fh:
            # Define a regular expresion to match the latest format
            lineRE = re.compile(r"\s*(?P<id>\d{1,}?)\s+\[(?P<tag>[\d_]+?)\]\s+\['?(?P<barcode>.+?)'?\]\s+(?P<outcome>\d)\s+\[(?P<msg>.*?)\]")
            
            result = {}
            for line in fh:
                line = line.replace('\n', '')
                if len(line) == 0:
                    continue
                    
                mtch = lineRE.search(line)
                if mtch is not None:
                    ## If it matches the new format
                    obsID = mtch.group('id')
                    opTag = mtch.group('tag')
                    drsuBarcode = mtch.group('barcode')
                    if drsuBarcode[:3] == 'Err':
                        try:
                            drsuBarcode = result[int(obsID)-1]['barcode']
                        except KeyError:
                            drsuBarcode = 'UNK'
                    obsOutcome = mtch.group('outcome')
                    msg = mtch.group('msg')
                    
                else:
                    ## Otherwise, I don't really know how the messages will look so we use this try...except
                    ## block should take care of the various situations.
                    try:
                        obsID, opTag, drsuBarcode, obsOutcome, msg = line.split(None, 4)
                        opTag = opTag.replace('[', '')
                        opTag = opTag.replace(']', '')
                        drsuBarcode = drsuBarcode.replace('[', '')
                        drsuBarcode = drsuBarcode.replace(']', '')
                        drsuBarcode = drsuBarcode.replace("'", '')
                    except ValueError:
                        try:
                            obsID, opTag, drsuBarcode, obsOutcome = line.split(None, 3)
                            msg = 'UNK'
                        except ValueError:
                            try:
                                obsID, opTag, obsOutcome = line.split(None, 2)
                                drsuBarcode = 'UNK'
                                obsOutcome = '-1'
                                msg = 'UNK'
                            except ValueError:
                                obsID, obsOutcome = line.split(None, 1)
                                drsuBarcode = 'UNK'
                                opTag = 'UNK'
                                msg = 'UNK'
                                
                obsID = int(obsID)
                obsOutcome = int(obsOutcome) if obsOutcome != 'Failed' else 1
                result[obsID] = {'tag': opTag, 'barcode': drsuBarcode, 
                                 'outcome': obsOutcome, 'msg': msg}
                
    # Return
    return result


def get_asp_configuration(tarname, which='beginning'):
    """
    Given an MCS meta-data tarball, extract the ASP MIB contained in it and return 
    a dictionary of values for the filter, AT1, AT2, and AT3.  The 'which'
    keyword is used to specify whether or not the configuration returned is at the
    beginning (default) or end of the session.
    
    .. versionadded:: 0.6.5
    """
    
    which = which.lower()
    if which not in ('beginning', 'begin', 'end'):
        raise ValueError(f"Unknown configuration time '{which}'")
        
    # Stub ASP configuration
    aspConfig = {'asp_filter':      [-1 for i in range(256)],
                 'asp_atten_1':     [-1 for i in range(256)],
                 'asp_atten_2':     [-1 for i in range(256)],
                 'asp_atten_3':     [-1 for i in range(256)]}
    
    with managed_mkdtemp(prefix='metadata-bundle-') as tempDir:
        # Find the .pag file and extract it
        tf = _open_tarball(tarname)
        mibs = []
        for ti in _get_members(tarname):
            if ti.name.find('_ASP_%s' % which[:5]) != -1:
                if ti.name[-4:] in ('.pag', '.gdb'):
                    mibs.append(ti)
                    
        if len(mibs) > 0:
            tf.extractall(path=tempDir, members=mibs)
            
            # Read in the right MIB
            aspMIB = {}
            for mib in mibs:
                if which[:5] == 'begin' and mib.name.find('_ASP_begin') == -1:
                    continue
                if which == 'end' and mib.name.find('_ASP_end') == -1:
                    continue
                    
                aspMIB = MIB.from_file(os.path.join(tempDir, mib.name))
                break
                
            # Extract the configuration
            for key in aspMIB.keys():
                values = [int(v) for v in key.split('.')]
                
                if values[0] == 3:
                    ## Filter
                    aspConfig['asp_filter'][values[1]-1] = int(aspMIB[key].value)
                    
                elif values[0] == 4:
                    ## Attenuators
                    if values[1] == 1:
                        ### AT1
                        aspConfig['asp_atten_1'][values[2]-1] = int(aspMIB[key].value)
                    elif values[1] == 2:
                        ### AT2
                        aspConfig['asp_atten_2'][values[2]-1] = int(aspMIB[key].value)
                    elif values[1] == 3:
                        ### AT3
                        aspConfig['asp_atten_3'][values[2]-1] = int(aspMIB[key].value)
                    else:
                        pass
                        
                else:
                    pass
                    
    return aspConfig


def get_asp_configuration_summary(tarname, which='beginning'):
    """
    Similar to get_asp_configuration, but returns only a single value for each
    of the four ASP paramters:  filter, AT, AT2, and AT3.  The values
    are based off the mode of the parameter.
    
    .. versionadded:: 0.6.5
    """
    
    # Get the full configuration
    aspConfig = get_asp_configuration(tarname, which=which)
    
    # Count
    count = {}
    for param in aspConfig.keys():
        count[param] = {}
        for ant in range(len(aspConfig[param])):
            value = aspConfig[param][ant]
            try:
                count[param][value] += 1
            except KeyError:
                count[param][value] = 1
                
    # Modes
    mode = {}
    for param in count.keys():
        best = 0
        mode[param] = 0
        
        for value in count[param].keys():
            num = count[param][value]
            if num > best:
                best = num
                mode[param] = value
                
    # Done
    return mode

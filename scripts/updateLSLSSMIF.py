#!/usr/bin/env python3

"""
Utility for updating/changing the LWA station SSMIF files in between LSL 
releases.
"""

import os
import re
import sys
import time
from urllib.request import urlopen
import hashlib
import argparse
import calendar
from datetime import datetime

from lsl.common.data_access import DataAccess
from lsl.common.progress import DownloadBar

from lsl.config import LSL_CONFIG

from lsl.misc import telemetry
telemetry.track_script()


# Regular expression for finding the SSMIF version from the file contents
versionRE = re.compile(r'(?P<year>\d{4}) (?P<month>[a-zA-Z]{3,4}) (?P<day>\d{1,2}) by (?P<author>.*)')


def _parse_index(index):
    """
    Parse the archive listing of SSMIF version and return a list of 
    filename/date tuples.
    """
    
    from xml.etree import ElementTree as ET
    from bs4 import BeautifulSoup
    
    # Find the table
    start = index.find('<table>')
    stop  = index.find('</table>')
    index = index[start:stop+8]
    
    # Clean it up in such a way that ElementTree can parse it
    myMassage = [(re.compile('<!([^--].*)>'), lambda match: '<!--' + match.group(1) + '-->'), 
                 (re.compile('<hr>'), lambda match: ''), 
                 (re.compile('&nbsp;'), lambda match: ' '), 
                 (re.compile('<a.*?>(.*)</a>'), lambda mtch: mtch.group(1))]
    for massage in myMassage:
        regex, replace = massage
        index = re.sub(regex, replace, index)
        
    soup = BeautifulSoup(index)
    index = soup.prettify()
    index = index.replace('<html>', '<?xml version="1.0" encoding="utf-8"?>')
    for tag in ('body', 'html'):
        index = index.replace(f"<{tag}>", '')
        index = index.replace(f"</{tag}>", '')
        
    # Parse it
    table = ET.XML(index)
    rows = iter(table)
    
    # Extract the SSMIF entries
    versions = []
    for row in rows:
        values = [col.text for col in row]
        if len(values) != 5:
            continue
            
        filename = values[1].lstrip().rstrip()
        if filename[:5] != 'SSMIF':
            continue
        if filename.find('CURRENT') != -1:
            continue
            
        date = filename.split('_', 1)[1]
        date = date.split('.', 1)[0]
        versions.append( (filename, date) )
        
    # Done
    return versions


def _compute_md5(fh, block_size=262144):
    """
    Compute the MD5 checksum of an open file handle.
    """
    
    m = hashlib.md5()
    while True:
        block = fh.read(block_size)
        if len(block) == 0:
            break
        m.update(block)
        
    return m.hexdigest()


def main(args):
    # Current LSL SSMIF
    if args.lwasv:
        _name = 'LWA-SV'
        _ssmif = 'lwasv-ssmif.txt'
        _url = "https://lda10g.alliance.unm.edu/metadata/lwasv/ssmif/"
    elif args.lwana:
        _name = 'LWA-NA'
        _ssmif = 'lwana-ssmif.txt'
        _url = "https://lda10g.alliance.unm.edu/metadata/lwana/ssmif/"
    else:
        _name = 'LWA1'
        _ssmif = 'lwa1-ssmif.txt'
        _url = "https://lda10g.alliance.unm.edu/metadata/lwa1/ssmif/"
        
    urlToDownload = None
    if args.revert:
        # Revert and Upgrade URL
        
        try:
            ## Retrieve the list
            try:
                with urlopen(_url, timeout=LSL_CONFIG.get('download.timeout')) as ah:
                    index = ah.read()
                    try:
                        index = index.decode()
                    except AttributeError:
                        pass
            except Exception as e:
                print(f"Error:  Cannot download SSMIF listing, {str(e)}")
                
            ## Parse
            versions = _parse_index(index)
            
            ## Prompt the user for the version to revert to
            for i,(filename,date) in enumerate(versions):
                print(f"{i}: {filename}")
            i = -1
            while i not in range(0, len(versions)):
                i = input("Enter SSMIF to revert to: ")
                try:
                    i = int(i)
                except ValueError:
                    print("-> Invalid value")
                    i = -1
            print(" ")
            
            ## Build the URL
            urlToDownload = f"{_url}/{versions[i][0]}"
        except Exception as e:
            print(f"Error:  Cannot process reversion, {str(e)}")
            
    elif args.update:
        # Update to the latest version
        
        urlToDownload = f"{_url}/SSMIF_CURRENT.txt"
        
    elif args.file is not None:
        # Use the specified file
        
        urlToDownload = os.path.abspath(args.file)
        
    # Put the new SSMIF in place
    if urlToDownload is not None:
        ## Retrieve
        try:
            print(f"Downloading {urlToDownload}")
            mtime = 0.0
            remote_size = 1
            with urlopen(urlToDownload, timeout=LSL_CONFIG.get('download.timeout')) as ah:
                remote_size = int(ah.headers["Content-Length"])
                mtime = ah.headers["Last-Modified"]
                
                mtime = datetime.strptime(mtime, "%a, %d %b %Y %H:%M:%S GMT")
                mtime = calendar.timegm(mtime.timetuple())
                
                pbar = DownloadBar(max=remote_size)
                while True:
                    new_data = ah.read(LSL_CONFIG.get('download.block_size'))
                    if len(new_data) == 0:
                        break
                    pbar.inc(len(new_data))
                    try:
                        newSSMIF += new_data
                    except NameError:
                        newSSMIF = new_data
                    sys.stdout.write(pbar.show()+'\r')
                    sys.stdout.flush()
                sys.stdout.write(pbar.show()+'\n')
                sys.stdout.flush()
        except Exception as e:
            raise RuntimeError(f"Cannot download SSMIF: {str(e)}")
            
        ## Save
        try:
            with DataAccess.open(_ssmif, 'wb') as fh:
                fh.write(newSSMIF)
            with DataAccess.open(_ssmif+"_meta", 'w') as fh:
                fh.write(f"created: {time.time():.0f}\n")
                fh.write(f"source: {urlToDownload}\n")
                fh.write(f"source size: {remote_size} B\n")
                fh.write(f"source last modified: {mtime:.0f}\n")
        except Exception as e:
            raise RuntimeError("Cannot %s SSMIF: %s" % ('update' if args.update else 'revert', str(e)))
            
    # Summarize the SSMIF
    ## Filesystem information
    size = DataAccess.getsize(_ssmif)
    mtime = datetime.utcfromtimestamp(DataAccess.stat(_ssmif)[8])
    age = datetime.utcnow() - mtime
    
    ## MD5 checksum
    with DataAccess.open(_ssmif, 'rb') as fh:
        md5 = _compute_md5(fh)
        
    ## SSMIF version (date)
    with DataAccess.open(_ssmif, 'r') as fh:
        lines = [fh.readline() for i in range(10)]
        
    version = None
    for line in lines:
        mtch = versionRE.search(line)
        if mtch is not None:
            try:
                version = datetime.strptime("%s %s %s" % (mtch.group('year'), mtch.group('month'), mtch.group('day')), "%Y %b %d")
            except:
                version = datetime.strptime("%s %s %s" % (mtch.group('year'), mtch.group('month'), mtch.group('day')), "%Y %B %d")
            break
            
    print("LSL %s SSMIF%s:" % (_name.upper(), ' (updated)' if args.update else '',))
    print(f"  Size: {size} bytes")
    print("  SSMIF Version: %s" % version.strftime("%Y %b %d"))
    print("  File Last Modified: %s (%i day%s ago)" % (mtime.strftime("%Y-%m-%d %H:%M:%S UTC"), age.days, 's' if age.days != 1 else ''))
    print(f"  MD5 Sum: {md5}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='update the internal LWA station SSMIFs used by LSL', 
            epilog='The -s/--lwasv option changes the behavior of all other options, e.g., -u/--update updates the LWA-SV SSMIF.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    sgroup = parser.add_mutually_exclusive_group(required=False)
    sgroup.add_argument('-s', '--lwasv', action='store_true', 
                        help='update LWA-SV instead of LWA1')
    sgroup.add_argument('-n', '--lwana', action='store_true', 
                        help='update LWA-NA instead of LWA1')
    parser.add_argument('-u', '--update', action='store_true', 
                        help='update the default LWA1 SSMIF')
    parser.add_argument('-r', '--revert', action='store_true', 
                        help='reveret the default LWA1 SSMIF to an older version')
    parser.add_argument('-f', '--file', type=str, 
                        help='update the default LWA1 SSMIF using the specified file')
    args = parser.parse_args()
    main(args)

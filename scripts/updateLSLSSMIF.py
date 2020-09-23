#!/usr/bin/env python

"""
Utility for updating/changing the LWA station SSMIF files in between LSL 
releases.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import re
import sys
try:
    from urllib2 import urlopen
except ImportError:
    from urllib.request import urlopen
import hashlib
import argparse
from datetime import datetime

from lsl.common.paths import DATA as dataPath
from lsl.common.progress import ProgressBarPlus

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
    from BeautifulSoup import BeautifulSoup
    
    # Find the table
    start = index.find('<table>')
    stop  = index.find('</table>')
    index = index[start:stop+8]
    
    # Clean it up in such a way that ElementTree can parse it
    myMassage = [(re.compile('<!([^--].*)>'), lambda match: '<!--' + match.group(1) + '-->'), 
            (re.compile('<hr>'), lambda match: ''), 
            (re.compile('&nbsp;'), lambda match: ' '), 
            (re.compile('<a.*?>(.*)</a>'), lambda mtch: mtch.group(1))]
    soup = BeautifulSoup(index, markupMassage=myMassage)
    index = soup.prettify()
    
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


def _compute_md5(filename, block_size=262144):
    """
    Compute the MD5 checksum of a file.
    """
    
    fh = open(filename, 'r')
    
    m = hashlib.md5()
    while True:
        block = fh.read(block_size)
        if len(block) == 0:
            break
        m.update(block)
        
    fh.close()
    
    return m.hexdigest()


def main(args):
    # Current LSL SSMIF
    if args.lwasv:
        _name = 'LWA-SV'
        _ssmif = os.path.join(dataPath, 'lwasv-ssmif.txt')
        _url = "https://lda10g.alliance.unm.edu/metadata/lwasv/ssmif/"
    else:
        _name = 'LWA1'
        _ssmif = os.path.join(dataPath, 'lwa1-ssmif.txt')
        _url = "https://lda10g.alliance.unm.edu/metadata/lwa1/ssmif/"
        
    urlToDownload = None
    if args.revert:
        # Revert and Upgrade URL
        
        try:
            ## Retrieve the list
            ah = urlopen(_url)
            index = ah.read()
            try:
                index = index.decode(encoding='ascii', errors='ignore')
            except AttributeError:
                pass
            ah.close()
            
            ## Parse
            versions = _parse_index(index)
            
            ## Prompt the user for the version to revert to
            for i,(filename,date) in enumerate(versions):
                print("%i: %s" % (i, filename))
            i = -1
            while i not in range(0, len(versions)):
                i = raw_input("Enter SSMIF to revert to: ")
                try:
                    i = int(i)
                except ValueError:
                    print("-> Invalid value")
                    i = -1
            print(" ")
            
            ## Build the URL
            urlToDownload = "%s/%s" % (_url, versions[i][0])
        except Exception as e:
            print("Error:  Cannot process reversion, %s" % str(e))
            
    elif args.update:
        # Update to the latest version
        
        urlToDownload = "%s/SSMIF_CURRENT.txt" % _url
        
    elif args.file is not None:
        # Use the specified file
        
        urlToDownload = os.path.abspath(args.file)
        
    # Put the new SSMIF in place
    if urlToDownload is not None:
        ## Retrieve
        try:
            print("Downloading %s" % urlToDownload)
            ah = urlopen(urlToDownload)
            meta = ah.info()
            pbar = ProgressBarPlus(max=int(meta.getheaders("Content-Length")[0]))
            while True:
                new_data = ah.read(32768)
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
            ah.close()
        except Exception as e:
            print("Error:  Cannot download SSMIF, %s" % str(e))
            
        ## Save
        try:
            fh = open(_ssmif, 'wb')
            fh.write(newSSMIF)
            fh.close()
        except Exception as e:
            print("Error:  Cannot %s SSMIF, %s" % ('update' if args.update else 'revert', str(e)))
            
    # Summarize the SSMIF
    ## Filesystem information
    size = os.path.getsize(_ssmif)
    mtime = datetime.utcfromtimestamp(os.stat(_ssmif)[8])
    age = datetime.utcnow() - mtime
    
    ## MD5 checksum
    md5 = _compute_md5(_ssmif)
    
    ## SSMIF version (date)
    fh = open(_ssmif, 'r')
    lines = [fh.readline() for i in range(10)]
    fh.close()
    
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
    print("  Size: %i bytes" % size)
    print("  SSMIF Version: %s" % version.strftime("%Y %b %d"))
    print("  File Last Modified: %s (%i day%s ago)" % (mtime.strftime("%Y-%m-%d %H:%M:%S UTC"), age.days, 's' if age.days != 1 else ''))
    print("  MD5 Sum: %s" % md5)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='update the internal LWA station SSMIFs used by LSL', 
            epilog='The -s/--lwasv option changes the behavior of all other options, e.g., -u/--update updates the LWA-SV SSMIF.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    parser.add_argument('-s', '--lwasv', action='store_true', 
                        help='update LWA-SV instead of LWA1')
    parser.add_argument('-u', '--update', action='store_true', 
                        help='update the default LWA1 SSMIF')
    parser.add_argument('-r', '--revert', action='store_true', 
                        help='reveret the default LWA1 SSMIF to an older version')
    parser.add_argument('-f', '--file', type=str, 
                        help='update the default LWA1 SSMIF using the specified file')
    args = parser.parse_args()
    main(args)

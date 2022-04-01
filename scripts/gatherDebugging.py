#!/usr/bin/env python

"""
Script to gather information about the Python interpreter, modules, 
C libraries, numpy installation, and LSL installation to help with 
debugging and install issues.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import re
import sys
import argparse
import platform
import subprocess

from lsl.misc import telemetry
telemetry.track_script()


def main(args):
    #
    # Clean the path
    #
    sys.path = sys.path[1:]

    #
    # Python interpreter
    #
    print("Executable path: %s" % sys.executable)
    print("Platform: %s" % sys.platform)
    print("Version: %s" % sys.version)
    print("API: %s" % sys.api_version)
    print("Bits: %s\nLinkage: %s" % platform.architecture())
    print(" ")
    
    #
    # Python Module Check
    #
    
    ## Required
    for mod in ('numpy', 'scipy', 'astropy', 'ephem', 'aipy', 'pytz'):
        try:
            exec "import %s" % mod
        except ImportError as e:
            if (str(e)).find('not found') != -1:
                print( "%s: not found" % mod)
            else:
                print("%s: WARNING import error '%s'" % (mod, str(e)))
        else:
            try:
                version = eval("%s.version.version" % mod)
            except AttributeError:
                try:
                    version = eval("%s.__version__" % mod)
                except AttributeError:
                    try:	
                        versionRE = re.compile(r'%s-(?P<version>[\d\.]+)-py.*' % mod)
                        mtch = versionRE.search(eval("%s.__file__" % mod))
                        version = mtch.group('version')
                    except:
                        version = "unknown"
            print("%s:  version %s" % (mod, version))
            
    ## Optional
    for mod in ('matplotlib', 'h5py', 'psrfits_utils'):
        try:
            exec "import %s" % mod
        except ImportError as e:
            if (str(e)).find('not found') != -1:
                print( "%s: not found" % mod)
            #else:
                #print("%s: WARNING import error '%s'" % (mod, str(e)))
        else:
            try:
                version = eval("%s.version.version" % mod)
            except AttributeError:
                try:
                    version = eval("%s.__version__" % mod)
                except AttributeError:
                    try:	
                        versionRE = re.compile(r'%s-(?P<version>[\d\.]+)-py.*' % mod)
                        mtch = versionRE.search(eval("%s.__file__" % mod))
                        version = mtch.group('version')
                    except:
                        version = "unknown"
            print("%s:  version %s" % (mod.capitalize(), version))
            
    
    print(" ")
    
    
    #
    # Library checks
    #
    
    ##  Via 'pkg-config'
    libsFound = []
    for pkgName in ('fftw3f',):
        try:
            pkgQuery = subprocess.Popen(['pkg-config', '--exists', pkgName])
            o, e = pkgQuery.communicate()
            try:
                o = o.decode()
                e = e.decode()
            except AttributeError:
                pass
            if pkgQuery.returncode == 0:
                pkgQuery = subprocess.Popen(['pkg-config', '--modversion', pkgName], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                o, e = pkgQuery.communicate()
                try:
                    o = o.decode()
                    e = e.decode()
                except AttributeError:
                    pass
                o = o.replace('\n', '')
                
                print("%s:  version %s" % (pkgName, o))
                libsFound.append( pkgName )
                
        except OSError:
            pass
            
    ## Via 'ldconfig'
    try:
        p = subprocess.Popen(['ldconfig', '-v'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        o, e = p.communicate()
        try:
            o = o.decode()
            e = e.decode()
        except AttributeError:
            pass
        o = o.split('\n')
        
        for lib in ('libfftw3f', 'libgdbm', 'librt'):
            libBaseName = lib.replace('lib', '')
            if libBaseName in libsFound:
                continue
                
            found = False
            currPath = None
            
            for line in o:
                if len(line) == 0:
                    continue
                elif line[0] != '\t':
                    currPath, junk = line.split(':', 1)
                    continue
                elif line.find(lib) != -1:
                    found = True
                    libFilename = line.split(None, 1)[0]
                    print("%s: found %s" % (lib, os.path.join(currPath, libFilename)))
            if not found:
                print("%s: WARNING - not found" % lib)
            
    except OSError:
        pass
    print(" ")
    
    #
    # Compiler check
    #
    import shutil
    import tempfile
    from distutils import sysconfig
    from distutils import ccompiler
    compiler = ccompiler.new_compiler()
    sysconfig.get_config_vars()
    sysconfig.customize_compiler(compiler)
    cc = compiler.compiler
    
    print("Compiler: %s" % cc[0])
    
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)
    
    fh = open('test.c', 'w')
    fh.write(r"""#include <omp.h>
#include <stdio.h>
int main(void) {
#pragma omp parallel
printf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads());
return 0;
}
""")
    fh.close()
    
    ccmd = []
    ccmd.extend( cc )
    ccmd.extend( ['-fopenmp', 'test.c', '-o test'] )
    if os.path.basename(cc[0]).find('gcc') != -1:
        ccmd.append( '-lgomp' )
    elif os.path.basename(cc[0]).find('clang') != -1:
        ccmd.extend( ['-L/opt/local/lib/libomp', '-lomp'] )
    openmp_support = "No"
    try:
        p = subprocess.Popen(ccmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        o, e = p.communicate()
        try:
            o = o.decode()
            e = e.decode()
        except AttributeError:
            pass
        openmp_support =" Yes" if p.returncode == 0 else "No"
    except subprocess.CalledProcessError:
        pass
    print("Compiler OpenMP Support: %s" % openmp_support)
    if openmp_support == 'Yes':
        o = o.split('\n')[:-1]
        for i in range(len(o)):
            o[i] = '  %s' % o[i]
        o = '\n'.join(o)
        e = e.split('\n')[:-1]
        for i in range(len(e)):
            e[i] = '  %s' % e[i]
        e = '\n'.join(e)
        
        print("Compiler OpenMP Test Command:")
        print("  %s" % ' '.join(ccmd))
        print("Compiler OpenMP Test Output:")
        print(o)
        print("Compiler OpenMP Test Errors:")
        print(e)
        
    os.chdir(curdir)
    shutil.rmtree(tmpdir)
    
    p = subprocess.Popen([cc[0], '-v'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    o, e = p.communicate()
    try:
        o = o.decode()
        e = e.decode()
    except AttributeError:
        pass
    e = e.split('\n')[:-1]
    for i in range(len(e)):
        e[i] = '  %s' % e[i]
    e = '\n'.join(e)
    print("Compiler Version:")
    print(e)
    print(" ")
    
    #
    # Numpy
    #
    try:
        import numpy
        nfp,junk = os.path.split(numpy.__file__)
        nfp = os.path.join(nfp, 'core', 'umath.so')
        nfp = os.path.realpath(nfp)

        p = subprocess.Popen(['file', nfp], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        o, e = p.communicate()
        try:
            o = o.decode()
            e = e.decode()
        except AttributeError:
            pass
        junk, numpyLinkage = o.split(None, 1)
        numpyLinkage = numpyLinkage.replace('\n', '')

        nfp, junk = os.path.split(numpy.__file__)
        print("Numpy Path: %s" % nfp)
        print("Numpy Version: %s" % numpy.version.version)
        print("Numpy Linkage: %s" % numpyLinkage)
    except ImportError as e:
        print("Numpy Import Error: %s" % str(e))
    print(" ")
    
    #
    # LSL
    #
    try:
        import lsl, lsl.version
        lfp,junk = os.path.split(lsl.__file__)
        lfp = os.path.join(lfp, 'correlator', '_core.so')
        lfp = os.path.realpath(lfp)

        p = subprocess.Popen(['file', lfp], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        o, e = p.communicate()
        try:
            o = o.decode()
            e = e.decode()
        except AttributeError:
            pass
        junk, lslLinkage = o.split(None, 1)
        lslLinkage = lslLinkage.replace('\n', '')

        lfp, junk = os.path.split(lsl.__file__)
        print("LSL Path: %s" % lfp)
        print("LSL Version: %s" % lsl.version.version)
        print("LSL Linkage: %s" % lslLinkage)
    except ImportError as e:
        print("LSL Import Error: %s" % str(e))
    print(" ")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='gather information about the Python interpreter, modules, C libraries, numpy installation, and LSL installation', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    args = parser.parse_args()
    main(args)

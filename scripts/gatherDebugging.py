#!/usr/bin/env python3

"""
Script to gather information about the Python interpreter, modules, 
C libraries, numpy installation, and LSL installation to help with 
debugging and install issues.
"""

import os
import re
import sys
import glob
import argparse
import platform
import subprocess
import importlib.util

try:
    from lsl.misc import telemetry
    telemetry.track_script()
except ImportError:
    pass


def main(args):
    #
    # Clean the path
    #
    sys.path = sys.path[1:]

    #
    # Python interpreter
    #
    print(f"Executable path: {sys.executable}")
    print(f"Platform: {sys.platform}")
    print(f"Version: {sys.version}")
    print(f"API: {sys.api_version}")
    print("Bits: %s\nLinkage: %s" % platform.architecture())
    print(" ")
    
    #
    # Python Module Check
    #
    
    ## Required
    for mod in ('numpy', 'scipy', 'astropy', 'ephem', 'aipy', 'pytz'):
        modSpec = importlib.util.find_spec(mod, [os.path.dirname(__file__)])
        if modSpec != None:
            modInfo = importlib.util.module_from_spec(modSpec)
            try:
                modSpec.loader.exec_module(modInfo)
            except ImportError:
                print(f"{mod}: not found")
                continue
            if hasattr(modSpec.loader, 'file'):
                modSpec.loader.file.close()
                
            try:
                version = modInfo.version.version
            except AttributeError:
                try:
                    version = modInfo.__version__
                except AttributeError:
                    try:	
                        versionRE = re.compile(r'%s-(?P<version>[\d\.]+)-py.*' % mod)
                        mtch = versionRE.search(modInfo.__file__)
                        version = mtch.group('version')
                    except (AttributeError, IndexError):
                        version = "unknown"
            print(f"{mod}:  version {version}")
        else:
            print(f"{mod}: not found")
            
    ## Optional
    for mod in ('casacore', 'matplotlib', 'h5py', 'psrfits_utils'):
        if modSpec != None:
            modInfo = importlib.util.module_from_spec(modSpec)
            modSpec.loader.exec_module(modInfo)
            if hasattr(modSpec.loader, 'file'):
                modSpec.loader.file.close()
                
            try:
                version = modInfo.version.version
            except AttributeError:
                try:
                    version = modInfo.__version__
                except AttributeError:
                    try:	
                        versionRE = re.compile(r'%s-(?P<version>[\d\.]+)-py.*' % mod)
                        mtch = versionRE.search(modInfo.__file__)
                        version = mtch.group('version')
                    except (AttributeError, IndexError):
                        version = "unknown"
            print(f"{mod}:  version {version}")
        else:
            print(f"{mod}: not found")
            
    print(" ")
    
    
    #
    # Library checks
    #
    
    ##  Via 'pkg-config'
    libsFound = []
    for pkgName in ('fftw3f',):
        try:
            subprocess.check_output(['pkg-config', '--exists', pkgName])
            o = subprocess.check_output(['pkg-config', '--modversion', pkgName], stderr=subprocess.DEVNULL)
            o = o.decode(encoding='ascii', errors='ignore').replace('\n', '')
            
            print(f"{pkgName}:  version {o}")
            libsFound.append( pkgName )
            
        except (OSError, subprocess.CalledProcessError):
            pass
            
    ## Via 'ldconfig'
    try:
        p = subprocess.run(['ldconfig', '-v'], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=False)
        o = p.stdout.decode().split('\n')
        
        for lib in ('libfftw3f', 'libgdbm', 'librt', 'libgsl'):
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
                    print(f"{lib}: found {os.path.join(currPath, libFilename)}")
            if not found:
                print(f"{lib}: WARNING - not found")
            
    except OSError:
        pass
    print(" ")
    
    #
    # Compiler check
    #
    try:
        import tempfile
        from setuptools import Distribution, Extension
        
        d = Distribution()
        b = d.get_command_obj('build_ext')
        b.finalize_options()
        b.extensions = [Extension('test', ['test.c'])]
        b.build_extensions = lambda: None
        b.run()

        cc = b.compiler
        print(f"Compiler: {cc.compiler[0]}")
        
        is_gcc = (os.path.basename(cc.compiler[0]).find('gcc') != -1)
        is_clang = (os.path.basename(cc.compiler[0]).find('clang') != -1)
        
        ipath = ''
        lpath = ''
        if is_clang:
            # Crude fix for clang + homebrew
            try:
                _ipath = subprocess.check_output(['find', '/opt/homebrew/Cellar/libomp', '-name', 'omp.h'])
                if _ipath != b'':
                    ipath = os.path.dirname(_ipath.decode())
            except subprocess.CalledProcessError:
                pass
            try:
                _lpath = subprocess.check_output(['find', '/opt/homebrew/Cellar/libomp', '-name', 'libomp.a'])
                if _lpath != b'':
                    lpath = os.path.dirname(_lpath.decode())
            except subprocess.CalledProcessError:
                pass
                
        outCFLAGS = []
        outLIBS = []
        if is_clang:
            if ipath != '':
                outCFLAGS.append( '-I'+ipath )
            outCFLAGS.append( '-Xclang' )
        outCFLAGS.append( '-fopenmp' )
        if is_gcc:
            outLIBS.append( '-lgomp' )
        elif is_clang:
            if lpath != '':
                outLIBS.append( '-L'+lpath )
            outLIBS.append( '-lomp' )
            
        curdir = os.getcwd()
        with tempfile.TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)
            
            with open('test.c', 'w') as fh:
                fh.write(r"""#include <omp.h>
#include <stdio.h>
int main(void) {
#pragma omp parallel
printf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads());
return 0;
}
""")
            
            openmp_support = "No"
            try:
                cc.compile(['test.c'], extra_postargs=outCFLAGS)
                cc.link_executable(['test.o'], 'test', extra_postargs=outLIBS)
                openmp_support = " Yes"
            except Exception:
                pass
            print(f"Compiler OpenMP Support: {openmp_support}")
            
        os.chdir(curdir)
        p = subprocess.Popen([cc.compiler[0], '-v'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _, e = p.communicate()
        e = e.decode(encoding='ascii', errors='ignore').split('\n')[:-1]
        for i in range(len(e)):
            e[i] = '  %s' % e[i]
        e = '\n'.join(e)
        print("Compiler Version:")
        print(e)
        print(" ")
    except ImportError as e:
        print(f"Cannot determine compiler OpenMP support: {str(e)}")
    print(" ")
    
    #
    # Numpy
    #
    try:
        import numpy
        nfp,junk = os.path.split(numpy.__file__)
        nfp = glob.glob(os.path.join(nfp, 'core', 'umath.*'))[0]
        nfp = os.path.realpath(nfp)

        o = subprocess.check_output(['file', nfp], stderr=subprocess.DEVNULL)
        o = o.decode(encoding='ascii', errors='ignore')
        junk, numpyLinkage = o.split(None, 1)
        numpyLinkage = numpyLinkage.replace('\n', '')

        nfp, junk = os.path.split(numpy.__file__)
        print(f"Numpy Path: {nfp}")
        print(f"Numpy Version: {numpy.version.version}")
        print(f"Numpy Linkage: {numpyLinkage}")
    except (OSError, subprocess.CalledProcessError):
        print(f"Numpy Path: {nfp}")
        print(f"Numpy Version: {numpy.version.version}")
        print("Numpy Linkage: Unknown")
    except ImportError as e:
        print(f"Numpy Import Error: {str(e)}")
    print(" ")
    
    #
    # LSL
    #
    try:
        import lsl, lsl.version
        lfp,junk = os.path.split(lsl.__file__)
        lfp = glob.glob(os.path.join(lfp, 'correlator', '_core.*'))[0]
        lfp = os.path.realpath(lfp)

        o = subprocess.check_output(['file', lfp], stderr=subprocess.DEVNULL)
        o = o.decode(encoding='ascii', errors='ignore')
        junk, lslLinkage = o.split(None, 1)
        lslLinkage = lslLinkage.replace('\n', '')

        lfp, junk = os.path.split(lsl.__file__)
        print(f"LSL Path: {lfp}")
        print(f"LSL Version: {lsl.version.version}")
        print(f"LSL Linkage: {lslLinkage}")
    except (OSError, subprocess.CalledProcessError):
        print(f"LSL Path: {lfp}")
        print(f"LSL Version: {lsl.version.version}")
        print(f"LSL Linkage: Unknown")
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

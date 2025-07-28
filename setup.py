
import os
import sys
import glob
import logging
import tempfile
import subprocess

from setuptools import setup, Distribution, Extension, find_namespace_packages
try:
    from setuptools.command.build import build
except ImportError:
    from distutils.command.build import build
from setuptools.command.build_ext import build_ext

log = logging.getLogger('__main__')

try:
    import numpy
except Exception as e:
    raise RuntimeError(f"numpy is required to run setup.py: {str(e)}")


def get_version():
    """Read the VERSION file and return the version number as a string."""

    with open('VERSION', 'r') as fh:
        version = fh.read().strip()
    return version


def get_description(filename):
    """Read in a README-type file and return the contents of the DESCRIPTION
    section."""

    desc = ''
    with open(filename, 'r') as fh:
        lines = fh.readlines()
        
    inDescription = False
    for line in lines:
        if line.find('DESCRIPTION') == 0:
            inDescription = True
            continue
        if line.find('REQUIREMENTS') == 0:
            inDescription = False
            break
        if inDescription:
            if line[:3] == '---':
                continue
            desc = ''.join([desc, line])

    return desc


def get_openmp():
    """Try to compile/link an example program to check for OpenMP support.
    
    Based on:
     1) http://stackoverflow.com/questions/16549893/programatically-testing-for-openmp-support-from-a-python-setup-script
     2) https://github.com/lpsinger/healpy/blob/6c3aae58b5f3281e260ef7adce17b1ffc68016f0/setup.py
     3) https://github.com/pypa/setuptools/issues/2806#issuecomment-961805789
    """
    
    d = Distribution()
    b = d.get_command_obj('build_ext')
    b.finalize_options()
    b.extensions = [Extension('test', ['test.c'])]
    b.build_extensions = lambda: None
    b.run()
    
    # Compiler selection/specific paths
    cc = b.compiler
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
            
    # Compiler and linker flags
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
            
        try:
            cc.compile(['test.c'], extra_postargs=outCFLAGS)
            cc.link_executable(['test.o'], 'test', extra_postargs=outLIBS)
        except Exception as e:
            print(f"WARNING:  OpenMP does not appear to be supported by {cc.compiler[0]}, disabling")
            print(f"WARNING: {str(e)}")
            outCFLAGS = []
            outLIBS = []
            
        os.chdir(curdir)
        
    return outCFLAGS, outLIBS


def get_gsl():
    """Use pkg-config (if installed) to figure out the C flags and linker flags
    needed to compile a C program with GSL.  If GSL cannot be found via
    pkg-config, some 'sane' values are returned."""
    
    try:
        subprocess.check_call(['pkg-config', 'gsl', '--exists'])
        
        p = subprocess.Popen(['pkg-config', 'gsl', '--modversion'], stdout=subprocess.PIPE)
        outVersion = p.communicate()[0].decode()
        outVersion = outVersion.rstrip().split()
        
        p = subprocess.Popen(['pkg-config', 'gsl', '--cflags'], stdout=subprocess.PIPE)
        outCFLAGS = p.communicate()[0].decode()
        outCFLAGS = outCFLAGS.rstrip().split()
        try:
            outCFLAGS = [str(v, 'utf-8') for v in outCFLAGS]
        except TypeError:
            pass
        
        p = subprocess.Popen(['pkg-config', 'gsl', '--libs'], stdout=subprocess.PIPE)
        outLIBS = p.communicate()[0].decode()
        outLIBS = outLIBS.rstrip().split()
        try:
            outLIBS = [str(v, 'utf-8') for v in outLIBS]
        except TypeError:
            pass
            
        if len(outVersion) > 0:
            print("Found GSL, version %s" % outVersion[0])
            
    except (OSError, subprocess.CalledProcessError):
        print("WARNING:  GSL cannot be found, using defaults")
        outCFLAGS = []
        outLIBS = ['-lgsl', '-lm']
        
    return outCFLAGS, outLIBS


def get_fftw():
    """Use pkg-config (if installed) to figure out the C flags and linker flags
    needed to compile a C program with single precision FFTW3.  If FFTW3 cannot 
    be found via pkg-config, some 'sane' values are returned."""
    
    try:
        subprocess.check_call(['pkg-config', 'fftw3f', '--exists'])
        
        p = subprocess.Popen(['pkg-config', 'fftw3f', '--modversion'], stdout=subprocess.PIPE)
        outVersion = p.communicate()[0].decode()
        outVersion = outVersion.strip().split()
        
        p = subprocess.Popen(['pkg-config', 'fftw3f', '--cflags'], stdout=subprocess.PIPE)
        outCFLAGS = p.communicate()[0].decode()
        outCFLAGS = outCFLAGS.rstrip().split()
        try:
            outCFLAGS = [str(v) for v in outCFLAGS]
        except TypeError:
            pass
        
        p = subprocess.Popen(['pkg-config', 'fftw3f', '--libs'], stdout=subprocess.PIPE)
        outLIBS = p.communicate()[0].decode()
        outLIBS = outLIBS.rstrip().split()
        try:
            outLIBS = [str(v) for v in outLIBS]
        except TypeError:
            pass
            
        if len(outVersion) > 0:
            print(f"Found FFTW3, version {outVersion[0]}")
            
    except (OSError, subprocess.CalledProcessError):
        print("WARNING:  single precision FFTW3 cannot be found, using defaults")
        outCFLAGS = []
        outLIBS = ['-lfftw3f', '-lm']
        
    return outCFLAGS, outLIBS


def write_version_info():
    """Write the version info to a module in LSL."""
    
    lslVersion = get_version()
    shortVersion = '.'.join(lslVersion.split('.')[:2])
    
    contents = f"""# This file is automatically generated by setup.py

import os
import glob
import hashlib

def _get_md5(filename, blockSize=262144):
    with open(filename, 'rb') as fh:
        m = hashlib.md5()
        while True:
            block = fh.read(blockSize)
            if len(block) == 0:
                break
            m.update(block)
            
    return m.hexdigest()

def get_fingerprint():
    \"\"\"
    Return a 'fingerprint' of the current LSL module state that is useful for
    tracking if the module has changed.
    \"\"\"
    
    filenames = []
    for ext in ('.py', '.so'):
        filenames.extend( glob.glob(os.path.join(os.path.dirname(__file__), f"*{{ext}}")) )
        filenames.extend( glob.glob(os.path.join(os.path.dirname(__file__), '*', f"*{{ext}}")) )
        
    m = hashlib.md5()
    for filename in filenames:
        m.update(_get_md5(filename).encode())
    return m.hexdigest()

version = '{lslVersion}'
full_version = '{lslVersion}'
short_version = '{shortVersion}'

"""
    
    with open('lsl/version/__init__.py', 'w') as fh:
        fh.write(contents)
        
    return True


# Get the FFTW3 flags/libs and manipulate the flags and libraries for 
# correlator._core appropriately.  This will, hopefully, fix the build
# problems on Mac
class lsl_build(build):
    user_options = build.user_options \
                   + [('with-gsl=', None, 'Installation path for GSL'),] \
                   + [('with-fftw=', None, 'Installation path for single precision FFTW3'),]
                   
    
    def initialize_options(self, *args, **kwargs):
        build.initialize_options(self, *args, **kwargs)
        self.with_gsl = None
        self.with_fftw = None
        
    def finalize_options(self, *args, **kwargs):
        build.finalize_options(self, *args, **kwargs)
        
        if self.distribution.has_ext_modules():
            ## Grab the 'build_ext' command
            beco = self.distribution.get_command_obj('build_ext')
            
            ## Grab the GSL flags
            if os.getenv('WITH_GSL', None) is not None:
                self.with_gsl = os.getenv('WITH_GSL')
            if self.with_gsl is not None:
                beco.with_gsl = self.with_gsl
            
            ## Grab the FFTW flags
            if os.getenv('WITH_FFTW', None) is not None:
                self.with_fftw = os.getenv('WITH_FFTW')
            if self.with_fftw is not None:
                beco.with_fftw = self.with_fftw


class lsl_build_ext(build_ext):
    user_options = build_ext.user_options \
                   + [('with-gls=', None, 'Installation path for GSL'),] \
                   + [('with-fftw=', None, 'Installation path for single precision FFTW3'),]
    
    def initialize_options(self, *args, **kwargs):
        build_ext.initialize_options(self, *args, **kwargs)
        self.with_gsl = None
        self.with_fftw = None
        
    def finalize_options(self, *args, **kwargs):
        build_ext.finalize_options(self, *args, **kwargs)
        self.verbose = True
        
        ## Grab the OpenMP flags
        openmpFlags, openmpLibs = get_openmp()
        
        ## Grab the GSL flags
        if os.getenv('WITH_GSL', None) is not None:
            self.with_gsl = os.getenv('WITH_GSL')
        if self.with_gsl is not None:
            gslFlags = ['-I%s/include' % self.with_gsl,]
            gslLibs = ['-L%s/lib' % self.with_gsl, '-lgsl']
        else:
            gslFlags, gslLibs = get_gsl()
        
        ## Grab the FFTW flags
        if os.getenv('WITH_FFTW', None) is not None:
            self.with_fftw = os.getenv('WITH_FFTW')
        if self.with_fftw is not None:
            fftwFlags = [f"-I{self.with_fftw}/include",]
            fftwLibs = [f"-L{self.with_fftw}/lib", '-lfftw3f']
        else:
            fftwFlags, fftwLibs = get_fftw()
            
        ## Update the extensions with the additional compilier/linker flags
        for ext in self.extensions:
            ### Compiler flags
            for cflags in (openmpFlags, gslFlags, fftwFlags):
                try:
                    ext.extra_compile_args.extend( cflags )
                except TypeError:
                    ext.extra_compile_args = cflags
            ### Linker flags
            for ldflags in (openmpLibs, gslLibs, fftwLibs):
                try:
                    ext.extra_link_args.extend( ldflags )
                except TypeError:
                    ext.extra_link_args = ldflags
                    
        ## HACK: Update the log verbosity - for some reason this gets set to 
        ##       WARN when I replace build_ext
        try:
            log.setLevel(min([logging.INFO, log.level]))
        except AttributeError:
            pass


coreExtraFlags = ['-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION', '--std=c++11']
coreExtraLibs = []


# Create the list of extension modules.  We do this here so that we can turn 
# off the DRSU direct module for non-linux system
ExtensionModules = [Extension('reader._gofast', ['lsl/reader/gofast.cpp', 'lsl/reader/tbn.cpp', 'lsl/reader/drx.cpp', 'lsl/reader/drx8.cpp', 'lsl/reader/drspec.cpp', 'lsl/reader/vdif.cpp', 'lsl/reader/tbf.cpp', 'lsl/reader/cor.cpp', 'lsl/reader/tbx.cpp'], include_dirs=[numpy.get_include()], libraries=['m'], extra_compile_args=coreExtraFlags, extra_link_args=coreExtraLibs),
            Extension('common._fir', ['lsl/common/fir.cpp'], include_dirs=[numpy.get_include()], libraries=['m'], extra_compile_args=coreExtraFlags, extra_link_args=coreExtraLibs),
            Extension('correlator._spec', ['lsl/correlator/spec.cpp'], include_dirs=[numpy.get_include()], libraries=['m'], extra_compile_args=coreExtraFlags, extra_link_args=coreExtraLibs), 
            Extension('correlator._stokes', ['lsl/correlator/stokes.cpp'], include_dirs=[numpy.get_include()], libraries=['m'], extra_compile_args=coreExtraFlags, extra_link_args=coreExtraLibs),
            Extension('correlator._core', ['lsl/correlator/core.cpp'], include_dirs=[numpy.get_include()], libraries=['m'], extra_compile_args=coreExtraFlags, extra_link_args=coreExtraLibs), 
            Extension('imaging._gridder', ['lsl/imaging/gridder.cpp'], include_dirs=[numpy.get_include()], libraries=['m'], extra_compile_args=coreExtraFlags, extra_link_args=coreExtraLibs), 
            Extension('sim._simfast', ['lsl/sim/simfast.cpp'], include_dirs=[numpy.get_include()], libraries=['m'], extra_compile_args=coreExtraFlags, extra_link_args=coreExtraLibs), 
            Extension('misc._wisdom', ['lsl/misc/wisdom.cpp'],include_dirs=[numpy.get_include()], libraries=['m'], extra_compile_args=coreExtraFlags, extra_link_args=coreExtraLibs), ]

# Update the version information
write_version_info()

setup(
    cmdclass = {'build': lsl_build, 'build_ext': lsl_build_ext}, 
    name = "lsl", 
    version = get_version(), 
    description = "LWA Software Library", 
    author = "Jayce Dowell", 
    author_email = "jdowell@unm.edu", 
    url = "https://github.com/lwa-project/lsl", 
    long_description = get_description('README.md'), 
    license = 'GPL',
    classifiers = ['Development Status :: 5 - Production/Stable',
                   'Intended Audience :: Developers',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: GNU General Public License (GPL)',
                   'Topic :: Scientific/Engineering :: Astronomy',
                   'Programming Language :: Python :: 3',
                   'Programming Language :: Python :: 3.8',
                   'Programming Language :: Python :: 3.9',
                   'Programming Language :: Python :: 3.10',
                   'Programming Language :: Python :: 3.11',
                   'Programming Language :: Python :: 3.12',
                   'Operating System :: MacOS :: MacOS X',
                   'Operating System :: POSIX :: Linux'],
    packages = find_namespace_packages(), 
    scripts = glob.glob('scripts/*.py'), 
    python_requires='>=3.8', 
    setup_requires = ['numpy>=1.7'], 
    install_requires = ['astropy>=5.2', 'jplephem', 'healpy', 'h5py', 'numpy>=1.7', 'scipy>=0.19', 'ephem>=3.7.5.3', 'aipy>=3.0.1', 'pytz>=2012c'],
    include_package_data = True,  
    ext_package = 'lsl', 
    ext_modules = ExtensionModules,
    zip_safe = False,  
    test_suite = "tests"
) 

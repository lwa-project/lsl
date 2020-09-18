import numpy

from lsl.complex.numpy_complex_int import complex_int8, complex_int16, complex_int32
from utils import get_include

__version__ = "0.1"
__all__ = ['complex_int8', 'complex_int16', 'complex_int32', 'get_include']

if numpy.__dict__.get('complex_int8') is None:
    numpy.complex_int8 = complex_int8
    numpy.typeDict['complex_int8'] = numpy.dtype(complex_int8)
    
if numpy.__dict__.get('complex_int16') is None:
    numpy.complex_int16 = complex_int16
    numpy.typeDict['complex_int16'] = numpy.dtype(complex_int16)

if numpy.__dict__.get('complex_int32') is None:
    numpy.complex_int326 = complex_int32
    numpy.typeDict['complex_int32'] = numpy.dtype(complex_int32)

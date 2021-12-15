#ifndef CORRELATOR_COMMON_H_INCLUDE_GUARD_
#define CORRELATOR_COMMON_H_INCLUDE_GUARD_

#include "Python.h"
#include <cmath>
#include <complex>
#include "numpy/arrayobject.h"
#include "numpy/npy_math.h"


/*
 Load in FFTW wisdom.  Based on the read_wisdom function in PRESTO.
*/

inline void read_wisdom(char *filename, PyObject *m) {
    int status = 0;
    FILE *wisdomfile;
    
    wisdomfile = fopen(filename, "r");
    if( wisdomfile != NULL ) {
        status = fftwf_import_wisdom_from_file(wisdomfile);
        fclose(wisdomfile);
    }
    PyModule_AddObject(m, "useWisdom", PyBool_FromLong(status));
}


/*
 Wrap around PyArray_ISCOMPLEX to deal with our ci8 data format
*/

#define PyArray_LSL_ISCOMPLEX(arr, minDim) (PyArray_ISCOMPLEX(arr) \
                                            || ((PyArray_TYPE(arr) == NPY_INT8) \
                                                && (PyArray_NDIM(arr) == (minDim+1))))


/*
 Wrap around PyArray_TYPE to deal with our ci8 data format
*/

#define LSL_CI8 262

#define PyArray_LSL_TYPE(arr, minDim) (((PyArray_TYPE(arr) == NPY_INT8) \
                                        && (PyArray_NDIM(arr) == (minDim+1))) \
                                       ? LSL_CI8 : PyArray_TYPE(arr))


/*
  Warp the Numpy PyArray_DATA macro so that it can deal with NULL values.
*/

#define PyArray_SAFE_DATA(arr) (arr != NULL ? PyArray_DATA(arr) : NULL)


/*
  Sinc function for use by the polyphase filter bank
*/

template<typename T>
inline T sinc(T x) {
    if(x == 0.0) {
        return 1.0;
    } else {
        return sin(x*NPY_PI)/(x*NPY_PI);
    }
}


/*
  Hanning window for use by the polyphase filter bank
*/

template<typename T>
inline T hanning(T x) {
    return 0.5 - 0.5*cos(x);
}


/*
  Hamming window for use by the polyphase filter bank
*/

template<typename T>
inline T hamming(T x) {
    return 0.53836 - 0.46164*cos(x); 
}


/*
  Number of PFB taps to use
 */

#define PFB_NTAP 4


/*
  Complex types
*/

typedef std::complex<float> Complex32;
typedef std::complex<double> Complex64;


/*
  Macro for 2*pi
*/

#define TPI (2*NPY_PI*Complex64(0,1))


/*
  Complex magnitude squared functions
*/

template<typename T>
inline T abs2(std::complex<T> z) {
    return z.real()*z.real() + z.imag()*z.imag();
}

#endif // CORRELATOR_COMMON_H_INCLUDE_GUARD_

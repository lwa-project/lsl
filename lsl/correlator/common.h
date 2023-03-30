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
  Warp the Numpy PyArray_DATA macro so that it can deal with NULL values.
*/

#define PyArray_SAFE_DATA(arr)   (arr != NULL ? PyArray_DATA(arr) : NULL)


/*
  Sinc function for use by the polyphase filter bank
*/

inline double sinc(double x) {
    if(x == 0.0) {
        return 1.0;
    } else {
        return sin(x*NPY_PI)/(x*NPY_PI);
    }
}

inline float sinc(float x) {
    if(x == 0.0) {
        return 1.0;
    } else {
        return sin(x*NPY_PI)/(x*NPY_PI);
    }
}


/*
  Hanning window for use by the polyphase filter bank
*/

inline double hanning(double x) {
    return 0.5 - 0.5*cos(x);
    
}

inline float hanning(float x) {
    return 0.5 - 0.5*cos(x);
}


/*
  Hamming window for use by the polyphase filter bank
*/

inline double hamming(double x) {
    return 0.543478261 - 0.456521739*cos(x);
    
}

inline float hamming(float x) {
    return 0.543478261 - 0.456521739*cos(x);
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

inline float abs2(Complex32 z) {
    return z.real()*z.real() + z.imag()*z.imag();
}
inline double abs2(Complex64 z) {
    return z.real()*z.real() + z.imag()*z.imag();
}

#endif // CORRELATOR_COMMON_H_INCLUDE_GUARD_

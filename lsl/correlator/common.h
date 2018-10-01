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

inline float sinc(float x) {
    if(x == 0.0) {
        return 1.0;
    } else {
        return sin(x*NPY_PI)/(x*NPY_PI);
    }
}
inline double sinc(double x) {
    if(x == 0.0) {
        return 1.0;
    } else {
        return sin(x*NPY_PI)/(x*NPY_PI);
    }
}


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
  fftwf_complex math helpers
*/


inline float real(fftwf_complex z) {
    return z[0];
}
inline float imag(fftwf_complex z) {
    return z[1];
}
inline float abs(fftwf_complex z) {
    return std::sqrt(z[0]*z[0] + z[1]*z[1]);
}
inline float abs2(fftwf_complex z) {
    return z[0]*z[0] + z[1]*z[1];
}


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
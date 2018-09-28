#ifndef CORRELATOR_COMMON_H_INCLUDE_GUARD_
#define CORRELATOR_COMMON_H_INCLUDE_GUARD_

#include "Python.h"
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
  Macro for 2*pi
*/

#define TPI (2*NPY_PI*_Complex_I)


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
  Complex magnitude squared functions
*/

inline float cabs2(float complex z) {
    return crealf(z)*crealf(z) + cimagf(z)*cimagf(z);
}
inline double cabs2(double complex z) {
    return creal(z)*creal(z) + cimag(z)*cimag(z);
}

#endif // CORRELATOR_COMMON_H_INCLUDE_GUARD_
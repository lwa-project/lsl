#pragma once

#include "Python.h"
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fftw3.h>
#include "numpy/arrayobject.h"
#include "numpy/npy_math.h"

/*
 64-byte aligned memory allocator/deallocator
*/

inline void* aligned64_malloc(size_t size) {
  void *ptr = NULL;
  int err = posix_memalign(&ptr, 64, size);
  if( err != 0 ) {
    return NULL;
  }
  return ptr;
}

inline void aligned64_free(void* ptr) {
  free(ptr);
}


/*
 FFTW precision selection
*/

#if defined(USE_FFTW_DOUBLE) && USE_FFTW_DOUBLE
    typedef double real_t;
    typedef std::complex<double> complex_t;
    typedef fftw_plan fftw_plan_t;
    typedef fftw_complex fftw_complex_t;
    #define NPY_REAL_T NPY_DOUBLE
    #define NPY_COMPLEX_T NPY_COMPLEX128

    #define FFTW_MALLOC fftw_malloc
    #define FFTW_FREE fftw_free
    #define FFTW_PLAN_MANY_DFT_R2C fftw_plan_many_dft_r2c
    #define FFTW_PLAN_MANY_DFT fftw_plan_many_dft
    #define FFTW_PLAN_DFT_1D fftw_plan_dft_1d
    #define FFTW_PLAN_DFT_2D fftw_plan_dft_2d
    #define FFTW_PLAN_DFT_R2C_1D fftw_plan_dft_r2c_1d
    #define FFTW_PLAN_R2R_1D fftw_plan_r2r_1d
    #define FFTW_EXECUTE fftw_execute
    #define FFTW_EXECUTE_DFT_R2C fftw_execute_dft_r2c
    #define FFTW_EXECUTE_DFT fftw_execute_dft
    #define FFTW_DESTROY_PLAN fftw_destroy_plan
    #define FFTW_IMPORT_SYSTEM_WISDOM fftw_import_system_wisdom
    #define FFTW_IMPORT_WISDOM_FROM_FILE fftw_import_wisdom_from_file
    #define FFTW_EXPORT_WISDOM_TO_FILE fftw_export_wisdom_to_file
#else
    typedef float real_t;
    typedef std::complex<float> complex_t;
    typedef fftwf_plan fftw_plan_t;
    typedef fftwf_complex fftw_complex_t;
    #define NPY_REAL_T NPY_FLOAT
    #define NPY_COMPLEX_T NPY_COMPLEX64

    #define FFTW_MALLOC fftwf_malloc
    #define FFTW_FREE fftwf_free
    #define FFTW_PLAN_MANY_DFT_R2C fftwf_plan_many_dft_r2c
    #define FFTW_PLAN_MANY_DFT fftwf_plan_many_dft
    #define FFTW_PLAN_DFT_1D fftwf_plan_dft_1d
    #define FFTW_PLAN_DFT_2D fftwf_plan_dft_2d
    #define FFTW_PLAN_DFT_R2C_1D fftwf_plan_dft_r2c_1d
    #define FFTW_PLAN_R2R_1D fftwf_plan_r2r_1d
    #define FFTW_EXECUTE fftwf_execute
    #define FFTW_EXECUTE_DFT_R2C fftwf_execute_dft_r2c
    #define FFTW_EXECUTE_DFT fftwf_execute_dft
    #define FFTW_DESTROY_PLAN fftwf_destroy_plan
    #define FFTW_IMPORT_SYSTEM_WISDOM fftwf_import_system_wisdom
    #define FFTW_IMPORT_WISDOM_FROM_FILE fftwf_import_wisdom_from_file
    #define FFTW_EXPORT_WISDOM_TO_FILE fftwf_export_wisdom_to_file
#endif


/*
 Load in FFTW wisdom.  Based on the read_wisdom function in PRESTO.
*/

inline void read_wisdom(char *filename, PyObject *m) {
    int status = 0;
    FILE *wisdomfile;
    
    wisdomfile = fopen(filename, "r");
    if( wisdomfile != NULL ) {
        status = FFTW_IMPORT_WISDOM_FROM_FILE(wisdomfile);
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
    return 0.54 - 0.46*cos(x); 
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

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
#define LSL_fft_plan fftw_plan
#define LSL_fft_rtype double
#define LSL_fft_ctype std::complex<double>
#define LSL_fft_complex fftw_complex
#define LSL_fft_np_rtype NPY_DOUBLE
#define LSL_fft_np_ctype NPY_COMPLEX128

#define LSL_fft_malloc fftw_malloc
#define LSL_fft_free fftw_free
#define LSL_fft_plan_many_dft_r2c fftw_plan_many_dft_r2c
#define LSL_fft_plan_many_dft fftw_plan_many_dft
#define LSL_fft_plan_dft_1d fftw_plan_dft_1d
#define LSL_fft_plan_dft_2d fftw_plan_dft_2d
#define LSL_fft_plan_dft_r2c_1d fftw_plan_dft_r2c_1d
#define LSL_fft_plan_r2r_1d fftw_plan_r2r_1d
#define LSL_fft_execute fftw_execute
#define LSL_fft_execute_dft_r2c fftw_execute_dft_r2c
#define LSL_fft_execute_dft fftw_execute_dft
#define LSL_fft_destroy_plan fftw_destroy_plan
#define LSL_fft_import_system_wisdom fftw_import_system_wisdom
#define LSL_fft_import_wisdom_from_file fftw_import_wisdom_from_file
#define LSL_fft_export_wisdom_to_file fftw_export_wisdom_to_file

#else
#define LSL_fft_plan fftwf_plan
#define LSL_fft_rtype float
#define LSL_fft_ctype std::complex<float>
#define LSL_fft_complex fftwf_complex
#define LSL_fft_np_rtype NPY_FLOAT
#define LSL_fft_np_ctype NPY_COMPLEX64

#define LSL_fft_malloc fftwf_malloc
#define LSL_fft_free fftwf_free
#define LSL_fft_plan_many_dft_r2c fftwf_plan_many_dft_r2c
#define LSL_fft_plan_many_dft fftwf_plan_many_dft
#define LSL_fft_plan_dft_1d fftwf_plan_dft_1d
#define LSL_fft_plan_dft_2d fftwf_plan_dft_2d
#define LSL_fft_plan_dft_r2c_1d fftwf_plan_dft_r2c_1d
#define LSL_fft_plan_r2r_1d fftwf_plan_r2r_1d
#define LSL_fft_execute fftwf_execute
#define LSL_fft_execute_dft_r2c fftwf_execute_dft_r2c
#define LSL_fft_execute_dft fftwf_execute_dft
#define LSL_fft_destroy_plan fftwf_destroy_plan
#define LSL_fft_import_system_wisdom fftwf_import_system_wisdom
#define LSL_fft_import_wisdom_from_file fftwf_import_wisdom_from_file
#define LSL_fft_export_wisdom_to_file fftwf_export_wisdom_to_file

#endif


/*
 Load in FFTW wisdom.  Based on the read_wisdom function in PRESTO.
*/

inline void read_wisdom(char *filename, PyObject *m) {
    int status = 0;
    FILE *wisdomfile;
    
    wisdomfile = fopen(filename, "r");
    if( wisdomfile != NULL ) {
        status = LSL_fft_import_wisdom_from_file(wisdomfile);
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

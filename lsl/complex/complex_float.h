#ifndef COMPLEX_COMPLEX_FLOAT_H_INCLUDE_GUARD_
#define COMPLEX_COMPLEX_FLOAT_H_INCLUDE_GUARD_

#ifdef __cplusplus
extern "C" {
#endif

#include <Python.h>
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>

typedef npy_cfloat complex_float32;
typedef npy_cdouble complex_float64;

static NPY_INLINE PyObject* PyComplex_FromComplexFloat32(complex_float32 c) {
    return PyComplex_FromDoubles(c.real, c.imag);
}

static NPY_INLINE PyObject* PyComplex_FromComplexFloat64(complex_float64 c) {
    return PyComplex_FromDoubles(c.real, c.imag);
}

#ifdef __cplusplus
}
#endif

#endif  //COMPLEX_COMPLEX_FLOAT_H_INCLUDE_GUARD_


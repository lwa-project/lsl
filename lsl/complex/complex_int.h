#ifndef COMPLEX_COMPLEX_INT_H_INCLUDE_GUARD_
#define COMPLEX_COMPLEX_INT_H_INCLUDE_GUARD_

#ifdef __cplusplus
extern "C" {
#endif

#include <Python.h>
#include "numpy/arrayobject.h"
#include "numpy/npy_math.h"

#include "complex_int8.h"
#include "complex_int16.h"
#include "complex_int32.h"

#define NPY_COMPLEX_INT8 256
#define NPY_COMPLEX_INT16 257
#define NPY_COMPLEX_INT32 258

static NPY_INLINE int import_complex_int(void) {
    // import_array();
    
    PyObject *module = PyImport_ImportModule("lsl.complex");
    if( module == NULL ) {
        PyErr_Warn(PyExc_RuntimeWarning, "Cannot load the LSL complex integer types");
        return -1;
    }
    Py_XDECREF(module);

    return 0;
}

void lsl_unpack_ci8(complex_int8 packed, signed char* real, signed char* imag);
void lsl_unpack_ci16(complex_int16 packed, signed char* real, signed char* imag);
void lsl_unpack_ci32(complex_int32 packed, short int* real, short int* imag);

#ifdef __cplusplus
}
#endif

#endif  //COMPLEX_COMPLEX_INT_H_INCLUDE_GUARD_



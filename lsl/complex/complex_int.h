#ifndef COMPLEX_COMPLEX_INT_H_INCLUDE_GUARD_
#define COMPLEX_COMPLEX_INT_H_INCLUDE_GUARD_

#ifdef __cplusplus
extern "C" {
#endif

#include <Python.h>
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>
#include <numpy/ufuncobject.h>

#include "complex_int8.h"
#include "complex_int16.h"
#include "complex_int32.h"

#define NPY_COMPLEX_INT8 256
#define NPY_COMPLEX_INT16 257
#define NPY_COMPLEX_INT32 258

static NPY_INLINE int import_complex_int(void) {
    // import_array();
    // import_umath();
    
    PyObject *module = PyImport_ImportModule("lsl.complex");
    if( module == NULL ) {
        PyErr_Warn(PyExc_RuntimeWarning, "Cannot load the LSL complex integer types");
        return -1;
    }
    Py_XDECREF(module);
    
    /* Fill the lookup tables */
    complex_int8_fillLUT();
    
    return 0;
}

#define PyTypeNum_ISCOMPLEX_INT(type) (((type) >= NPY_COMPLEX_INT8) \
                                       && ((type) <= NPY_COMPLEX_INT32))

#define PyDataType_ISCOMPLEX_INT(obj) PyTypeNum_ISCOMPLEX_INT(((PyArray_Descr*)(obj))->type_num)

#define PyArray_ISCOMPLEX_INT(obj) PyTypeNum_ISCOMPLEX_INT(PyArray_TYPE(obj))

#define LSLTypeNum_ISCOMPLEX(type) ((PyTypeNum_ISCOMPLEX(type)) \
                                    || (PyTypeNum_ISCOMPLEX_INT(type)))

#define LSLDataType_ISCOMPLEX(obj) ((PyDataType_ISCOMPLEX(obj)) \
                                    || (PyDataType_ISCOMPLEX_INT(obj)))

#define LSLArray_ISCOMPLEX(obj) ((PyArray_ISCOMPLEX(obj)) \
                                 || (PyArray_ISCOMPLEX_INT(obj)))

#ifdef __cplusplus
}
#endif

#endif  //COMPLEX_COMPLEX_INT_H_INCLUDE_GUARD_



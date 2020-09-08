#define NPY_NO_DEPRECATED_API NPY_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>
#include <numpy/ufuncobject.h>
#include "structmember.h"

#include "../common/py3_compat.h"
#include "complex_int.h"

/*
The supported complex integer types
*/

// Complex 4-bit + 4-bit integer
#include "numpy_complex_int8.c"

// Complex 8-bit + 8-bit integer
#include "numpy_complex_int16.c"

// Complex 16-bit + 16-bit integer
#include "numpy_complex_int32.c"


/*
Module Setup - Function Definitions and Documentation
*/

static PyMethodDef ComplexIntMethods[] = {
    {NULL, NULL, 0, NULL}
};

PyDoc_STRVAR(ComplexIntDoc, \
"Module for representing complex integer data as NumPy complex integer types.\n\
The types supported are:\n\
 * complex_int8  - 4+4-bit complex integer [-8, +7]\n\
 * complex_int16 - 8+8-bit complex integer [-128, +127]\n\
 * complex_int32 - 16+16-bit complex integer [-32768, +32767]");


/*
Module Setup - Initialization
*/

MOD_INIT(numpy_complex_int) {
    PyObject *m, *all, *numpy, *numpy_dict;
    
    Py_Initialize();
    
    // Module definitions and functions
    MOD_DEF(m, "numpy_complex_int", ComplexIntMethods, ComplexIntDoc);
    if( m == NULL ) {
        return MOD_ERROR_VAL;
    }
    
    // Make sure NumPy is initialized
    numpy = PyImport_ImportModule("numpy");
    numpy_dict = PyModule_GetDict(numpy);
    import_array();
    import_umath();
    
    // Register the complexi8 array scalar type
    int complexi8Num = create_complex_int8(m, numpy_dict);
    if( complexi8Num == -2 ) {
        PyErr_Print();
        PyErr_SetString(PyExc_SystemError, "could not initialize PyComplexI8ArrType_Type");
        return MOD_ERROR_VAL;
    } else if( complexi8Num == -1 ) {
        return MOD_ERROR_VAL;
    }
    
    // Register the complexi16 array scalar type
    int complexi16Num = create_complex_int16(m, numpy_dict);
    if( complexi16Num == -2 ) {
        PyErr_Print();
        PyErr_SetString(PyExc_SystemError, "could not initialize PyComplexI16ArrType_Type");
        return MOD_ERROR_VAL;
    } else if( complexi16Num == -1 ) {
        return MOD_ERROR_VAL;
    }
    
    // Register the complexi32 array scalar type
    int complexi32Num = create_complex_int32(m, numpy_dict);
    if( complexi32Num == -2 ) {
        PyErr_Print();
        PyErr_SetString(PyExc_SystemError, "could not initialize PyComplexI32ArrType_Type");
        return MOD_ERROR_VAL;
    } else if( complexi32Num == -1 ) {
        return MOD_ERROR_VAL;
    }
    
    // Version and revision information
    PyModule_AddObject(m, "__version__", PyString_FromString("0.1"));
    
    // Module listing
    all = PyList_New(0);
    PyList_Append(all, PyString_FromString("complex_int8"));
    PyList_Append(all, PyString_FromString("complex_int16"));
    PyList_Append(all, PyString_FromString("complex_int32"));
    PyModule_AddObject(m, "__all__", all);
    
    return MOD_SUCCESS_VAL(m);
}

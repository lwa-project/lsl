#ifndef PY3_COMPAT_H_INCLUDE_GUARD_
#define PY3_COMPAT_H_INCLUDE_GUARD_

#include "Python.h"

/*
  Python3 Compatiability
*/

#if PY_MAJOR_VERSION >= 3
    #define PyCapsule_Type PyCObject_Type
    #define PyInt_AsLong PyLong_AsLong
    #define PyInt_FromLong PyLong_FromLong
    #define PyString_GET_SIZE PyBytes_GET_SIZE
    #define PyString_AS_STRING PyBytes_AS_STRING
    #define PyString_FromString PyUnicode_FromString
inline char* PyString_AsString(PyObject *ob) {
    PyObject *enc;
    char *cstr;
    enc = PyUnicode_AsEncodedString(ob, "utf-8", "Error");
    if( enc == NULL ) {
        PyErr_Format(PyExc_ValueError, "Cannot encode string");
        return NULL;
    }
    cstr = PyBytes_AsString(enc);
    Py_XDECREF(enc);
    return cstr;
}
    #define MOD_ERROR_VAL NULL
    #define MOD_SUCCESS_VAL(val) val
    #define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
    #define MOD_DEF(ob, name, methods, doc) \
       static struct PyModuleDef moduledef = { \
          PyModuleDef_HEAD_INIT, name, doc, -1, methods, }; \
       ob = PyModule_Create(&moduledef);
#else
    #define MOD_ERROR_VAL
    #define MOD_SUCCESS_VAL(val)
    #define MOD_INIT(name) PyMODINIT_FUNC init##name(void)
    #define MOD_DEF(ob, name, methods, doc) \
       ob = Py_InitModule3(name, methods, doc);
#endif

#endif // PY3_COMPAT_H_INCLUDE_GUARD_
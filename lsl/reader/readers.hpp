#pragma once

#include "Python.h"

#if defined(__linux__)
/* Linux */
#include <byteswap.h>
#elif defined(__APPLE__) && defined(__MACH__)
/* OSX */
#include <libkern/OSByteOrder.h>
#define __bswap_16 OSSwapInt16
#define __bswap_32 OSSwapInt32
#define __bswap_64 OSSwapInt64
#endif

/*
 Python3 compatibility
*/
#if PY_MAJOR_VERSION >= 3
    #define PyCapsule_Type PyCObject_Type
    #define PyString_FromString PyUnicode_FromString
    #define PyString_GET_SIZE PyBytes_GET_SIZE
    #define PyString_AS_STRING PyBytes_AS_STRING
    #define PyInt_AsLong PyLong_AsLong
    #define PyInt_FromLong PyLong_FromLong 
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


/*
  Minimum function for two values
*/

#define min(a,b) (((a) < (b)) ? (a) : (b))


/*
  Exceptions for the Go Fast! Readers
*/

// gofast.cpp
extern PyObject *SyncError;
extern PyObject *EOFError;


/*
  LWA1/LWA-SV Look up tables
*/

// gofast.cpp
extern int16_t tbw4LUT[256][2];
extern int8_t  tbnLUT[256];
extern int8_t  drxLUT[256][2];
extern int8_t  tbfLUT[256][2];


/* 
  Support Functions
*/

// gofast.cpp
extern int validSync5C(uint32_t);
// vdif.cpp
extern void initVDIFLUTs(void);

/*
 ADP COR mode channel information
*/

#define COR_NCHAN 72

/*
  Reader Functions and Documentation
*/

// tbw.cpp
extern PyObject *read_tbw(PyObject*, PyObject*);
extern char read_tbw_doc[];
// tbn.cpp
extern PyObject *read_tbn_cf32(PyObject*, PyObject*);
extern char read_tbn_cf32_doc[];
extern PyObject *read_tbn_ci8(PyObject*, PyObject*);
extern char read_tbn_ci8_doc[];
// drx.cpp
extern PyObject *read_drx_cf32(PyObject*, PyObject*);
extern char read_drx_cf32_doc[];
extern PyObject *read_drx_ci8(PyObject*, PyObject*);
extern char read_drx_ci8_doc[];
// drspec.cpp
extern PyObject *read_drspec(PyObject*, PyObject*);
extern char read_drspec_doc[];

// vdif.cpp
extern PyObject *read_vdif_f32(PyObject*, PyObject*, PyObject*);
extern char read_vdif_f32_doc[];
extern PyObject *read_vdif_i8(PyObject*, PyObject*, PyObject*);
extern char read_vdif_i8_doc[];

// tbf.cpp
extern PyObject *read_tbf_cf32(PyObject*, PyObject*);
extern char read_tbf_cf32_doc[];
extern PyObject *read_tbf_ci8(PyObject*, PyObject*);
extern char read_tbf_ci8_doc[];
// cor.cpp
extern PyObject *read_cor(PyObject*, PyObject*);
extern char read_cor_doc[];

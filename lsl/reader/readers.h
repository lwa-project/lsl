#ifndef __READERS_H
#define __READERS_H

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

// gofast.c
extern PyObject *SyncError;
extern PyObject *EOFError;


/*
  LWA1/LWA-SV Look up tables
*/

// gofast.c
extern short int tbw4LUT[256][2];
extern float tbnLUT[256];
extern float drxLUT[256][2];
extern float tbfLUT[256][2];
extern float ovroLUT[256][2];
extern float drx8LUT[256];


/* 
  Support Functions
*/

// gofast.c
extern int validSync5C(unsigned int);
// vdif.c
extern void initVDIFLUTs(void);

/*
 ADP COR mode channel information
*/

#define COR_NCHAN 72

/*
  Reader Functions and Documentation
*/

// tbw.c
extern PyObject *read_tbw(PyObject*, PyObject*);
extern char read_tbw_doc[];
// tbn.c
extern PyObject *read_tbn(PyObject*, PyObject*);
extern char read_tbn_doc[];
// drx.c
extern PyObject *read_drx(PyObject*, PyObject*);
extern char read_drx_doc[];
// drspec.c
extern PyObject *read_drspec(PyObject*, PyObject*);
extern char read_drspec_doc[];

// vdif.c
extern PyObject *read_vdif(PyObject*, PyObject*, PyObject*);
extern char read_vdif_doc[];

// tbf.c
extern PyObject *read_tbf(PyObject*, PyObject*);
extern char read_tbf_doc[];
// cor.c
extern PyObject *read_cor(PyObject*, PyObject*);
extern char read_cor_doc[];

// ovro.c
extern PyObject *read_ovro_spec(PyObject*, PyObject*);
extern char read_ovro_spec_doc[];

#endif	// __READERS_H

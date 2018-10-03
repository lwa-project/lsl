#include "Python.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#define PY_ARRAY_UNIQUE_SYMBOL gofast_ARRAY_API
#include "numpy/arrayobject.h"

#include "readers.h"


/*
Exceptions for the Go Fast! Readers
*/

PyObject *SyncError;
PyObject *EOFError;


/*
Validate a collection of Mark 5C sync words.  Return 1 is all are
valid.
*/

int validSync5C(unsigned int syncWord) {
    int valid = 0;
    
    if( syncWord == 0x5CDEC0DE ) {
        valid = 1;
    }
    
    return valid;
}


/* 
Look-up Tables
*/
short int tbw4LUT[256][2];
float tbnLUT[256];
float drxLUT[256][2];
float tbfLUT[256][2];

static void initLWALUTs(void) {
    // Look-up table inialization function from the VDIFIO library
    
    int i,j;
    short int t;
    
    //TBW - 4 bit
    for(i=0; i<256; i++) {
        for(j=0; j<2; j++) {
            t = (i >> 4*(1-j)) & 15;
            tbw4LUT[i][j] = t;
            tbw4LUT[i][j] -= ((t&8)<<1);
        }
    }
    
    // TBN
    for(i=0; i<256; i++) {
        tbnLUT[i] = i;
        tbnLUT[i] -= ((i&128)<<1);
    }
    
    // DRX & TBF
    for(i=0; i<256; i++) {
        for(j=0; j<2; j++) {
            t = (i >> 4*(1-j)) & 15;
            drxLUT[i][j] = t;
            drxLUT[i][j] -= ((t&8)<<1);
            
            tbfLUT[i][j] = t;
            tbfLUT[i][j] -= ((t&8)<<1);
        }
    }
}


/*
Module Setup - Function Definitions and Documentation
*/

static PyMethodDef GoFastMethods[] = {
    {"readTBW",    (PyCFunction) readTBW,    METH_VARARGS,               readTBW_doc   }, 
    {"readTBN",    (PyCFunction) readTBN,    METH_VARARGS,               readTBN_doc   }, 
    {"readDRX",    (PyCFunction) readDRX,    METH_VARARGS,               readDRX_doc   }, 
    {"readDRSpec", (PyCFunction) readDRSpec, METH_VARARGS,               readDRSpec_doc},
    {"readVDIF",   (PyCFunction) readVDIF,   METH_VARARGS|METH_KEYWORDS, readVDIF_doc  }, 
    {"readTBF",    (PyCFunction) readTBF,    METH_VARARGS,               readTBF_doc   }, 
    {"readCOR",    (PyCFunction) readCOR,    METH_VARARGS,               readCOR_doc   }, 
    {NULL,         NULL,                     0,                          NULL          }
};

PyDoc_STRVAR(GoFast_doc, "Go Fast! (TM) - TBW, TBN, DRX, DR Spectrometer, and VDIF readers written in C");


/*
Module Setup - Initialization
*/

#if PY_MAJOR_VERSION >= 3
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

MOD_INIT(_gofast) {
    PyObject *m, *dict1, *dict2;
    
    Py_Initialize();
    
    // Initialize the look-up tables
    initLWALUTs();
    initVDIFLUTs();
    
    // Module definitions and functions
    MOD_DEF(m, "_gofast", GoFastMethods, GoFast_doc);
    if( m == NULL ) {
        return MOD_ERROR_VAL;
    }
    import_array();
    
    // Exceptions
    
    //   1.  SyncError -> similar to lsl.reader.errors.SyncError
    dict1 = (PyObject *) PyDict_New();
    if(dict1 == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create exception dictionary");
        Py_XDECREF(dict1);
        Py_XDECREF(m);
        return MOD_ERROR_VAL;
    }
    PyDict_SetItemString(dict1, "__doc__", \
        PyString_FromString("Exception raised when a reader encounters an error with one or more of the four sync. words."));
    SyncError = PyErr_NewException("_gofast.SyncError", PyExc_IOError, dict1);
    Py_INCREF(SyncError);
    PyModule_AddObject(m, "SyncError", SyncError);
    
    //    2. EOFError -> similar to lsl.reader.errors.EOFError
    dict2 = (PyObject *) PyDict_New();
    if(dict2 == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create exception dictionary");
        Py_XDECREF(dict1);
        Py_XDECREF(SyncError);
        Py_XDECREF(dict2);
        Py_XDECREF(m);
        return MOD_ERROR_VAL;
    }
    PyDict_SetItemString(dict2, "__doc__", \
        PyString_FromString("Exception raised when a reader encounters the end-of-file while reading."));
    EOFError = PyErr_NewException("_gofast.EOFError", PyExc_IOError, dict2);
    Py_INCREF(EOFError);
    PyModule_AddObject(m, "EOFError", EOFError);
    
    // Version and revision information
    PyModule_AddObject(m, "__version__", PyString_FromString("0.8"));
    PyModule_AddObject(m, "__revision__", PyString_FromString("$Rev$"));
    
    // Correlator channel count
    PyModule_AddObject(m, "NCHAN_COR", PyInt_FromLong(COR_NCHAN));
    
    return MOD_SUCCESS_VAL(m);
}

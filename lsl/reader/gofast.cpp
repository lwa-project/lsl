#include "Python.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#define PY_ARRAY_UNIQUE_SYMBOL gofast_ARRAY_API
#include "numpy/arrayobject.h"

#include "readers.hpp"


/*
  Exceptions for the Go Fast! Readers
*/

PyObject *SyncError;
PyObject *EOFError;


/*
  Validate a collection of Mark 5C sync words.  Return 1 is all are
  valid.
*/

int validSync5C(uint32_t syncWord) {
    int valid = 0;
    
    if( syncWord == 0x5CDEC0DE ) {
        valid = 1;
    }
    
    return valid;
}


/* 
  Look-up Tables
*/
int16_t tbw4LUT[256][2];
int8_t  tbnLUT[256];
int8_t  drx8LUT[256];
int8_t  drxLUT[256][2];
int8_t  tbfLUT[256][2];

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
    
    // TBN & DRX8
    for(i=0; i<256; i++) {
        tbnLUT[i] = i;
        tbnLUT[i] -= ((i&128)<<1);
        
        drx8LUT[i] = i;
        drx8LUT[i] -= ((i&128)<<1);
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
    {"read_tbw",      (PyCFunction) read_tbw,       METH_VARARGS,               read_tbw_doc      }, 
    {"read_tbn",      (PyCFunction) read_tbn_cf32,  METH_VARARGS,               read_tbn_cf32_doc }, 
    {"read_tbn_ci8",  (PyCFunction) read_tbn_ci8,   METH_VARARGS,               read_tbn_ci8_doc  }, 
    {"read_drx",      (PyCFunction) read_drx_cf32,  METH_VARARGS,               read_drx_cf32_doc }, 
    {"read_drx_ci8",  (PyCFunction) read_drx_ci8,   METH_VARARGS,               read_drx_ci8_doc  }, 
    {"read_drx8",     (PyCFunction) read_drx8_cf32, METH_VARARGS,               read_drx8_cf32_doc}, 
    {"read_drx8_ci8", (PyCFunction) read_drx8_ci8,  METH_VARARGS,               read_drx8_ci8_doc }, 
    {"read_drspec",   (PyCFunction) read_drspec,    METH_VARARGS,               read_drspec_doc   },
    {"read_vdif",     (PyCFunction) read_vdif_f32,  METH_VARARGS|METH_KEYWORDS, read_vdif_f32_doc }, 
    {"read_vdif_i8",  (PyCFunction) read_vdif_i8,   METH_VARARGS|METH_KEYWORDS, read_vdif_i8_doc  }, 
    {"read_tbf",      (PyCFunction) read_tbf_cf32,  METH_VARARGS,               read_tbf_cf32_doc }, 
    {"read_tbf_ci8",  (PyCFunction) read_tbf_ci8,   METH_VARARGS,               read_tbf_ci8_doc  }, 
    {"read_cor",      (PyCFunction) read_cor,       METH_VARARGS,               read_cor_doc      }, 
    {NULL,            NULL,                         0,                          NULL              }
};

PyDoc_STRVAR(GoFast_doc, \
"Go Fast! (TM) - TBW, TBN, DRX, DR Spectrometer, VDIF, TBF, and COR readers\n\
written in C++");


/*
  Module Setup - Initialization
*/

MOD_INIT(_gofast) {
    PyObject *m, *all, *dict1, *dict2;
    
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
    PyModule_AddObject(m, "__version__", PyString_FromString("0.9"));
    
    // Correlator channel count
    PyModule_AddObject(m, "NCHAN_COR", PyInt_FromLong(COR_NCHAN));
    
    // Function listings
    all = PyList_New(0);
    PyList_Append(all, PyString_FromString("read_tbw"));
    PyList_Append(all, PyString_FromString("read_tbn"));
    PyList_Append(all, PyString_FromString("read_tbn_ci8"));
    PyList_Append(all, PyString_FromString("read_drx"));
    PyList_Append(all, PyString_FromString("read_drx_ci8"));
    PyList_Append(all, PyString_FromString("read_drx8"));
    PyList_Append(all, PyString_FromString("read_drx8_ci8"));
    PyList_Append(all, PyString_FromString("read_drspec"));
    PyList_Append(all, PyString_FromString("read_vdif"));
    PyList_Append(all, PyString_FromString("read_tbf"));
    PyList_Append(all, PyString_FromString("read_tbf_ci8"));
    PyList_Append(all, PyString_FromString("read_cor"));
    PyList_Append(all, PyString_FromString("SyncError"));
    PyList_Append(all, PyString_FromString("EOFError"));
    PyList_Append(all, PyString_FromString("NCHAN_COR"));
    PyModule_AddObject(m, "__all__", all);
    
    return MOD_SUCCESS_VAL(m);
}

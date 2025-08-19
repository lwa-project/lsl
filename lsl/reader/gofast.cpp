#include "Python.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <complex>

#define PY_ARRAY_UNIQUE_SYMBOL gofast_ARRAY_API
#include "numpy/arrayobject.h"

#include "readers.hpp"

/*
  Exceptions for the Go Fast! Readers
*/

PyObject *SyncError;
PyObject *EOFError;


/* 
  Look-up Tables
*/
int8_t  drx8LUT[256];
int8_t  drxLUT[256][2];
int8_t  tbxLUT[256][2];

static int luts_loaded = 0;

static void initLWALUTs(void) {
    // Look-up table inialization function from the VDIFIO library
    
    int i,j;
    short int t;
    
    // DRX8
    for(i=0; i<256; i++) {
        drx8LUT[i] = i;
        drx8LUT[i] -= ((i&128)<<1);
    }
    
    // DRX & TBX
    for(i=0; i<256; i++) {
        for(j=0; j<2; j++) {
            t = (i >> 4*(1-j)) & 15;
            drxLUT[i][j] = t;
            drxLUT[i][j] -= ((t&8)<<1);
            
            tbxLUT[i][j] = t;
            tbxLUT[i][j] -= ((t&8)<<1);
        }
    }
}


/*
  Module Setup - Function Definitions and Documentation
*/

static PyMethodDef gofast_methods[] = {
    {"read_drx",      (PyCFunction) read_drx_cf32,  METH_VARARGS,               read_drx_cf32_doc }, 
    {"read_drx_ci8",  (PyCFunction) read_drx_ci8,   METH_VARARGS,               read_drx_ci8_doc  }, 
    {"read_drx8",     (PyCFunction) read_drx8_cf32, METH_VARARGS,               read_drx8_cf32_doc}, 
    {"read_drx8_ci8", (PyCFunction) read_drx8_ci8,  METH_VARARGS,               read_drx8_ci8_doc }, 
    {"read_drspec",   (PyCFunction) read_drspec,    METH_VARARGS,               read_drspec_doc   },
    {"read_vdif",     (PyCFunction) read_vdif_f32,  METH_VARARGS|METH_KEYWORDS, read_vdif_f32_doc }, 
    {"read_vdif_i8",  (PyCFunction) read_vdif_i8,   METH_VARARGS|METH_KEYWORDS, read_vdif_i8_doc  }, 
    {"read_cor",      (PyCFunction) read_cor,       METH_VARARGS,               read_cor_doc      },
    {"read_tbx",      (PyCFunction) read_tbx_cf32,  METH_VARARGS,               read_tbx_cf32_doc }, 
    {"read_tbx_ci8",  (PyCFunction) read_tbx_ci8,   METH_VARARGS,               read_tbx_ci8_doc  }, 
    {"read_cor",      (PyCFunction) read_cor,       METH_VARARGS,               read_cor_doc      },
    {NULL,            NULL,                         0,                          NULL              }
};

PyDoc_STRVAR(gofast_doc, \
"Go Fast! (TM) - DRX, DR Spectrometer, VDIF, COR, and TBX readers\n\
written in C++");


/*
  Module Setup - Initialization
*/

static int gofast_exec(PyObject *module) {
    import_array();
    
    if( !luts_loaded ) {
        // Initialize the look-up tables
        initLWALUTs();
        initVDIFLUTs();
        luts_loaded = 1;
    } else {
        PyErr_SetString(PyExc_ImportError, "cannot load module more than once per process");
        return -1;
    }
    
    // Exceptions
    //   1.  SyncError -> similar to lsl.reader.errors.SyncError
    PyObject* dict1 = PyDict_New();
    if(dict1 == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create exception dictionary");
        Py_XDECREF(dict1);
        return -1;
    }
    PyDict_SetItemString(dict1, "__doc__", \
        PyUnicode_FromString("Exception raised when a reader encounters an error with one or more of the four sync. words."));
    SyncError = PyErr_NewException("_gofast.SyncError", PyExc_IOError, dict1);
    Py_INCREF(SyncError);
    PyModule_AddObject(module, "SyncError", SyncError);
    
    //    2. EOFError -> similar to lsl.reader.errors.EOFError
    PyObject* dict2 = PyDict_New();
    if(dict2 == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create exception dictionary");
        Py_XDECREF(dict1);
        Py_XDECREF(SyncError);
        Py_XDECREF(dict2);
        return -1;
    }
    PyDict_SetItemString(dict2, "__doc__", \
        PyUnicode_FromString("Exception raised when a reader encounters the end-of-file while reading."));
    EOFError = PyErr_NewException("_gofast.EOFError", PyExc_IOError, dict2);
    Py_INCREF(EOFError);
    PyModule_AddObject(module, "EOFError", EOFError);
    
    // Version and revision information
    PyModule_AddObject(module, "__version__", PyUnicode_FromString("1.1"));
    
    // Function listings
    PyObject* all = PyList_New(0);
    PyList_Append(all, PyUnicode_FromString("read_drx"));
    PyList_Append(all, PyUnicode_FromString("read_drx_ci8"));
    PyList_Append(all, PyUnicode_FromString("read_drx8"));
    PyList_Append(all, PyUnicode_FromString("read_drx8_ci8"));
    PyList_Append(all, PyUnicode_FromString("read_drspec"));
    PyList_Append(all, PyUnicode_FromString("read_vdif"));
    PyList_Append(all, PyUnicode_FromString("read_cor"));
    PyList_Append(all, PyUnicode_FromString("read_tbx"));
    PyList_Append(all, PyUnicode_FromString("read_tbx_ci8"));
    PyList_Append(all, PyUnicode_FromString("SyncError"));
    PyList_Append(all, PyUnicode_FromString("EOFError"));
    PyModule_AddObject(module, "__all__", all);
    return 0;
}

static PyModuleDef_Slot gofast_slots[] = {
    {Py_mod_exec, (void *)&gofast_exec},
    {0,           NULL}
};

static PyModuleDef gofast_def = {
    PyModuleDef_HEAD_INIT,    /* m_base */
    "_gofast",                /* m_name */
    gofast_doc,               /* m_doc */
    0,                        /* m_size */
    gofast_methods,           /* m_methods */
    gofast_slots,             /* m_slots */
    NULL,                     /* m_traverse */
    NULL,                     /* m_clear */
    NULL,                     /* m_free */
};

PyMODINIT_FUNC PyInit__gofast(void) {
    return PyModuleDef_Init(&gofast_def);
}

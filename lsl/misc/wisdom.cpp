#include "Python.h"
#include <cmath>
#include <complex>
#include <iostream>
#include <fftw3.h>

#ifdef _OPENMP
    #include <omp.h>
    
    // OpenMP scheduling method
    #ifndef OMP_SCHEDULER
    #define OMP_SCHEDULER dynamic
	#endif
#endif

#include "numpy/arrayobject.h"
#include "numpy/npy_math.h"

#include "../common/py3_compat.h"
#include "../correlator/common.hpp"


#define MAXTRANSFORM 262144


/*
 buildWisdom - Build up LSL-specific FFTW wisdom.  Based on makewisdom.c 
               included with PRESTO.
*/

static PyObject *buildWisdom(PyObject *self, PyObject *args) {
    PyObject *output;
    int fftlen;
    FILE *fh;
    fftw_plan_t plan;
    complex_t *inout;
    real_t *inR;
    char *filename;
    
    if(!PyArg_ParseTuple(args, "s", &filename)) {
        PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
        return NULL;
    }
    
    Py_BEGIN_ALLOW_THREADS
    
    if(!FFTW_IMPORT_SYSTEM_WISDOM()) {
        std::cout << "Warning: system wisdom file not found, continuing" << std::endl;
    }
    
    // Real to complex - powers of 2
    fftlen = 2;
    while(fftlen <= MAXTRANSFORM) {
        // Setup
        inR = (real_t*) FFTW_MALLOC(sizeof(real_t) * 2*fftlen);
        inout = (complex_t*) FFTW_MALLOC(sizeof(complex_t) * (fftlen+1));
        
        // Forward
        plan = FFTW_PLAN_DFT_R2C_1D(2*fftlen, \
                                    inR, \
                                    reinterpret_cast<fftw_complex_t*>(inout), \
                                    FFTW_PATIENT);
        FFTW_DESTROY_PLAN(plan);
        
        // Teardown
        FFTW_FREE(inR);
        FFTW_FREE(inout);
        
        // Next
        fftlen *= 2;
    }
    
    // Real to complex - powers of 10
    fftlen = 10;
    while(fftlen <= MAXTRANSFORM) {
        // Setup
        inR = (real_t*) FFTW_MALLOC(sizeof(real_t) * 2*fftlen);
        inout = (complex_t*) FFTW_MALLOC(sizeof(complex_t) * (fftlen+1));
        
        // Forward
        plan = FFTW_PLAN_DFT_R2C_1D(2*fftlen, \
                                    inR, \
                                    reinterpret_cast<fftw_complex_t*>(inout), \
                                    FFTW_PATIENT);
        FFTW_DESTROY_PLAN(plan);
        
        // Teardown
        FFTW_FREE(inR);
        FFTW_FREE(inout);
        
        // Next
        fftlen *= 10;
    }
    
    // Complex in-place - powers of 2
    fftlen = 2;
    while(fftlen <= MAXTRANSFORM) {
        // Setup
        inout = (complex_t*) FFTW_MALLOC(sizeof(complex_t) * fftlen);
        
        // Forward
        plan = FFTW_PLAN_DFT_1D(fftlen, \
                                reinterpret_cast<fftw_complex_t*>(inout), \
                                reinterpret_cast<fftw_complex_t*>(inout), \
                                FFTW_FORWARD, FFTW_PATIENT);
        FFTW_DESTROY_PLAN(plan);
        
        // Backward
        plan = FFTW_PLAN_DFT_1D(fftlen, \
                                reinterpret_cast<fftw_complex_t*>(inout), \
                                reinterpret_cast<fftw_complex_t*>(inout), \
                                FFTW_BACKWARD, FFTW_PATIENT);
        FFTW_DESTROY_PLAN(plan);
        
        // Teardown
        FFTW_FREE(inout);
        
        // Next
        fftlen *= 2;
    }
    
    // Complex in-place - powers of 10
    fftlen = 10;
    while(fftlen <= MAXTRANSFORM) {
        // Setup
        inout = (complex_t*) FFTW_MALLOC(sizeof(complex_t) * fftlen);
        
        // Forward
        plan = FFTW_PLAN_DFT_1D(fftlen, \
                                reinterpret_cast<fftw_complex_t*>(inout), \
                                reinterpret_cast<fftw_complex_t*>(inout), \
                                FFTW_FORWARD, FFTW_PATIENT);
        FFTW_DESTROY_PLAN(plan);
        
        // Backward
        plan = FFTW_PLAN_DFT_1D(fftlen, \
                                reinterpret_cast<fftw_complex_t*>(inout), \
                                reinterpret_cast<fftw_complex_t*>(inout), \
                                FFTW_BACKWARD, FFTW_PATIENT);
        FFTW_DESTROY_PLAN(plan);
        
        // Teardown
        FFTW_FREE(inout);
        
        // Next
        fftlen *= 10;
    }
    
    // Save the wisdom
    fh = fopen(filename, "w");
    if(ferror(fh)) {
        PyErr_Format(PyExc_IOError, "An error occurred while opening the file for writing");
        return NULL;
    }
    FFTW_EXPORT_WISDOM_TO_FILE(fh);
    fclose(fh);
    
    Py_END_ALLOW_THREADS
    
    output = Py_True;
    Py_INCREF(output);
    return output;
}

PyDoc_STRVAR(buildWisdom_doc, \
"Build a LSL-specific wisdom file.\n\
\n\
Input arguments are:\n\
 * filename: filename to write the wisdom to\n\
\n\
Outputs are:\n\
 True if the wisdom was built successfully\n\
");


static PyMethodDef WisdomMethods[] = {
    {"buildWisdom", (PyCFunction) buildWisdom, METH_VARARGS, buildWisdom_doc }, 
    {NULL,          NULL,                      0,            NULL            }
};

PyDoc_STRVAR(wisdom_doc, \
"Extension to build LSL-specific FFTW wisdom.\n\
\n\
The functions defined in this module are:\n\
 * buildWisdom - Build aLSL wisdom file\n\
\n\
See the individual functions for more details.");


/*
  Module Setup - Initialization
*/

MOD_INIT(_wisdom) {
    PyObject *m;
    
    Py_Initialize();
    
    // Module definitions and functions
    MOD_DEF(m, "_wisdom", WisdomMethods, wisdom_doc);
    if( m == NULL ) {
        return MOD_ERROR_VAL;
    }
    import_array();
    
    // Version and revision information
    PyModule_AddObject(m, "__version__", PyString_FromString("0.4"));
    
    return MOD_SUCCESS_VAL(m);
}

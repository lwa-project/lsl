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
#include "../correlator/common.h"


#define MAXTRANSFORM 262144


/*
 buildWisdom - Build up LSL-specific FFTW wisdom.  Based on makewisdom.c 
               included with PRESTO.
*/

static PyObject *buildWisdom(PyObject *self, PyObject *args) {
    PyObject *output;
    int fftlen;
    FILE *fh;
    fftwf_plan plan;
    Complex32 *inout;
    float *inR;
    char *filename;
    
    if(!PyArg_ParseTuple(args, "s", &filename)) {
        PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
        return NULL;
    }
    
    Py_BEGIN_ALLOW_THREADS
    
    if(!fftwf_import_system_wisdom()) {
        std::cout << "Warning: system wisdom file not found, continuing" << std::endl;
    }
    
    // Real to complex - powers of 2
    fftlen = 2;
    while(fftlen <= MAXTRANSFORM) {
        // Setup
        inR = (float*) fftwf_malloc(sizeof(float) * 2*fftlen);
        inout = (Complex32*) fftwf_malloc(sizeof(Complex32) * (fftlen+1));
        
        // Forward
        plan = fftwf_plan_dft_r2c_1d(2*fftlen, \
                                     inR, \
                                     reinterpret_cast<fftwf_complex*>(inout), \
                                     FFTW_PATIENT);
        fftwf_destroy_plan(plan);
        
        // Teardown
        fftwf_free(inR);
        fftwf_free(inout);
        
        // Next
        fftlen *= 2;
    }
    
    // Real to complex - powers of 10
    fftlen = 10;
    while(fftlen <= MAXTRANSFORM) {
        // Setup
        inR = (float*) fftwf_malloc(sizeof(float) * 2*fftlen);
        inout = (Complex32*) fftwf_malloc(sizeof(Complex32) * (fftlen+1));
        
        // Forward
        plan = fftwf_plan_dft_r2c_1d(2*fftlen, \
                                     inR, \
                                     reinterpret_cast<fftwf_complex*>(inout), \
                                     FFTW_PATIENT);
        fftwf_destroy_plan(plan);
        
        // Teardown
        fftwf_free(inR);
        fftwf_free(inout);
        
        // Next
        fftlen *= 10;
    }
    
    // Complex in-place - powers of 2
    fftlen = 2;
    while(fftlen <= MAXTRANSFORM) {
        // Setup
        inout = (Complex32*) fftwf_malloc(sizeof(Complex32) * fftlen);
        
        // Forward
        plan = fftwf_plan_dft_1d(fftlen, \
                                 reinterpret_cast<fftwf_complex*>(inout), \
                                 reinterpret_cast<fftwf_complex*>(inout), \
                                 FFTW_FORWARD, FFTW_PATIENT);
        fftwf_destroy_plan(plan);
        
        // Backward
        plan = fftwf_plan_dft_1d(fftlen, \
                                 reinterpret_cast<fftwf_complex*>(inout), \
                                 reinterpret_cast<fftwf_complex*>(inout), \
                                 FFTW_BACKWARD, FFTW_PATIENT);
        fftwf_destroy_plan(plan);
        
        // Teardown
        fftwf_free(inout);
        
        // Next
        fftlen *= 2;
    }
    
    // Complex in-place - powers of 10
    fftlen = 10;
    while(fftlen <= MAXTRANSFORM) {
        // Setup
        inout = (Complex32*) fftwf_malloc(sizeof(Complex32) * fftlen);
        
        // Forward
        plan = fftwf_plan_dft_1d(fftlen, \
                                 reinterpret_cast<fftwf_complex*>(inout), \
                                 reinterpret_cast<fftwf_complex*>(inout), \
                                 FFTW_FORWARD, FFTW_PATIENT);
        fftwf_destroy_plan(plan);
        
        // Backward
        plan = fftwf_plan_dft_1d(fftlen, \
                                 reinterpret_cast<fftwf_complex*>(inout), \
                                 reinterpret_cast<fftwf_complex*>(inout), \
                                 FFTW_BACKWARD, FFTW_PATIENT);
        fftwf_destroy_plan(plan);
        
        // Teardown
        fftwf_free(inout);
        
        // Next
        fftlen *= 10;
    }
    
    // Save the wisdom
    fh = fopen(filename, "w");
    if(ferror(fh)) {
        PyErr_Format(PyExc_IOError, "An error occurred while opening the file for writing");
        return NULL;
    }
    fftwf_export_wisdom_to_file(fh);
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

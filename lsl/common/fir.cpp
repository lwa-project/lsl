#include "Python.h"
#include <cmath>
#include <limits>
#include <complex>
#include <fftw3.h>
#include <cstdlib>

#ifdef _OPENMP
    #include <omp.h>
    
    // OpenMP scheduling method
    #ifndef OMP_SCHEDULER
    #define OMP_SCHEDULER guided
    #endif
#endif

#include "numpy/arrayobject.h"
#include "numpy/npy_math.h"

#include "../correlator/common.hpp"
#include "../correlator/pool.hpp"

/*
   64-bit aligned memory pool
*/

static Aligned64BufferPool& apool = get_aligned64_buffer_pool("fir");


/*
applyFIR - Given a pointer to a 16-bit integer data stream, the number of data samples to 
            process, a pointer a set of 16-bit FIR coefficients, the number of taps, and a
            pointer to a 32-bit float output stream, apply the FIR coefficents.
*/

template<typename InType, typename OutType>
void applyFIR(long nSamps,
              long nTaps,
              InType const *data, 
              InType const *coeff, 
              OutType *output) {
    long i, j;
    
    memset(output, 0, sizeof(OutType)*nSamps);
    
    for(i=0; i<nTaps; i++) {
        for(j=0; j<=i; j++) {
            *(output + i) += (OutType) *(coeff + j) * *(data + i - j);
        }
        *(output + i) /= std::numeric_limits<InType>::max();
    }
    
    for(i=nTaps; i<nSamps; i++) {
        for(j=0; j<nTaps; j++) {
            *(output + i) += (OutType) *(coeff + j) * *(data + i - j);
        }
        *(output + i) /= std::numeric_limits<InType>::max();
    }
}


/*
applyFIRDelayed - Like applyFIR but also applies a sample delay before filtering.
*/

template<typename InType, typename OutType>
void applyFIRDelayed(long nSamps,
                     long nTaps,
                     long sampleDelay,
                     InType const *data, 
                     InType const *coeff, 
                     OutType *output) {
    long i, j;
    
    memset(output, 0, sizeof(OutType)*nSamps);
    
    for(i=0; i<sampleDelay; i++) {
        *(output + i) = 0;
    }
    
    for(i=sampleDelay; i<(nTaps+sampleDelay); i++) {
        for(j=0; j<=(i-sampleDelay); j++) {
            *(output + i) += (OutType) *(coeff + j) * *(data + i - sampleDelay - j);
        }
        *(output + i) /= std::numeric_limits<InType>::max();
    }
    
    for(i=(nTaps+sampleDelay); i<nSamps; i++) {
        for(j=0; j<nTaps; j++) {
            *(output + i) += (OutType) *(coeff + j) * *(data + i - sampleDelay - j);
        }
        *(output + i) /= std::numeric_limits<InType>::max();
    }
}


/*
integerFIR - Function for taking a numpy.int16 data set and applying a FIR
            filter.  This function is similar to scipy.signal.lfilter but
            is more memory efficient since it doesn't require floats.
*/

static PyObject *integerFIR(PyObject *self, PyObject *args, PyObject *kwds) {
    PyObject *signals, *filter, *output;
    PyArrayObject *data=NULL, *coeff=NULL, *dataF=NULL;
    
    long nSamps, nTaps;
    
    if(!PyArg_ParseTuple(args, "OO", &signals, &filter)) {
        PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
        goto fail;
    }

    // Bring the data into C and make it usable
    data  = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_INT16, 1, 1);
    coeff = (PyArrayObject *) PyArray_ContiguousFromObject(filter, NPY_INT16, 1, 1);
    if( data == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input signals array to 1-D int16");
        goto fail;
    }
    if( coeff == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input filter array to 1-D int16");
        goto fail;
    }
    
    // Get sample and tap counts
    nSamps = (long) PyArray_DIM(data, 0);
    nTaps  = (long) PyArray_DIM(coeff, 0);
    
    // Create the output data holders
    npy_intp dims[1];
    dims[0] = (npy_intp) nSamps;
    dataF = (PyArrayObject*) PyArray_ZEROS(1, dims, NPY_INT32, 0);
    if(dataF == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array");
        goto fail;
    }
    
    Py_BEGIN_ALLOW_THREADS
    
    // Go
    short int *a, *b;
    float *c;
    a = (short int *) PyArray_DATA(data);
    b = (short int *) PyArray_DATA(coeff);
    c = (float *) PyArray_DATA(dataF);
    applyFIR(nSamps, nTaps, a, b, c);
    
    Py_END_ALLOW_THREADS
    
    output = Py_BuildValue("O", PyArray_Return(dataF));
    
    Py_XDECREF(data);
    Py_XDECREF(coeff);
    Py_XDECREF(dataF);
    
    return output;
    
fail:
    Py_XDECREF(data);
    Py_XDECREF(coeff);
    Py_XDECREF(dataF);
    
    return NULL;
}

PyDoc_STRVAR(integerFIR_doc, \
"Given a 1-D numpy.int16 array of data values and a numpy.int16 array of FIR\n\
coefficients, apply the coefficients to the data.\n\
\n\
Inputs arguments are:\n\
 * data: 1-D numpy.int16 array of data\n\
 * coeffs: 1-D numpy.int16 array of FIR coefficients\n\
\n\
Outputs:\n\
 * result: 1-D numpy.float32 array of the filtered data\n\
");


/*
integerFIRDelayed - Function similar to integerFIR but adds in a FIFO delay segment.
*/

static PyObject *integerFIRDelayed(PyObject *self, PyObject *args, PyObject *kwds) {
    PyObject *signals, *filter, *output;
    PyArrayObject *data=NULL, *coeff=NULL, *dataF=NULL;
    
    long nSamps, nTaps, sampleDelay;
    
    if(!PyArg_ParseTuple(args, "OOl", &signals, &filter, &sampleDelay)) {
        PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
        goto fail;
    }

    // Bring the data into C and make it usable
    data  = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_INT16, 1, 1);
    coeff = (PyArrayObject *) PyArray_ContiguousFromObject(filter, NPY_INT16, 1, 1);
    if( data == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input signals array to 1-D int16");
        goto fail;
    }
    if( coeff == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input filter array to 1-D int16");
        goto fail;
    }
    
    // Get sample and tap counts
    nSamps = (long) PyArray_DIM(data, 0);
    nTaps  = (long) PyArray_DIM(coeff, 0);
    
    // Create the output data holders
    npy_intp dims[1];
    dims[0] = (npy_intp) nSamps;
    dataF = (PyArrayObject*) PyArray_ZEROS(1, dims, NPY_FLOAT32, 0);
    if(dataF == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array");
        goto fail;
    }
    
    Py_BEGIN_ALLOW_THREADS
    
    // Go
    short int *a, *b;
    float *c;
    a = (short int *) PyArray_DATA(data);
    b = (short int *) PyArray_DATA(coeff);
    c = (float *) PyArray_DATA(dataF);
    applyFIRDelayed(nSamps, nTaps, sampleDelay, a, b, c);
    
    Py_END_ALLOW_THREADS
    
    output = Py_BuildValue("O", PyArray_Return(dataF));
    
    Py_XDECREF(data);
    Py_XDECREF(coeff);
    Py_XDECREF(dataF);
    
    return output;
    
fail:
    Py_XDECREF(data);
    Py_XDECREF(coeff);
    Py_XDECREF(dataF);
    
    return NULL;
}

PyDoc_STRVAR(integerFIRDelayed_doc, \
"Given a 1-D numpy.int16 array of data values, a numpy.int16 array of FIR\n\
coefficients, and a FIFO sample delay, delay the signal and apply the\n\
coefficients to the data.\n\
\n\
Inputs arguments are:\n\
 * data: 1-D numpy.int16 array of data\n\
 * coeffs: 1-D numpy.int16 array of FIR coefficients\n\
 * sampleDelay: interger number of samples to delay the signal (must be >=0)\n\
\n\
Outputs:\n\
 * result: 1-D numpy.float32 array of the delayed and filtered data\n\
");


static PyObject *integerBeamformer(PyObject *self, PyObject *args, PyObject *kwds) {
    PyObject *signals, *filters, *courses, *fines, *gains, *output;
    PyArrayObject *data=NULL, *filter=NULL, *course=NULL, *fine=NULL, *gain=NULL, *dataFX=NULL, *dataFY=NULL;
    
    long i, j, k, nStand, nSamps, nFilts, nTaps, maxCourse;
    
    if(!PyArg_ParseTuple(args, "OOOOO", &signals, &filters, &courses, &fines, &gains)) {
        PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
        return NULL;
    }

    // Bring the data into C and make it usable
    data   = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_INT16, 2, 2);
    filter = (PyArrayObject *) PyArray_ContiguousFromObject(filters, NPY_INT16, 3, 3);
    course = (PyArrayObject *) PyArray_ContiguousFromObject(courses, NPY_INT16, 1, 1);
    fine   = (PyArrayObject *) PyArray_ContiguousFromObject(fines,   NPY_INT16, 1, 1);
    gain   = (PyArrayObject *) PyArray_ContiguousFromObject(gains,   NPY_INT16, 2, 2);
    if( data == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input signals array to 2-D int16");
        goto fail;
    }
    if( filter == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input filters array to 3-D int16");
        goto fail;
    }
    if( course == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input courses array to 1-D int16");
        goto fail;
    }
    if( fine == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input fines array to 1-D int16");
        goto fail;
    }
    if( gain == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input gains array to 2-D int16");
        goto fail;
    }
    
    // Check data dimensions
    if( PyArray_DIM(data, 0) != PyArray_DIM(filter, 0) ) {
        PyErr_Format(PyExc_RuntimeError, "signals and FIR filters have different input counts");
        goto fail;
    }
    if( PyArray_DIM(data, 0) != PyArray_DIM(course, 0) ) {
        PyErr_Format(PyExc_RuntimeError, "signals and course delays have different input counts");
        goto fail;
    }
    if( PyArray_DIM(data, 0) != PyArray_DIM(fine, 0) ) {
        PyErr_Format(PyExc_RuntimeError, "signals and find delays have different input counts");
        goto fail;
    }
    if( PyArray_DIM(data, 0)/2 != PyArray_DIM(gain, 0) ) {
        PyErr_Format(PyExc_RuntimeError, "signals and gains have different input counts");
        goto fail;
    }
    if( PyArray_DIM(gain, 1) != 4 ) {
        PyErr_Format(PyExc_RuntimeError, "seconds dimension of gains must be four");
        goto fail;
    }
    
    // Get sample and tap counts
    nStand = (long) PyArray_DIM(data, 0);
    nSamps = (long) PyArray_DIM(data, 1);
    nFilts = (long) PyArray_DIM(filter, 1);
    nTaps  = (long) PyArray_DIM(filter, 2);
    
    maxCourse = 0;
    short int *c;
    c = (short int *) PyArray_DATA(course);
    for(i=0; i<nStand/2; i++) {
        k = 2*i;
        if( *(c + k) > maxCourse ) {
            maxCourse = (long) *(c + k);
        }
        k = 2*i + 1;
        if( *(c + k) > maxCourse ) {
            maxCourse = (long) *(c + k);
        }
    }
    maxCourse += nTaps/2 - 1;
    
    // Create the output data holders
    npy_intp dims[1];
    dims[0] = (npy_intp) (nSamps - maxCourse);
    dataFX = (PyArrayObject*) PyArray_ZEROS(1, dims, NPY_FLOAT32, 0);
    if(dataFX == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array for X polarization");
        goto fail;
    }
    
    dataFY = (PyArrayObject*) PyArray_ZEROS(1, dims, NPY_FLOAT32, 0);
    if(dataFY == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array for Y polarization");
        goto fail;
    }
    
    Py_BEGIN_ALLOW_THREADS
    
    // Go
    short int *d, *f, *w, *g;
    float *x, *y;
    d = (short int *) PyArray_DATA(data);
    f = (short int *) PyArray_DATA(filter);
    c = (short int *) PyArray_DATA(course);
    w = (short int *) PyArray_DATA(fine);
    g = (short int *) PyArray_DATA(gain);
    x = (float *) PyArray_DATA(dataFX);
    y = (float *) PyArray_DATA(dataFY);
    
    float *t1, *t2, *tX, *tY;
    
    #ifdef _OPENMP
        #pragma omp parallel default(shared) private(j, k, t1, t2, tX, tY)
    #endif
    {
        #ifdef _OPENMP
            #pragma omp for schedule(OMP_SCHEDULER)
        #endif
        for(i=0; i<nStand/2; i++) {
            /* 
            Skip stands with all gains of zero.
            */
            if( *(g + 4*i + 0) == 0 && *(g + 4*i + 1) == 0 && *(g + 4*i + 2) == 0 && *(g + 4*i + 3) == 0 ) {
                continue;
            }
            
            /*
            Allocate temporary arrays
                * t1 -> output of integer delay and FIR filter for X pol.
                * t2 -> output of integer delay and FIR filter for Y pol.
                * tX -> temporary beam for X pol.
                * tY -> temporary beam for Y pol.
            */
            t1 = (float *) apool.acquire<float>(nSamps);
            t2 = (float *) apool.acquire<float>(nSamps);
            tX = (float *) apool.acquire<float>(nSamps);
            tY = (float *) apool.acquire<float>(nSamps);
            
            /*
            Delay + FIR
            */
            k = 2*i;
            applyFIRDelayed(nSamps, nTaps, *(c+k), (d+k*nSamps), (f + k*nFilts*nTaps + *(w+k)*nTaps), t1);
            k = 2*i + 1;
            applyFIRDelayed(nSamps, nTaps, *(c+k), (d+k*nSamps), (f + k*nFilts*nTaps + *(w+k)*nTaps), t2);
            
            for(j=0; j<nSamps; j++) {
                *(tX + j) = *(t1 + j) * *(g + 4*i + 0) / 32767.0;
                *(tY + j) = *(t1 + j) * *(g + 4*i + 1) / 32767.0;
                
                *(tX + j) += *(t2 + j) * *(g + 4*i + 2) / 32767.0;
                *(tY + j) += *(t2 + j) * *(g + 4*i + 3) / 32767.0;
            }
            
            /*
            Add the partial sums to the full beam
            */
            #ifdef _OPENMP
                #pragma omp critical
            #endif
            {
                for(j=0; j<(nSamps-maxCourse); j++) {
                    *(x + j) += *(tX + j + maxCourse);
                    *(y + j) += *(tY + j + maxCourse);
                }
            }
            
            /*
            Cleanup
            */
//             free(t1);
//             free(t2);
//             free(tX);
//             free(tY);
        }
    }
    
    Py_END_ALLOW_THREADS
    
    output = Py_BuildValue("(OO)", PyArray_Return(dataFX), PyArray_Return(dataFY));
    
    Py_XDECREF(data);
    Py_XDECREF(filter);
    Py_XDECREF(course);
    Py_XDECREF(fine);
    Py_XDECREF(gain);
    Py_XDECREF(dataFX);
    Py_XDECREF(dataFY);
    
    return output;
    
fail:
    Py_XDECREF(data);
    Py_XDECREF(filter);
    Py_XDECREF(course);
    Py_XDECREF(fine);
    Py_XDECREF(gain);
    Py_XDECREF(dataFX);
    Py_XDECREF(dataFY);
    
    return NULL;
}

PyDoc_STRVAR(integerBeamformer_doc, \
"Given a 2-D numpy.int16 array (stands by samples) of data values, 3-D array of FIR\n\
filter coefficients (stands by filters by taps), a 1-D numpy.int16 array of course\n\
(FIFO) delays, a 1-D numpy.int16 array of fine delay FIR filters, and a 2-D array of\n\
gains (stands by [XX, XY, YX, YY]), apply the delays and sum the signals.\n\
\n\
Inputs arguments are:\n\
 * data: 2-D numpy.int16 array of data (stands by samples)\n\
 * coeffs: 3-D numpy.int16 array of FIR coefficients (stands by filters by taps)\n\
 * course: 1-D numpy.int16 array of FIFO delays in samples\n\
 * fine: 1-D numpy.int16 array of which FIR filter to apply for fine delay\n\
 * gain: 2-D numpy.int16 arry of gains (stands by [XX, XY, YX, YY]), where XX is the X\n\
         contribution to the output X pol., XY is the X contribution to the output Y\n\
         pol., YX is the Y contribtion to the output X pol., and YY is the Y contribu-\n\
         tion to the output Y pol.\n\
\n\
Outputs:\n\
 * results: two element tuple (output X, outpuY) of the beamformer sum.  Each element is\n\
            a 1-D numpy.float32 array.\n\
\n\
.. note::\n\
\tThe structure of data is assumed to be that the polarizations are ordered, e.g., the X\n\
\tpolarization of stand 1 is immediately followed by the Y polarization.\n\
");


static PyMethodDef fir_methods[] = {
    {"integer16",         (PyCFunction) integerFIR,        METH_VARARGS, integerFIR_doc       }, 
    {"integer16Delayed",  (PyCFunction) integerFIRDelayed, METH_VARARGS, integerFIRDelayed_doc}, 
    {"integerBeamformer", (PyCFunction) integerBeamformer, METH_VARARGS, integerBeamformer_doc}, 
    {NULL,                NULL,              0,            NULL                 }
};

PyDoc_STRVAR(fir_doc, \
"This module contains a collection of function to speed up FIR filtering of TBW\n\
data (represented as numpy.int16 arrays) and the SoftwareDP class.  The funtions\n\
provided in this module are:\n\
 * integer16: Apply a FIR filter to numpy.int16 data,\n\
 * integer16Delayed: Apply a FIFO delay and a FIR filter to numpy.int16 data, and\n\
 * integerBeamformer: Software implementation of the DP beamformer.\n\
");


/*
Module Setup - Initialization
*/

static int fir_exec(PyObject *module) {
    import_array();
    
    // Version and revision information
    PyModule_AddObject(module, "__version__", PyUnicode_FromString("0.2"));
    
    // Function listings
    PyObject *all = PyList_New(0);
    PyList_Append(all, PyUnicode_FromString("integer16"));
    PyList_Append(all, PyUnicode_FromString("integer16Delayed"));
    PyList_Append(all, PyUnicode_FromString("integerBeamformer"));
    PyModule_AddObject(module, "__all__", all);
    return 0;
}

static PyModuleDef_Slot fir_slots[] = {
    {Py_mod_exec, (void *)&fir_exec},
    {0,           NULL}
};

static PyModuleDef fir_def = {
    PyModuleDef_HEAD_INIT,    /* m_base */
    "_fir",                   /* m_name */
    fir_doc,                  /* m_doc */
    0,                        /* m_size */
    fir_methods,              /* m_methods */
    fir_slots,                /* m_slots */
    NULL,                     /* m_traverse */
    NULL,                     /* m_clear */
    NULL,                     /* m_free */
};

PyMODINIT_FUNC PyInit__fir(void) {
    return PyModuleDef_Init(&fir_def);
}

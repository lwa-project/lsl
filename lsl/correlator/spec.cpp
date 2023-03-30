#include "Python.h"
#include <cmath>
#include <complex>
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
#include "common.h"
#include "blas.h"

/*
  Holder for window function callback
*/

static PyObject *windowFunc = NULL;


/*
  FFT Functions
    1. FPSD - FFT a real or complex-valued collection of signals
*/

template<typename InType, typename OutType>
void compute_psd_real(long nStand,
                      long nSamps,
                      long nFFT,
                      int nChan, 
                      int Overlap, 
                      int Clip, 
                      InType const* data,
                      double const* window,
                      OutType* psd) {
    // Setup
    long i, j, k;
    
    Py_BEGIN_ALLOW_THREADS
    
    // Create the FFTW plan                          
    float *inP, *in;                          
    Complex32 *outP, *out;
    inP = (float*) fftwf_malloc(sizeof(float) * 2*nChan);
    outP = (Complex32*) fftwf_malloc(sizeof(Complex32) * (nChan+1));
    fftwf_plan p;
    p = fftwf_plan_dft_r2c_1d(2*nChan, \
                              inP, \
                              reinterpret_cast<fftwf_complex*>(outP), \
                              FFTW_ESTIMATE);
    
    // Data indexing and access
    long secStart;
    
    // Time-domain blanking control
    double cleanFactor;
    long nActFFT;
    
    #ifdef _OPENMP
        #pragma omp parallel default(shared) private(in, out, i, j, k, secStart, cleanFactor, nActFFT)
    #endif
    {
        in = (float*) fftwf_malloc(sizeof(float) * 2*nChan);
        out = (Complex32*) fftwf_malloc(sizeof(Complex32) * (nChan+1));
        
        #ifdef _OPENMP
            #pragma omp for schedule(OMP_SCHEDULER)
        #endif
        for(i=0; i<nStand; i++) {
            nActFFT = 0;
            
            for(j=0; j<nFFT; j++) {
                cleanFactor = 1.0;
                secStart = nSamps * i + 2*nChan*j/Overlap;
                
                for(k=0; k<2*nChan; k++) {
                    in[k] = (float) *(data + secStart + k);
                    
                    if( Clip && fabs(in[k]) >= Clip ) {
                        cleanFactor = 0.0;
                    }
                    
                    if( window != NULL ) {
                        in[k] *= *(window + k);
                    }
                }
                
                fftwf_execute_dft_r2c(p, \
                                      in, \
                                      reinterpret_cast<fftwf_complex*>(out));
                
                for(k=0; k<nChan; k++) {
                    *(psd + nChan*i + k) += cleanFactor*abs2(out[k]);
                }
                
                nActFFT += (long) cleanFactor;
            }
            
            // Scale FFTs
            blas_scal(nChan, 1.0/(2*nChan*nActFFT), (psd + i*nChan), 1);
        }
        
        fftwf_free(in);
        fftwf_free(out);
    }
    fftwf_destroy_plan(p);
    fftwf_free(inP);
    fftwf_free(outP);
    
    Py_END_ALLOW_THREADS
}


template<typename InType, typename OutType>
void compute_pfb_real(long nStand,
                      long nSamps,
                      long nFFT,
                      int nChan,
                      int nTap,
                      int Overlap, 
                      int Clip, 
                      InType const* data,
                      double const* window,
                      OutType* psd) {
    // Setup
    long i, j, k, l;
    
    Py_BEGIN_ALLOW_THREADS
    
    // Create the FFTW plan                          
    float *inP, *in;                          
    Complex32 *outP, *out;
    inP = (float*) fftwf_malloc(sizeof(float) * 2*nChan*nTap);
    outP = (Complex32*) fftwf_malloc(sizeof(Complex32) * (nChan+1)*nTap);
    fftwf_plan p;
    int n[] = {2*nChan,};
    p = fftwf_plan_many_dft_r2c(1, n, nTap, \
                                inP, NULL, 1, 2*nChan, \
                                reinterpret_cast<fftwf_complex*>(outP), NULL, 1, nChan+1, \
                                FFTW_ESTIMATE);
    
    // Data indexing and access
    long secStart;
    
    // Time-domain blanking control
    double cleanFactor;
    long nActFFT;
    
    #ifdef _OPENMP
        #pragma omp parallel default(shared) private(in, out, i, j, k, l, secStart, cleanFactor, nActFFT)
    #endif
    {
        in = (float*) fftwf_malloc(sizeof(float) * 2*nChan*nTap);
        out = (Complex32*) fftwf_malloc(sizeof(Complex32) * (nChan+1)*nTap);
        
        #ifdef _OPENMP
            #pragma omp for schedule(OMP_SCHEDULER)
        #endif
        for(i=0; i<nStand; i++) {
            nActFFT = 0;
            
            for(j=0; j<nFFT; j++) {
                cleanFactor = 1.0;
                secStart = nSamps * i + 2*nChan*j/Overlap;
                
                for(k=0; k<2*nChan*nTap; k++) {
                    if( secStart - 2*nChan*(nTap-1) + k < nSamps*i ) {
                        in[k] = 0.0;
                    } else {
                        in[k] = (float) *(data + secStart - 2*nChan*(nTap-1) + k);
                    }
                    
                    if( Clip && fabs(in[k]) >= Clip ) {
                        cleanFactor = 0.0;
                    }
                    
                    if( window != NULL ) {
                      in[k] *= *(window + k);
                    }
                }
                
                fftwf_execute_dft_r2c(p, \
                                      in, \
                                      reinterpret_cast<fftwf_complex*>(out));
                
                for(l=1; l<nTap; l++) {
                    for(k=0; k<nChan; k++) {
                        out[k] += out[k+l*(nChan+1)];
                    }
                }
                
                for(k=0; k<nChan; k++) {
                    *(psd + nChan*i + k) += cleanFactor*abs2(out[k]);
                }
                
                nActFFT += (long) cleanFactor;
            }
            
            // Scale FFTs
            blas_scal(nChan, 1.0/(2*nChan*nActFFT), (psd + i*nChan), 1);
        }
        
        fftwf_free(in);
        fftwf_free(out);
    }
    fftwf_destroy_plan(p);
    fftwf_free(inP);
    fftwf_free(outP);
    
    Py_END_ALLOW_THREADS
}


template<typename InType, typename OutType>
void compute_psd_complex(long nStand,
                         long nSamps,
                         long nFFT,
                         int nChan, 
                         int Overlap, 
                         int Clip, 
                         InType const* data,
                         double const* window,
                         OutType* psd) {
    // Setup
    long i, j, k;
    
    Py_BEGIN_ALLOW_THREADS
    
    // Create the FFTW plan
    Complex32 *inP, *in;
    inP = (Complex32*) fftwf_malloc(sizeof(Complex32) * nChan);
    fftwf_plan p;
    p = fftwf_plan_dft_1d(nChan, \
                          reinterpret_cast<fftwf_complex*>(inP), \
                          reinterpret_cast<fftwf_complex*>(inP), \
                          FFTW_FORWARD, FFTW_ESTIMATE);
    
    // Data indexing and access
    long secStart;
    double* temp2;
    // Time-domain blanking control
    double cleanFactor;
    long nActFFT;
    
    #ifdef _OPENMP
        #pragma omp parallel default(shared) private(in, i, j, k, secStart, cleanFactor, nActFFT, temp2)
    #endif
    {
        in = (Complex32*) fftwf_malloc(sizeof(Complex32) * nChan);
        temp2 = (double*) malloc(sizeof(double) * (nChan/2+nChan%2));
        
        #ifdef _OPENMP
            #pragma omp for schedule(OMP_SCHEDULER)
        #endif
        for(i=0; i<nStand; i++) {
            nActFFT = 0;
            
            for(j=0; j<nFFT; j++) {
                cleanFactor = 1.0;
                secStart = nSamps * i + nChan*j/Overlap;
                
                for(k=0; k<nChan; k++) {
                    in[k] = Complex32(*(data + 2*secStart + 2*k + 0), \
                                      *(data + 2*secStart + 2*k + 1));
                    
                    if( Clip && abs(in[k]) >= Clip ) {
                        cleanFactor = 0.0;
                    }
                    
                    if( window != NULL ) {
                        in[k] *= *(window + k);
                    }
                }
                
                fftwf_execute_dft(p, \
                                  reinterpret_cast<fftwf_complex*>(in), \
                                  reinterpret_cast<fftwf_complex*>(in));
                
                for(k=0; k<nChan; k++) {
                    *(psd + nChan*i + k) += cleanFactor*abs2(in[k]);
                }
                
                nActFFT += (long) cleanFactor;
            }
            
            // Shift FFTs
            memcpy(temp2, psd + nChan*i, sizeof(OutType)*(nChan/2+nChan%2));
            memmove(psd + nChan*i, psd + nChan*i+nChan/2+nChan%2, sizeof(OutType)*nChan/2);
            memcpy(psd + nChan*i+nChan/2, temp2, sizeof(OutType)*(nChan/2+nChan%2));
            
            // Scale FFTs
            blas_scal(nChan, 1.0/(nActFFT*nChan), (psd + i*nChan), 1);
        }
        
        fftwf_free(in);
        free(temp2);
    }
    fftwf_destroy_plan(p);
    fftwf_free(inP);
    
    Py_END_ALLOW_THREADS
}


template<typename InType, typename OutType>
void compute_pfb_complex(long nStand,
                         long nSamps,
                         long nFFT,
                         int nChan,
                         int nTap,
                         int Overlap, 
                         int Clip, 
                         InType const* data,
                         double const* window,
                         OutType* psd) {
    // Setup
    long i, j, k, l;
    
    Py_BEGIN_ALLOW_THREADS
    
    // Create the FFTW plan
    Complex32 *inP, *in;
    inP = (Complex32*) fftwf_malloc(sizeof(Complex32) * nChan*nTap);
    fftwf_plan p;
    int n[] = {nChan,};
    p = fftwf_plan_many_dft(1, n, nTap, \
                            reinterpret_cast<fftwf_complex*>(inP), NULL, 1, nChan, \
                            reinterpret_cast<fftwf_complex*>(inP), NULL, 1, nChan, \
                            FFTW_FORWARD, FFTW_ESTIMATE);
    
    // Data indexing and access
    long secStart;
    double* temp2;
    // Time-domain blanking control
    double cleanFactor;
    long nActFFT;
    
    #ifdef _OPENMP
        #pragma omp parallel default(shared) private(in, i, j, k, l, secStart, cleanFactor, nActFFT, temp2)
    #endif
    {
        in = (Complex32*) fftwf_malloc(sizeof(Complex32) * nChan*nTap);
        temp2 = (double*) malloc(sizeof(double) * (nChan/2+nChan%2));
        
        #ifdef _OPENMP
            #pragma omp for schedule(OMP_SCHEDULER)
        #endif
        for(i=0; i<nStand; i++) {
            nActFFT = 0;
            
            for(j=0; j<nFFT; j++) {
                cleanFactor = 1.0;
                secStart = nSamps * i + nChan*j/Overlap;
                
                for(k=0; k<nChan*nTap; k++) {
                    if( secStart - nChan*(nTap-1) + k < nSamps*i ) {
                        in[k] = 0.0;
                    } else {
                        in[k] = Complex32(*(data + 2*secStart - 2*nChan*(nTap-1) + 2*k + 0), \
                                          *(data + 2*secStart - 2*nChan*(nTap-1) + 2*k + 1));
                    }
                    
                    if( Clip && abs(in[k]) >= Clip ) {
                        cleanFactor = 0.0;
                    }
                    
                    if( window != NULL ) {
                      in[k] *= *(window + k);
                    }
                }
                
                fftwf_execute_dft(p, \
                                  reinterpret_cast<fftwf_complex*>(in), \
                                  reinterpret_cast<fftwf_complex*>(in));
                
                for(l=1; l<nTap; l++) {
                    for(k=0; k<nChan; k++) {
                        in[k] += in[k+l*nChan];
                    }
                }
                
                for(k=0; k<nChan; k++) {
                    *(psd + nChan*i + k) += cleanFactor*abs2(in[k]);
                }
                
                nActFFT += (long) cleanFactor;
            }
            
            // Shift FFTs
            memcpy(temp2, psd + nChan*i, sizeof(OutType)*(nChan/2+nChan%2));
            memmove(psd + nChan*i, psd + nChan*i+nChan/2+nChan%2, sizeof(OutType)*nChan/2);
            memcpy(psd + nChan*i+nChan/2, temp2, sizeof(OutType)*(nChan/2+nChan%2));
            
            // Scale FFTs
            blas_scal(nChan, 1.0/(nActFFT*nChan), (psd + i*nChan), 1);
        }
        
        fftwf_free(in);
        free(temp2);
    }
    fftwf_destroy_plan(p);
    fftwf_free(inP);
    
    Py_END_ALLOW_THREADS
}


static PyObject *FPSD(PyObject *self, PyObject *args, PyObject *kwds) {
    PyObject *signals, *signalsF, *window=Py_None, *arglist, *windowValue;
    PyArrayObject *data=NULL, *dataF=NULL, *windowData=NULL;
    int isReal;
    int nChan = 64;
    int Overlap = 1;
    int Clip = 0;
    
    long nStand, nSamps, nFFT;
    
    char const* kwlist[] = {"signals", "LFFT", "overlap", "clip_level", "window", NULL};
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|iiiO:set_callback", const_cast<char **>(kwlist), &signals, &nChan, &Overlap, &Clip, &window)) {
        PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
        goto fail;
    } else {
        if(!PyCallable_Check(window) && window != Py_None) {
            PyErr_Format(PyExc_TypeError, "window must be a callable function or None");
            goto fail;
        }
        Py_XINCREF(window);
        Py_XDECREF(windowFunc);
        windowFunc = window;
    }
    
    // Bring the data into C and make it usable
    data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, 
                                                        PyArray_TYPE((PyArrayObject *) signals), 
                                                        2, 2);
    if( data == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input array signals as a 2-D array");
        goto fail;
    }
    
    // Get the properties of the data
    nStand = (long) PyArray_DIM(data, 0);
    nSamps = (long) PyArray_DIM(data, 1);
    isReal = 1 - PyArray_ISCOMPLEX(data);
    
    // Calculate the windowing function
    if( windowFunc != Py_None ) {
        arglist = Py_BuildValue("(i)", (1+isReal)*nChan);
        windowValue = PyObject_CallObject(windowFunc, arglist);
        windowData = (PyArrayObject *) PyArray_ContiguousFromObject(windowValue, NPY_DOUBLE, 1, 1);
        Py_DECREF(arglist);
        Py_DECREF(windowValue);
    }
    
    // Find out how large the output array needs to be and initialize it
    nFFT = nSamps / (((1+isReal)*nChan)/Overlap) - ((1+isReal)*nChan)/(((1+isReal)*nChan)/Overlap) + 1;
    npy_intp dims[2];
    dims[0] = (npy_intp) nStand;
    dims[1] = (npy_intp) nChan;
    dataF = (PyArrayObject*) PyArray_ZEROS(2, dims, NPY_DOUBLE, 0);
    if(dataF == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array");
        goto fail;
    }
    
#define LAUNCH_PSD_REAL(IterType) \
        compute_pfb_real<IterType>(nStand, nSamps, nFFT, nChan, 1, Overlap, Clip, \
                                   (IterType*) PyArray_DATA(data), \
                                   (double*) PyArray_SAFE_DATA(windowData), \
                                   (double*) PyArray_DATA(dataF))
#define LAUNCH_PSD_COMPLEX(IterType) \
        compute_pfb_complex<IterType>(nStand, nSamps, nFFT, nChan, 1, Overlap, Clip, \
                                      (IterType*) PyArray_DATA(data), \
                                      (double*) PyArray_SAFE_DATA(windowData), \
                                      (double*) PyArray_DATA(dataF))
    
    switch( PyArray_TYPE(data) ){
        case( NPY_INT8       ): LAUNCH_PSD_REAL(int8_t);    break;
        case( NPY_INT16      ): LAUNCH_PSD_REAL(int16_t);   break;
        case( NPY_INT32      ): LAUNCH_PSD_REAL(int);       break;
        case( NPY_INT64      ): LAUNCH_PSD_REAL(long);      break;
        case( NPY_FLOAT32    ): LAUNCH_PSD_REAL(float);     break;
        case( NPY_FLOAT64    ): LAUNCH_PSD_REAL(double);    break;
        case( NPY_COMPLEX64  ): LAUNCH_PSD_COMPLEX(float);  break;
        case( NPY_COMPLEX128 ): LAUNCH_PSD_COMPLEX(double); break;
        default: PyErr_Format(PyExc_RuntimeError, "Unsupport input data type"); goto fail;
    }
        
#undef LAUNCH_PSD_REAL
#undef LAUNCH_PSD_COMPLEX
    
    signalsF = Py_BuildValue("O", PyArray_Return(dataF));
    
    Py_XDECREF(data);
    Py_XDECREF(windowData);
    Py_XDECREF(dataF);

    return signalsF;
    
fail:
    Py_XDECREF(data);
    Py_XDECREF(windowData);
    Py_XDECREF(dataF);
    
    return NULL;
}

PyDoc_STRVAR(FPSD_doc, \
"Perform a series of Fourier transforms with windows on data to get the\n\
PSD.\n\
\n\
Input arguments are:\n\
 * signals: 2-D numpy (stands by samples) array of data to FFT\n\
\n\
Input keywords are:\n\
 * LFFT: number of FFT channels to make (default=64)\n\
 * overlap: number of overlapped FFTs to use (default=1)\n\
 * window: Callable Python function for generating the window or None for\n\
           no window\n\
 * clip_level: count value of 'bad' data.  FFT windows with instantaneous\n\
               powers greater than or equal to this value greater are zeroed.  \n\
               Setting the ClipLevel to zero disables time-domain blanking\n\
\n\
Outputs:\n\
 * psd: 2-D numpy.double (stands by channels) of PSD data\n\
");


static PyObject *PFBPSD(PyObject *self, PyObject *args, PyObject *kwds) {
    PyObject *signals, *signalsF, *window=Py_None, *arglist, *windowValue;
    PyArrayObject *data=NULL, *dataF=NULL, *windowData=NULL;
    int isReal;
    int nChan = 64;
    int nTap = 4;
    int Overlap = 1;
    int Clip = 0;
    
    long nStand, nSamps, nFFT;
    
    char const* kwlist[] = {"signals", "LFFT", "overlap", "clip_level", "window", NULL};
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|iiiO:set_callback", const_cast<char **>(kwlist), &signals, &nChan, &Overlap, &Clip, &window)) {
        PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
        goto fail;
    } else {
        if(!PyCallable_Check(window) && window != Py_None) {
            PyErr_Format(PyExc_TypeError, "window must be a callable function or None");
            goto fail;
        }
        Py_XINCREF(window);
        Py_XDECREF(windowFunc);
        windowFunc = window;
    }
    
    // Bring the data into C and make it usable
    data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, 
                                                        PyArray_TYPE((PyArrayObject *) signals), 
                                                        2, 2);
    if( data == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input array signals as a 2-D array");
        goto fail;
    }
    
    // Get the properties of the data
    nStand = (long) PyArray_DIM(data, 0);
    nSamps = (long) PyArray_DIM(data, 1);
    isReal = 1 - PyArray_ISCOMPLEX(data);
    
    // Calculate the windowing function for the PFB
    double *pfb;
    pfb = (double*) malloc(sizeof(float) * nChan*nTap);
    for(int i=0; i<(1+isReal)*nChan*nTap; i++) {
        *(pfb + i) = sinc((i - (1+isReal)*nChan*nTap/2.0 + 0.5)/((1+isReal)*nChan));
        *(pfb + i) *= hamming(2*NPY_PI*i/((1+isReal)*nChan*nTap));
    }
    
    // Find out how large the output array needs to be and initialize it
    nFFT = nSamps / (((1+isReal)*nChan)/Overlap) - ((1+isReal)*nChan)/(((1+isReal)*nChan)/Overlap) + 1;
    npy_intp dims[2];
    dims[0] = (npy_intp) nStand;
    dims[1] = (npy_intp) nChan;
    dataF = (PyArrayObject*) PyArray_ZEROS(2, dims, NPY_DOUBLE, 0);
    if(dataF == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array");
        free(pfb);
        goto fail;
    }
    
#define LAUNCH_PFB_REAL(IterType) \
        compute_pfb_real<IterType>(nStand, nSamps, nFFT, nChan, nTap, Overlap, Clip, \
                                   (IterType*) PyArray_DATA(data), \
                                   pfb, \
                                   (double*) PyArray_DATA(dataF))
#define LAUNCH_PFB_COMPLEX(IterType) \
        compute_pfb_complex<IterType>(nStand, nSamps, nFFT, nChan, nTap, Overlap, Clip, \
                                      (IterType*) PyArray_DATA(data), \
                                      pfb, \
                                      (double*) PyArray_DATA(dataF))
    
    switch( PyArray_TYPE(data) ){
        case( NPY_INT8       ): LAUNCH_PFB_REAL(int8_t);    break;
        case( NPY_INT16      ): LAUNCH_PFB_REAL(int16_t);   break;
        case( NPY_INT32      ): LAUNCH_PFB_REAL(int);       break;
        case( NPY_INT64      ): LAUNCH_PFB_REAL(long);      break;
        case( NPY_FLOAT32    ): LAUNCH_PFB_REAL(float);     break;
        case( NPY_FLOAT64    ): LAUNCH_PFB_REAL(double);    break;
        case( NPY_COMPLEX64  ): LAUNCH_PFB_COMPLEX(float);  break;
        case( NPY_COMPLEX128 ): LAUNCH_PFB_COMPLEX(double); break;
        default: PyErr_Format(PyExc_RuntimeError, "Unsupport input data type"); goto fail;
    }
        
#undef LAUNCH_PFB_REAL
#undef LAUNCH_PFB_COMPLEX
    
    free(pfb);
    
    signalsF = Py_BuildValue("O", PyArray_Return(dataF));
    
    Py_XDECREF(data);
    Py_XDECREF(windowData);
    Py_XDECREF(dataF);

    return signalsF;
    
fail:
    Py_XDECREF(data);
    Py_XDECREF(windowData);
    Py_XDECREF(dataF);
    
    return NULL;
}

PyDoc_STRVAR(PFBPSD_doc, \
"Perform a series of polyphase filter bank transforms (4-tap plus a\n\
Hanning window) on data to get the PSD.\n\
\n\
Input arguments are:\n\
 * signals: 2-D numpy (stands by samples) array of data to FFT\n\
\n\
Input keywords are:\n\
 * LFFT: number of FFT channels to make (default=64)\n\
 * overlap: number of overlapped FFTs to use (default=1)\n\
 * window: Callable Python function for generating the window or None for\n\
           no window\n\
 * clip_level: count value of 'bad' data.  FFT windows with instantaneous\n\
               powers greater than or equal to this value greater are zeroed.  \n\
               Setting the ClipLevel to zero disables time-domain blanking\n\
\n\
Outputs:\n\
 * psd: 2-D numpy.double (stands by channels) of PSD data\n\
");


/*
  Module Setup - Function Definitions and Documentation
*/

static PyMethodDef SpecMethods[] = {
    {"FPSD",   (PyCFunction) FPSD,   METH_VARARGS|METH_KEYWORDS, FPSD_doc },
    {"PFBPSD", (PyCFunction) PFBPSD, METH_VARARGS|METH_KEYWORDS, PFBPSD_doc},
    {NULL,      NULL,                0,                          NULL     }
};

PyDoc_STRVAR(spec_doc, \
"Extension to take timeseries data and convert it to the frequency domain.\n\
\n\
The functions defined in this module are:\n\
 * FPSD -  FFT and integrate function for computing a series of overlapped\n\
           Fourier transforms for a real-valued (TBW) or complex-valued (TBN\n\
           and DRX) signals from a collection of stands all at once.\n\
 * PFBPSD - Similar to FPSD, but using a 4-tap + Hanning windowed polyphase\n\
            filter bank.\n\
\n\
See the inidividual functions for more details.\n\
\n\
.. versionchanged:: 1.3.0\n\
\tMerged FPSDR and FPSDC into the new FPSD\n\
\n\
.. versionchanged:: 1.0.1\n\
\tRemoved the polyphase filterbank versions of the four core functions.\n\
");


/*
  Module Setup - Initialization
*/

MOD_INIT(_spec) {
    char filename[256];
    PyObject *m, *all, *pModule, *pDataPath=NULL;
    
    Py_Initialize();
    
    // Module definitions and functions
    MOD_DEF(m, "_spec", SpecMethods, spec_doc);
    if( m == NULL ) {
        return MOD_ERROR_VAL;
    }
    import_array();
    
    // Version and revision information
    PyModule_AddObject(m, "__version__", PyString_FromString("0.7"));
    
    // Function listings
    all = PyList_New(0);
    PyList_Append(all, PyString_FromString("FPSD"));
    PyList_Append(all, PyString_FromString("PFBPSD"));
    PyModule_AddObject(m, "__all__", all);
    
    // LSL FFTW Wisdom
    pModule = PyImport_ImportModule("lsl.common.paths");
    if( pModule != NULL ) {
        pDataPath = PyObject_GetAttrString(pModule, "WISDOM");
        if( pDataPath != NULL ) {
            sprintf(filename, "%s/fftwf_wisdom.txt", PyString_AsString(pDataPath));
            read_wisdom(filename, m);
        }
    } else {
        PyErr_Warn(PyExc_RuntimeWarning, "Cannot load the LSL FFTWF wisdom");
    }
    Py_XDECREF(pDataPath);
    Py_XDECREF(pModule);
    
    return MOD_SUCCESS_VAL(m);
}

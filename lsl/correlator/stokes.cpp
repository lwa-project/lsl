#include "Python.h"
#include <cmath>
#include <cstdint>
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

#include "common.hpp"
#include "blas.hpp"


/*
  Holder for window function callback
*/

static PyObject *windowFunc = NULL;


/*
  FFT Functions
    1. FPSD - FFT a real or complex-valued collection of signals
*/

template<typename InType, typename OutType>
void compute_stokes_real(long nStand,
                         long nSamps,
                         long nFFT,
                         int nChan,
                         int nTap,
                         int Overlap,
                         int Clip,
                         InType const* dataX,
                         InType const* dataY,
                         double const* window,
                         OutType* psd) {
    // Setup
    long i, j, k, l;
    
    Py_BEGIN_ALLOW_THREADS
    
    // Create the FFTW plan                          
    float *inP, *inX, *inY;                          
    Complex32 *outP, *outX, *outY;
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
        #pragma omp parallel default(shared) private(inX, inY, outX, outY, i, j, k, l, secStart, cleanFactor, nActFFT)
    #endif
    {
        inX = (float*) fftwf_malloc(sizeof(float) * 2*nChan*nTap);
        inY = (float*) fftwf_malloc(sizeof(float) * 2*nChan*nTap);
        outX = (Complex32*) fftwf_malloc(sizeof(Complex32) * (nChan+1)*nTap);
        outY = (Complex32*) fftwf_malloc(sizeof(Complex32) * (nChan+1)*nTap);
        
        #ifdef _OPENMP
            #pragma omp for schedule(OMP_SCHEDULER)
        #endif
        for(i=0; i<nStand; i++) {
            nActFFT = 0;
            
            for(j=0; j<nFFT; j++) {
                cleanFactor = 1.0;
                secStart = nSamps * i + 2*nChan*j/Overlap;
                
                for(k=0; k<2*nChan*nTap; k+=2) {
                    if( secStart - 2*nChan*(nTap-1) + k < nSamps*i ) {
                        inX[k] = 0.0;
                        inY[k] = 0.0;
                    } else {
                        inX[k] = (float) *(dataX + secStart - 2*nChan*(nTap-1) + k);
                        inY[k] = (float) *(dataY + secStart - 2*nChan*(nTap-1) + k);
                    }
                    if( secStart - 2*nChan*(nTap-1) + k + 1 < nSamps*i ) {
                        inX[k+1] = 0.0;
                        inY[k+1] = 0.0;
                    } else {
                        inX[k+1] = (float) *(dataX + secStart - 2*nChan*(nTap-1) + k + 1);
                        inY[k+1] = (float) *(dataY + secStart - 2*nChan*(nTap-1) + k + 1);
                    }
                    
                    if( Clip && (   fabs(inX[k]) >= Clip || fabs(inY[k]) >= Clip \
                                 || fabs(inX[k+1]) >= Clip || fabs(inY[k+1]) >= Clip) ) {
                        cleanFactor = 0.0;
                    }
                    
                    if( window != NULL ) {
                        inX[k] *= *(window + k);
                        inY[k] *= *(window + k);
                        inX[k+1] *= *(window + k + 1);
                        inY[k+1] *= *(window + k + 1);
                    }
                }
                
                fftwf_execute_dft_r2c(p, \
                                      inX, \
                                      reinterpret_cast<fftwf_complex*>(outX));
                fftwf_execute_dft_r2c(p, \
                                      inY, \
                                      reinterpret_cast<fftwf_complex*>(outY));
                
                for(l=1; l<nTap; l++) {
                    for(k=0; k<nChan; k++) {
                        outX[k] += outX[k+l*(nChan+1)];
                        outY[k] += outY[k+l*(nChan+1)];
                    }
                }
                
                for(k=0; k<nChan; k++) {
                    // I
                    *(psd + 0*nChan*nStand + nChan*i + k) += cleanFactor*abs2(outX[k]);
                    *(psd + 0*nChan*nStand + nChan*i + k) += cleanFactor*abs2(outY[k]);
                    
                    // Q
                    *(psd + 1*nChan*nStand + nChan*i + k) += cleanFactor*abs2(outX[k]);
                    *(psd + 1*nChan*nStand + nChan*i + k) -= cleanFactor*abs2(outY[k]);
                    
                    // U
                    *(psd + 2*nChan*nStand + nChan*i + k) += 2*cleanFactor*outX[k].real()*outY[k].real();
                    *(psd + 2*nChan*nStand + nChan*i + k) += 2*cleanFactor*outX[k].imag()*outY[k].imag();
                    
                    // V
                    *(psd + 3*nChan*nStand + nChan*i + k) +=2*cleanFactor*outX[k].imag()*outY[k].real();
                    *(psd + 3*nChan*nStand + nChan*i + k) -=2*cleanFactor*outX[k].real()*outY[k].imag();
                }
                
                nActFFT += (long) cleanFactor;
            }
            
            // Scale FFTs
            for(j=0; j<4; j++) {
                blas_scal(nChan, 1.0/(2*nChan*nActFFT), (psd + j*nChan*nStand + nChan*i), 1);
            }
        }
        
        fftwf_free(inX);
        fftwf_free(inY);
        fftwf_free(outX);
        fftwf_free(outY);
    }
    fftwf_destroy_plan(p);
    fftwf_free(inP);
    fftwf_free(outP);
    
    Py_END_ALLOW_THREADS
}


template<typename InType, typename OutType>
void compute_stokes_complex(long nStand,
                            long nSamps,
                            long nFFT,
                            int nChan,
                            int nTap,
                            int Overlap,
                            int Clip,
                            InType const* dataX,
                            InType const* dataY,
                            double const* window,
                            OutType* psd) {
    // Setup
    long i, j, k, l;
    
    Py_BEGIN_ALLOW_THREADS
    
    // Create the FFTW plan
    Complex32 *inP, *inX, *inY;
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
        #pragma omp parallel default(shared) private(inX, inY, i, j, k, l, secStart, cleanFactor, nActFFT, temp2)
    #endif
    {
        inX = (Complex32*) fftwf_malloc(sizeof(Complex32) * nChan*nTap);
        inY = (Complex32*) fftwf_malloc(sizeof(Complex32) * nChan*nTap);
        temp2 = (double*) aligned64_malloc(sizeof(double) * (nChan/2+nChan%2));
        
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
                        inX[k] = 0.0;
                        inY[k] = 0.0;
                    } else {
                        inX[k] = Complex32(*(dataX + 2*secStart - 2*nChan*(nTap-1) + 2*k + 0), \
                                           *(dataX + 2*secStart - 2*nChan*(nTap-1) + 2*k + 1));
                        inY[k] = Complex32(*(dataY + 2*secStart - 2*nChan*(nTap-1) + 2*k + 0), \
                                           *(dataY + 2*secStart - 2*nChan*(nTap-1) + 2*k + 1));
                    }
                    
                    if( Clip && ( abs(inX[k]) >= Clip || abs(inY[k]) >= Clip ) ) {
                        cleanFactor = 0.0;
                    }
                    
                    if( window != NULL ) {
                        inX[k] *= *(window + k);
                        inY[k] *= *(window + k);
                    }
                }
                
                fftwf_execute_dft(p, \
                                  reinterpret_cast<fftwf_complex*>(inX), \
                                  reinterpret_cast<fftwf_complex*>(inX));
                fftwf_execute_dft(p,  \
                                  reinterpret_cast<fftwf_complex*>(inY), \
                                  reinterpret_cast<fftwf_complex*>(inY));
                
                for(l=1; l<nTap; l++) {
                    for(k=0; k<nChan; k++) {
                        inX[k] += inX[k+l*nChan];
                        inY[k] += inY[k+l*nChan];
                    }
                }
                
                for(k=0; k<nChan; k++) {
                    // I
                    *(psd + 0*nChan*nStand + nChan*i + k) += cleanFactor*abs2(inX[k]);
                    *(psd + 0*nChan*nStand + nChan*i + k) += cleanFactor*abs2(inY[k]);
                    
                    // Q
                    *(psd + 1*nChan*nStand + nChan*i + k) += cleanFactor*abs2(inX[k]);
                    *(psd + 1*nChan*nStand + nChan*i + k) -= cleanFactor*abs2(inY[k]);
                    
                    // U
                    *(psd + 2*nChan*nStand + nChan*i + k) += 2*cleanFactor*inX[k].real()*inY[k].real();
                    *(psd + 2*nChan*nStand + nChan*i + k) += 2*cleanFactor*inX[k].imag()*inY[k].imag();
                    
                    // V
                    *(psd + 3*nChan*nStand + nChan*i + k) +=2*cleanFactor*inX[k].imag()*inY[k].real();
                    *(psd + 3*nChan*nStand + nChan*i + k) -=2*cleanFactor*inX[k].real()*inY[k].imag();
                }
                
                nActFFT += (long) cleanFactor;
            }
            
            for(j=0; j<4; j++) {
                // Shift FFTs
                memcpy(temp2, (psd + j*nChan*nStand + nChan*i), sizeof(OutType)*(nChan/2+nChan%2));
                memmove((psd + j*nChan*nStand + nChan*i), (psd + j*nChan*nStand + nChan*i)+nChan/2+nChan%2, sizeof(OutType)*nChan/2);
                memcpy((psd + j*nChan*nStand + nChan*i)+nChan/2, temp2, sizeof(OutType)*(nChan/2+nChan%2));
                
                // Scale FFTs
                blas_scal(nChan, 1.0/(nActFFT*nChan), (psd + j*nChan*nStand + nChan*i), 1);
            }
        }
        
        fftwf_free(inX);
        fftwf_free(inY);
        aligned64_free(temp2);
    }
    fftwf_destroy_plan(p);
    fftwf_free(inP);
    
    Py_END_ALLOW_THREADS
}


static PyObject *FPSD(PyObject *self, PyObject *args, PyObject *kwds) {
    PyObject *signalsX, *signalsY, *signalsF, *window=Py_None, *arglist, *windowValue;
    PyArrayObject *dataX=NULL, *dataY=NULL, *dataF=NULL, *windowData=NULL;
    int isReal;
    int nChan = 64;
    int Overlap = 1;
    int Clip = 0;
    
    long nStand, nSamps, nFFT;
    
    char const* kwlist[] = {"signalsX", "signalsY", "LFFT", "overlap", "clip_level", "window", NULL};
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "OO|iiiO:set_callback", const_cast<char **>(kwlist), &signalsX, &signalsY, &nChan, &Overlap, &Clip, &window)) {
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
    dataX = (PyArrayObject *) PyArray_ContiguousFromObject(signalsX, 
                                                        PyArray_TYPE((PyArrayObject *) signalsX), 
                                                        2, 3);
    dataY = (PyArrayObject *) PyArray_ContiguousFromObject(signalsY, 
                                                        PyArray_TYPE((PyArrayObject *) signalsX), 
                                                        2, 3);
    if( dataX == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input array signalsX as a 2-D array");
        goto fail;
    }
    if( dataY == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input array signalsY as a 2-D array");
        goto fail;
    }
    
    // Get the properties of the data
    nStand = (long) PyArray_DIM(dataX, 0);
    nSamps = (long) PyArray_DIM(dataX, 1);
    isReal = 1 - PyArray_LSL_ISCOMPLEX(dataX, 2);
    if( PyArray_NDIM(dataX) == 3 && PyArray_TYPE(dataX) != NPY_INT8 ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input signalsX array as a 2-D array");
        goto fail;
    }
    if( PyArray_NDIM(dataY) == 3 && PyArray_TYPE(dataY) != NPY_INT8 ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input signalsY array as a 2-D array");
        goto fail;
    }
    
    // Make sure the dimensions of X and Y agree
    if( PyArray_DIM(dataY, 0) != nStand ) {
        PyErr_Format(PyExc_RuntimeError, "X and Y signals have different stand counts");
        goto fail;
    }
    if( PyArray_DIM(dataY, 1) != nSamps ) {
        PyErr_Format(PyExc_RuntimeError, "X and Y signals have different sample counts");
        goto fail;
    }
    
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
    npy_intp dims[3];
    dims[0] = (npy_intp) 4;
    dims[1] = (npy_intp) nStand;
    dims[2] = (npy_intp) nChan;
    dataF = (PyArrayObject*) PyArray_ZEROS(3, dims, NPY_DOUBLE, 0);
    if(dataF == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array");
        goto fail;
    }
    
#define LAUNCH_PSD_REAL(IterType) \
        compute_stokes_real<IterType>(nStand, nSamps, nFFT, nChan, 1, Overlap, Clip, \
                                      (IterType*) PyArray_DATA(dataX), \
                                      (IterType*) PyArray_DATA(dataY), \
                                      (double*) PyArray_SAFE_DATA(windowData), \
                                      (double*) PyArray_DATA(dataF))
#define LAUNCH_PSD_COMPLEX(IterType) \
        compute_stokes_complex<IterType>(nStand, nSamps, nFFT, nChan, 1, Overlap, Clip, \
                                         (IterType*) PyArray_DATA(dataX), \
                                         (IterType*) PyArray_DATA(dataY), \
                                         (double*) PyArray_SAFE_DATA(windowData), \
                                         (double*) PyArray_DATA(dataF))
    
    switch( PyArray_LSL_TYPE(dataX, 2) ){
        case( NPY_INT8       ): LAUNCH_PSD_REAL(int8_t);    break;
        case( NPY_INT16      ): LAUNCH_PSD_REAL(int16_t);   break;
        case( NPY_INT32      ): LAUNCH_PSD_REAL(int);       break;
        case( NPY_INT64      ): LAUNCH_PSD_REAL(long);      break;
        case( NPY_FLOAT32    ): LAUNCH_PSD_REAL(float);     break;
        case( NPY_FLOAT64    ): LAUNCH_PSD_REAL(double);    break;
        case( LSL_CI8        ): LAUNCH_PSD_COMPLEX(int8_t); break;
        case( NPY_COMPLEX64  ): LAUNCH_PSD_COMPLEX(float);  break;
        case( NPY_COMPLEX128 ): LAUNCH_PSD_COMPLEX(double); break;
        default: PyErr_Format(PyExc_RuntimeError, "Unsupport input data type"); goto fail;
    }
        
#undef LAUNCH_PSD_REAL
#undef LAUNCH_PSD_COMPLEX
    
    
    signalsF = Py_BuildValue("O", PyArray_Return(dataF));
    
    Py_XDECREF(dataX);
    Py_XDECREF(dataY);
    Py_XDECREF(windowData);
    Py_XDECREF(dataF);
    
    return signalsF;
    
fail:
    Py_XDECREF(dataX);
    Py_XDECREF(dataY);
    Py_XDECREF(windowData);
    Py_XDECREF(dataF);
    
    return NULL;
}

PyDoc_STRVAR(FPSD_doc, \
"Perform a series of Fourier transforms with windows on data to get the\n\
PSD for the four Stokes parameters: I, Q, U, and V.\n\
\n\
Input arguments are:\n\
 * signals: 2-D numpy.int16 (stands by samples) array of data to FFT\n\
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
 * psd: 3-D numpy.double (Stokes parameter (I,Q,U,V) by stands by channels)\n\
        of PSD data\n\
");


static PyObject *PFBPSD(PyObject *self, PyObject *args, PyObject *kwds) {
    PyObject *signalsX, *signalsY, *signalsF, *window=Py_None;
    PyArrayObject *dataX=NULL, *dataY=NULL, *dataF=NULL;
    int isReal;
    int nChan = 64;
    int nTap = PFB_NTAP;
    int Overlap = 1;
    int Clip = 0;
    
    long nStand, nSamps, nFFT;
    double *pfb = NULL;
    
    char const* kwlist[] = {"signalsX", "signalsY", "LFFT", "overlap", "clip_level", "window", NULL};
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "OO|iiiO:set_callback", const_cast<char **>(kwlist), &signalsX, &signalsY, &nChan, &Overlap, &Clip, &window)) {
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
    dataX = (PyArrayObject *) PyArray_ContiguousFromObject(signalsX, 
                                                        PyArray_TYPE((PyArrayObject *) signalsX), 
                                                        2, 3);
    dataY = (PyArrayObject *) PyArray_ContiguousFromObject(signalsY, 
                                                        PyArray_TYPE((PyArrayObject *) signalsX), 
                                                        2, 3);
    if( dataX == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input array signalsX as a 2-D array");
        goto fail;
    }
    if( dataY == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input array signalsY as a 2-D array");
        goto fail;
    }
    
    // Get the properties of the data
    nStand = (long) PyArray_DIM(dataX, 0);
    nSamps = (long) PyArray_DIM(dataX, 1);
    isReal = 1 - PyArray_LSL_ISCOMPLEX(dataX, 2);
    if( PyArray_NDIM(dataX) == 3 && PyArray_TYPE(dataX) != NPY_INT8 ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input signalsX array as a 2-D array");
        goto fail;
    }
    if( PyArray_NDIM(dataY) == 3 && PyArray_TYPE(dataY) != NPY_INT8 ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input signalsY array as a 2-D array");
        goto fail;
    }
    
    // Make sure the dimensions of X and Y agree
    if( PyArray_DIM(dataY, 0) != nStand ) {
        PyErr_Format(PyExc_RuntimeError, "X and Y signals have different stand counts");
        goto fail;
    }
    if( PyArray_DIM(dataY, 1) != nSamps ) {
        PyErr_Format(PyExc_RuntimeError, "X and Y signals have different sample counts");
        goto fail;
    }
    
    // Calculate the windowing function for the PFB
    pfb = (double*) aligned64_malloc(sizeof(double) * (1+isReal)*nChan*nTap);
    for(int i=0; i<(1+isReal)*nChan*nTap; i++) {
        *(pfb + i) = sinc((i - (1+isReal)*nChan*nTap/2.0 + 0.5)/((1+isReal)*nChan));
        *(pfb + i) *= hamming(2*NPY_PI*i/((1+isReal)*nChan*nTap));
    }
    
    // Find out how large the output array needs to be and initialize it
    nFFT = nSamps / (((1+isReal)*nChan)/Overlap) - ((1+isReal)*nChan)/(((1+isReal)*nChan)/Overlap) + 1;
    npy_intp dims[3];
    dims[0] = (npy_intp) 4;
    dims[1] = (npy_intp) nStand;
    dims[2] = (npy_intp) nChan;
    dataF = (PyArrayObject*) PyArray_ZEROS(3, dims, NPY_DOUBLE, 0);
    if(dataF == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array");
        goto fail;
    }
    
#define LAUNCH_PFB_REAL(IterType) \
        compute_stokes_real<IterType>(nStand, nSamps, nFFT, nChan, nTap, Overlap, Clip, \
                                      (IterType*) PyArray_DATA(dataX), \
                                      (IterType*) PyArray_DATA(dataY), \
                                      pfb, \
                                      (double*) PyArray_DATA(dataF))
#define LAUNCH_PFB_COMPLEX(IterType) \
        compute_stokes_complex<IterType>(nStand, nSamps, nFFT, nChan, nTap, Overlap, Clip, \
                                         (IterType*) PyArray_DATA(dataX), \
                                         (IterType*) PyArray_DATA(dataY), \
                                         pfb, \
                                         (double*) PyArray_DATA(dataF))
    
    switch( PyArray_LSL_TYPE(dataX, 2) ){
        case( NPY_INT8       ): LAUNCH_PFB_REAL(int8_t);    break;
        case( NPY_INT16      ): LAUNCH_PFB_REAL(int16_t);   break;
        case( NPY_INT32      ): LAUNCH_PFB_REAL(int);       break;
        case( NPY_INT64      ): LAUNCH_PFB_REAL(long);      break;
        case( NPY_FLOAT32    ): LAUNCH_PFB_REAL(float);     break;
        case( NPY_FLOAT64    ): LAUNCH_PFB_REAL(double);    break;
        case( LSL_CI8        ): LAUNCH_PFB_COMPLEX(int8_t); break;
        case( NPY_COMPLEX64  ): LAUNCH_PFB_COMPLEX(float);  break;
        case( NPY_COMPLEX128 ): LAUNCH_PFB_COMPLEX(double); break;
        default: PyErr_Format(PyExc_RuntimeError, "Unsupport input data type"); goto fail;
    }
        
#undef LAUNCH_PFB_REAL
#undef LAUNCH_PFB_COMPLEX
    
    aligned64_free(pfb);
    
    signalsF = Py_BuildValue("O", PyArray_Return(dataF));
    
    Py_XDECREF(dataX);
    Py_XDECREF(dataY);
    Py_XDECREF(dataF);
    
    return signalsF;
    
fail:
    if( pfb != NULL ) {
        aligned64_free(pfb);
    }
    Py_XDECREF(dataX);
    Py_XDECREF(dataY);
    Py_XDECREF(dataF);
    
    return NULL;
}

PyDoc_STRVAR(PFBPSD_doc, \
"Perform a series of polyphase filter bank transforms (4-tap plus a\n\
Hanning window) on data to get the PSD for the four Stokes parameters: \n\
I, Q, U, and V.\n\
\n\
Input arguments are:\n\
 * signals: 2-D numpy.int16 (stands by samples) array of data to FFT\n\
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
 * psd: 3-D numpy.double (Stokes parameter (I,Q,U,V) by stands by channels)\n\
        of PSD data\n\
");


/*
  Cross-Multiplication And Accumulation Function ("X Engines")
    1. XEngine3 - XMAC two collections of signals
*/


template<typename InType, typename OutType>
void compute_xengine_three(long nStand,
                           long nChan,
                           long nFFT,
                           long nBL,
                           InType const* dataX,
                           InType const* dataY,
                           unsigned char const* validX,
                           unsigned char const* validY,
                           OutType* dataA) {
    // Setup
    Py_BEGIN_ALLOW_THREADS
    
    // Mapper for baseline number to stand 1, stand 2
    long s1, s2, mapper[nBL][2];
    long k = 0;
    for(s1=0; s1<nStand; s1++) {
        for(s2=s1; s2<nStand; s2++) {
            mapper[k][0] = s1;
            mapper[k++][1] = s2;
        }
    }
    
    // Cross-multiplication and accumulation
    long bl, c, f;
    OutType tempVis1, tempVis2;
    
    // Time-domain blanking control
    long nActVis;
    
    #ifdef _OPENMP
        #pragma omp parallel default(shared) private(c, f, nActVis, tempVis1, tempVis2)
    #endif
    {
        #ifdef _OPENMP
            #pragma omp for schedule(OMP_SCHEDULER)
        #endif
        for(bl=0; bl<nBL; bl++) {
            nActVis = 0;
            for(f=0; f<nFFT; f++) {
                nActVis += (long) (*(validX+ mapper[bl][0]*nFFT + f) & *(validY + mapper[bl][1]*nFFT + f));
            }
            
            for(c=0; c<nChan; c++) {
                // I
                blas_dotc_sub(nFFT, (dataX + mapper[bl][1]*nChan*nFFT + c*nFFT), 1, (dataX + mapper[bl][0]*nChan*nFFT + c*nFFT), 1, &tempVis1);
                blas_dotc_sub(nFFT, (dataY + mapper[bl][1]*nChan*nFFT + c*nFFT), 1, (dataY + mapper[bl][0]*nChan*nFFT + c*nFFT), 1, &tempVis2);
                *(dataA + 0*nBL*nChan + bl*nChan + c) = (tempVis1 + tempVis2) / (float) nActVis;
                
                // Q
                *(dataA + 1*nBL*nChan + bl*nChan + c) = (tempVis1 - tempVis2) / (float) nActVis;
                
                // U
                blas_dotc_sub(nFFT, (dataY + mapper[bl][1]*nChan*nFFT + c*nFFT), 1, (dataX + mapper[bl][0]*nChan*nFFT + c*nFFT), 1, &tempVis1);
                blas_dotc_sub(nFFT, (dataX + mapper[bl][0]*nChan*nFFT + c*nFFT), 1, (dataY + mapper[bl][1]*nChan*nFFT + c*nFFT), 1, &tempVis2);
                *(dataA + 2*nBL*nChan + bl*nChan + c) = (tempVis1 + tempVis2) / (float) nActVis;
                
                // V
                *(dataA + 3*nBL*nChan + bl*nChan + c) = (tempVis1 - tempVis2) / (float) nActVis / OutType(0,1);
            }
        }
    }
    
    Py_END_ALLOW_THREADS
    
}


static PyObject *XEngine3(PyObject *self, PyObject *args) {
    PyObject *signalsX, *signalsY, *sigValidX, *sigValidY, *output;
    PyArrayObject *dataX=NULL, *dataY=NULL, *validX=NULL, *validY=NULL, *vis=NULL;
    long nStand, nChan, nFFT, nBL;	
    
    if(!PyArg_ParseTuple(args, "OOOO", &signalsX, &signalsY, &sigValidX, &sigValidY)) {
        PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
        goto fail;
    }
    
    // Bring the data into C and make it usable
    dataX = (PyArrayObject *) PyArray_ContiguousFromObject(signalsX, 
                                                        PyArray_TYPE((PyArrayObject *) signalsX), 
                                                        3, 3);
    dataY = (PyArrayObject *) PyArray_ContiguousFromObject(signalsY, 
                                                        PyArray_TYPE((PyArrayObject *) signalsX), 
                                                        3, 3);
    validX = (PyArrayObject *) PyArray_ContiguousFromObject(sigValidX, NPY_UINT8, 2, 2);
    validY = (PyArrayObject *) PyArray_ContiguousFromObject(sigValidY, NPY_UINT8, 2, 2);
    if( dataX == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input signalsX array as a 3-D array");
        goto fail;
    }
    if( dataY == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input signalsY array as a 3-D array");
        goto fail;
    }
    if( validX == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input sigValidX array to 2-D uint8");
        goto fail;
    }
    if( validY == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input sigValidY array to 2-D uint8");
        goto fail;
    }
    
    // Get channel count and number of FFTs stored
    nStand = (long) PyArray_DIM(dataX, 0);
    nChan = (long) PyArray_DIM(dataX, 1);
    nFFT = (long) PyArray_DIM(dataX, 2);
    nBL = (nStand+1)*nStand/2;
    
    // Create the output visibility array and fill with zeros
    npy_intp dims[3];
    dims[0] = (npy_intp) 4;
    dims[1] = (npy_intp) nBL;
    dims[2] = (npy_intp) nChan;
    vis = (PyArrayObject*) PyArray_ZEROS(3, dims, NPY_COMPLEX64, 0);
    if(vis == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array");
        goto fail;
    }
    
    #define LAUNCH_XENGINE_THREE(IterType) \
        compute_xengine_three<IterType>(nStand, nChan, nFFT, nBL, \
                                        (IterType *) PyArray_DATA(dataX), \
                                        (IterType *) PyArray_DATA(dataY), \
                                        (unsigned char *) PyArray_DATA(validX), \
                                        (unsigned char *) PyArray_DATA(validY), \
                                        (Complex32 *) PyArray_DATA(vis))
    
    switch( PyArray_TYPE(dataX) ){
        case( NPY_COMPLEX64  ): LAUNCH_XENGINE_THREE(Complex32); break;
        case( NPY_COMPLEX128 ): LAUNCH_XENGINE_THREE(Complex64); break;
        default: PyErr_Format(PyExc_RuntimeError, "Unsupport input data type"); goto fail;
    }
    
    #undef LAUNCH_XENGINE_THREE
    
    output = Py_BuildValue("O", PyArray_Return(vis));
    
    Py_XDECREF(dataX);
    Py_XDECREF(dataY);
    Py_XDECREF(validX);
    Py_XDECREF(validY);
    Py_XDECREF(vis);
    
    return output;
    
fail:
    Py_XDECREF(dataX);
    Py_XDECREF(dataY);
    Py_XDECREF(validX);
    Py_XDECREF(validY);
    Py_XDECREF(vis);
    
    return NULL;
}

PyDoc_STRVAR(XEngine3_doc, \
"Perform all XMACs for a data stream out of the F engine using OpenMP that\n\
creates the four Stokes parameters: I, Q, U, and V.\n\
\n\
Input arguments are:\n\
 * fsignalsX: 3-D numpy.complex64 (stand by channels by FFT_set) array of FFTd\n\
              data from an F engine.\n\
 * fsignalsY: 3-D numpy.complex64 (stand by channels by FFT_set) array of FFTd\n\
              data from an F engine.\n\
 * sigValidX: 1-D numpy.uint8 (FFT_set) array of whether or not the FFT_set is\n\
              valid (1) or not (0) for the first signal.\n\
 * sigValidY: 1-D numpy.uint8 (FFT_set) array of whether or not the FFT_set is\n\
              valid (1) or not (0) for the second signal.\n\
\n\
Ouputs:\n\
  * visibility: 3-D numpy.cdouble (Stokes parameter (I,Q,U,V by baseline by\n\
                channel) array of cross-correlated and averaged visibility data.\n\
");


/*
  Module Setup - Function Definitions and Documentation
*/

static PyMethodDef stokes_methods[] = {
    {"FPSD",     (PyCFunction) FPSD,     METH_VARARGS|METH_KEYWORDS, FPSD_doc    }, 
    {"PFBPSD",   (PyCFunction) PFBPSD,   METH_VARARGS|METH_KEYWORDS, PFBPSD_doc  }, 
    {"XEngine3", (PyCFunction) XEngine3, METH_VARARGS,               XEngine3_doc}, 
    {NULL,       NULL,                   0,                          NULL        }
};

PyDoc_STRVAR(stokes_doc, \
"Extension to take X and Y timeseries data and create the four Stokes\n\
parameters.\n\
\n\
The functions defined in this module are:\n\
 * FPSDR -    FFT and integrate function for computing a series of overlapped\n\
              Fourier transforms for a real-valued (TBW) or complex-valued (TBN\n\
              and DRX) signals from a collection of stands all at once.\n\
 * PFBPSD -   Similar to FPSD, but using a 4-tap + Hanning windowed polyphase\n\
              filter bank.\n\
 * XEngine3 - Similar to XEngine2, but works with all linear polarization products at\n\
              once to compute all four Stokes parameters.\n\
\n\
Also included is an X-Engine for use with the lsl.correlator._core module to\n\
perform cross-correlations for the stokes parameters.\n\
\n\
See the inidividual functions for more details.");


/*
  Module Setup - Initialization
*/

static int stokes_exec(PyObject *module) {
    import_array();
    
    // Version and revision information
    PyModule_AddObject(module, "__version__", PyUnicode_FromString("0.4"));
    
    // Function listings
    PyObject* all = PyList_New(0);
    PyList_Append(all, PyUnicode_FromString("FPSD"));
    PyList_Append(all, PyUnicode_FromString("PFBPSD"));
    PyList_Append(all, PyUnicode_FromString("XEngine3"));
    PyModule_AddObject(module, "__all__", all);
    
    // LSL FFTW Wisdom
    PyObject* pModule = PyImport_ImportModule("lsl.common.paths");
    if( pModule != NULL ) {
        PyObject* pDataPath = PyObject_GetAttrString(pModule, "WISDOM");
        if( pDataPath != NULL ) {
            char filename[256];
            sprintf(filename, "%s/fftwf_wisdom.txt", PyString_AsString(pDataPath));
            read_wisdom(filename, module);
        }
        Py_XDECREF(pDataPath);
    } else {
        PyErr_Warn(PyExc_RuntimeWarning, "Cannot load the LSL FFTWF wisdom");
    }
    Py_XDECREF(pModule);
    return 0;
}

static PyModuleDef_Slot stokes_slots[] = {
    {Py_mod_exec, (void *)&stokes_exec},
    {0,           NULL}
};

static PyModuleDef stokes_def = {
    PyModuleDef_HEAD_INIT,    /* m_base */
    "_stokes",                /* m_name */
    stokes_doc,               /* m_doc */
    0,                        /* m_size */
    stokes_methods,           /* m_methods */
    stokes_slots,             /* m_slots */
    NULL,                     /* m_traverse */
    NULL,                     /* m_clear */
    NULL,                     /* m_free */
};

PyMODINIT_FUNC PyInit__stokes(void) {
    return PyModuleDef_Init(&stokes_def);
}

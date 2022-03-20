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
Function to compute the interger and fractional delays for a set of inputs
*/

long compute_delay_components(long nStand, 
                              long nChan, 
                              double SampleRate, 
                              double const* delay,
                              long* fifo, 
                              double* frac) {
    long i, j;
    long fifoMax;
    double minDelay;
    
    // Find the minimum delay
#ifdef USE_TRUE_MINIMUM_DELAY
    minDelay = 1e9;
    for(i=0; i<nStand; i++) {
        for(j=0; j<nChan; j++) {
            if( *(delay + nChan*i + j) < minDelay ) {
                minDelay = *(delay + nChan*i + j);
            }
        }
    }
#else
    minDelay = 0.0;
#endif
    
    // Compute the FIFO and fractional delays
    fifoMax = 0.0;
    for(i=0; i<nStand; i++) {
        *(fifo + i) = lround( (*(delay + nChan*i + nChan/2) - minDelay) * SampleRate );
        if( *(fifo + i) > fifoMax) {
            fifoMax = *(fifo + i);
        }
        
        for(j=0; j<nChan; j++) {
            *(frac + nChan*i + j) = (*(delay + nChan*i + j) - minDelay) - (double) *(fifo + i)/SampleRate;
        }
    }
    
    return fifoMax;
}


/*
  Function to build the phase rotator
*/

template<typename OutType>
void compute_phase_rotator(long nStand, 
                           long nChan, 
                           double SampleRate, 
                           long ChanRef,
                           float norm,
                           double const* freq,
                           const long* fifo, 
                           const double* frac,
                           OutType* rot) {
    long ij, i, j;
    
    #ifdef _OPENMP
        #pragma omp parallel default(shared) private(i, j)
    #endif
    {
        #ifdef _OPENMP
            #pragma omp for schedule(OMP_SCHEDULER)
        #endif
        for(ij=0; ij<nStand*nChan; ij++) {
            i = ij / nChan;
            j = ij % nChan;
            *(rot + nChan*i + j)  = exp(TPI * *(freq + j) * *(frac + nChan*i + j));
            *(rot + nChan*i + j) *= exp(TPI * *(freq + ChanRef) / SampleRate * (double) *(fifo + i));
            *(rot + nChan*i + j) /= sqrt(norm);
        }
    }
}


/*
FFT Functions ("F-engines")
    1. FEngine - FFT a real or complex-valued collection of signals
*/

template<typename InType, typename OutType>
void compute_fengine_real(long nStand,
                          long nSamps,
                          long nFFT,
                          int nChan, 
                          int Overlap, 
                          int Clip, 
                          double SampleRate,
                          InType const* data,
                          double const* freq,
                          long const* fifo,
                          double const* frac,
                          double const* window,
                          OutType* fdomain,
                          unsigned char* valid) {
    // Setup
    long ij, i, j, k;
    
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
    
    // Pre-compute the phase rotation and scaling factor
    OutType* rot;
    rot = (OutType*) malloc(sizeof(OutType) * nStand*nChan);
    compute_phase_rotator(nStand, nChan, SampleRate, 0, 2*nChan, freq, fifo, frac, rot);
    
    #ifdef _OPENMP
        #pragma omp parallel default(shared) private(in, out, i, j, k, secStart, cleanFactor)
    #endif
    {
        in = (float*) fftwf_malloc(sizeof(float) * 2*nChan);
        out = (Complex32*) fftwf_malloc(sizeof(Complex32) * (nChan+1));
        
        #ifdef _OPENMP
            #pragma omp for schedule(OMP_SCHEDULER)
        #endif
        for(ij=0; ij<nStand*nFFT; ij++) {
            i = ij / nFFT;
            j = ij % nFFT;
            
            cleanFactor = 1.0;
            secStart = *(fifo + i) + nSamps*i + 2*nChan*j/Overlap;
            
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
                *(fdomain + nChan*nFFT*i + nFFT*k + j)  = (float) cleanFactor * out[k];
                *(fdomain + nChan*nFFT*i + nFFT*k + j) *= *(rot + nChan*i + k);
            }
            
            *(valid + nFFT*i + j) = (unsigned char) cleanFactor;
        }
        
        fftwf_free(in);
        fftwf_free(out);
    }
    free(rot);
    
    fftwf_destroy_plan(p);
    fftwf_free(inP);
    fftwf_free(outP);
    
    Py_END_ALLOW_THREADS
    
}


template<typename InType, typename OutType>
void compute_pfbengine_real(long nStand,
                            long nSamps,
                            long nFFT,
                            int nChan, 
                            int Overlap, 
                            int Clip, 
                            double SampleRate,
                            InType const* data,
                            double const* freq,
                            long const* fifo,
                            double const* frac,
                            double const* window,
                            OutType* fdomain,
                            unsigned char* valid) {
    // Setup
    long ij, i, j, k, l;
    
    Py_BEGIN_ALLOW_THREADS
    
    // Create the FFTW plan     
    float *inP, *in;                          
    Complex32 *outP, *out;
    inP = (float*) fftwf_malloc(sizeof(float) * 2*nChan*PFB_NTAP);
    outP = (Complex32*) fftwf_malloc(sizeof(Complex32) * (nChan+1)*PFB_NTAP);
    fftwf_plan p;
    int n[] = {2*nChan,};
    p = fftwf_plan_many_dft_r2c(1, n, PFB_NTAP, \
                                inP, NULL, 1, 2*nChan, \
                                reinterpret_cast<fftwf_complex*>(outP), NULL, 1, nChan+1, \
                                FFTW_ESTIMATE);
    
    // Filter bank
    float *pfb;
    pfb = (float*) malloc(sizeof(float) * 2*nChan*PFB_NTAP);
    for(i=0; i<2*nChan*PFB_NTAP; i++) {
        *(pfb + i) = sinc((i - 2.0*nChan*PFB_NTAP/2.0 + 0.5)/(2.0*nChan));
        *(pfb + i) *= hanning(2*NPY_PI*i/(2*nChan*PFB_NTAP-1));
    }
    
    // Data indexing and access
    long secStart;
    
    // Time-domain blanking control
    double cleanFactor;
    
    // Pre-compute the phase rotation and scaling factor
    OutType* rot;
    rot = (OutType*) malloc(sizeof(OutType) * nStand*nChan);
    compute_phase_rotator(nStand, nChan, SampleRate, 0, 2*nChan, freq, fifo, frac, rot);
    
    #ifdef _OPENMP
        #pragma omp parallel default(shared) private(in, out, i, j, k, l, secStart, cleanFactor)
    #endif
    {
        in = (float*) fftwf_malloc(sizeof(float) * 2*nChan*PFB_NTAP);
        out = (Complex32*) fftwf_malloc(sizeof(Complex32) * (nChan+1)*PFB_NTAP);
        
        #ifdef _OPENMP
            #pragma omp for schedule(OMP_SCHEDULER)
        #endif
        for(ij=0; ij<nStand*nFFT; ij++) {
            i = ij / nFFT;
            j = ij % nFFT;
            
            cleanFactor = 1.0;
            secStart = *(fifo + i) + nSamps*i + 2*nChan*j/Overlap;
            
            for(k=0; k<2*nChan*PFB_NTAP; k++) {
                if( secStart - 2*nChan*(PFB_NTAP-1) + k < nSamps*i ) {
                    in[k] = 0.0;
                } else {
                    in[k] = (float) *(data + secStart - 2*nChan*(PFB_NTAP-1) + k);
                }
                
                if( Clip && fabs(in[k]) >= Clip ) {
                    cleanFactor = 0.0;
                }
                
                in[k] *= *(pfb + k);
            }
            
            fftwf_execute_dft_r2c(p, \
                                  in, \
                                  reinterpret_cast<fftwf_complex*>(out));
            
            for(l=1; l<PFB_NTAP; l++) { 
                for(k=0; k<nChan; k++) {
                    out[k] += out[k+l*(nChan+1)];
                }
            }
            
            for(k=0; k<nChan; k++) {
                *(fdomain + nChan*nFFT*i + nFFT*k + j)  = (float) cleanFactor * out[k];
                *(fdomain + nChan*nFFT*i + nFFT*k + j) *= *(rot + nChan*i + k);
            }
            
            *(valid + nFFT*i + j) = (unsigned char) cleanFactor;
        }
        
        fftwf_free(in);
        fftwf_free(out);
    }
    free(rot);
    
    free(pfb);
    fftwf_destroy_plan(p);
    fftwf_free(inP);
    fftwf_free(outP);
    
    Py_END_ALLOW_THREADS
    
}


template<typename InType, typename OutType>
void compute_fengine_complex(long nStand,
                             long nSamps,
                             long nFFT,
                             int nChan, 
                             int Overlap, 
                             int Clip, 
                             double SampleRate,
                             InType const* data,
                             double const* freq,
                             long const* fifo,
                             double const* frac,
                             double const* window,
                             OutType* fdomain,
                             unsigned char* valid) {
    // Setup
    long ij, i, j, k;
    
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
    
    // Time-domain blanking control
    double cleanFactor;
    
    // Pre-compute the phase rotation and scaling factor
    OutType* rot;
    rot = (OutType*) malloc(sizeof(OutType) * nStand*nChan);
    compute_phase_rotator(nStand, nChan, SampleRate, nChan/2, nChan, freq, fifo, frac, rot);
    
    #ifdef _OPENMP
        #pragma omp parallel default(shared) private(in, i, j, k, secStart, cleanFactor)
    #endif
    {
        in = (Complex32*) fftwf_malloc(sizeof(Complex32) * nChan);
        
        #ifdef _OPENMP
            #pragma omp for schedule(OMP_SCHEDULER)
        #endif
        for(ij=0; ij<nStand*nFFT; ij++) {
            i = ij / nFFT;
            j = ij % nFFT;
            
            cleanFactor = 1.0;
            secStart = *(fifo + i) + nSamps*i + nChan*j/Overlap;
            
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
            
            for(k=0; k<nChan/2; k++) {
                *(fdomain + nChan*nFFT*i + nFFT*k + j)  = (float) cleanFactor * in[k+nChan/2+nChan%2];
                *(fdomain + nChan*nFFT*i + nFFT*k + j) *= *(rot + nChan*i + k);
            }
            for(k=nChan/2; k<nChan; k++) {
                *(fdomain + nChan*nFFT*i + nFFT*k + j)  = (float) cleanFactor * in[k-nChan/2];
                *(fdomain + nChan*nFFT*i + nFFT*k + j) *= *(rot + nChan*i + k);
            }
            
            *(valid + nFFT*i + j) = (unsigned char) cleanFactor;
        }
        
        fftwf_free(in);
    }
    free(rot);
    
    fftwf_destroy_plan(p);
    fftwf_free(inP);
    
    Py_END_ALLOW_THREADS
}


template<typename InType, typename OutType>
void compute_pfbengine_complex(long nStand,
                               long nSamps,
                               long nFFT,
                               int nChan, 
                               int Overlap, 
                               int Clip, 
                               double SampleRate,
                               InType const* data,
                               double const* freq,
                               long const* fifo,
                               double const* frac,
                               double const* window,
                               OutType* fdomain,
                               unsigned char* valid) {
    // Setup
    long ij, i, j, k, l;
    
    Py_BEGIN_ALLOW_THREADS
    
    // Create the FFTW plan
    Complex32 *inP, *in;
    inP = (Complex32*) fftwf_malloc(sizeof(Complex32) * nChan*PFB_NTAP);
    fftwf_plan p;
    int n[] = {nChan,};
    p = fftwf_plan_many_dft(1, n, PFB_NTAP, \
                          reinterpret_cast<fftwf_complex*>(inP), NULL, 1, nChan, \
                          reinterpret_cast<fftwf_complex*>(inP), NULL, 1, nChan, \
                          FFTW_FORWARD, FFTW_ESTIMATE);
    
    // Filter bank
    float *pfb;
    pfb = (float*) malloc(sizeof(float) * nChan*PFB_NTAP);
    for(i=0; i<nChan*PFB_NTAP; i++) {
        *(pfb + i) = sinc((i - nChan*PFB_NTAP/2.0 + 0.5)/nChan);
        *(pfb + i) *= hanning(2*NPY_PI*i/(nChan*PFB_NTAP-1));
    }
    
    // Data indexing and access
    long secStart;
    
    // Time-domain blanking control
    double cleanFactor;
    
    // Pre-compute the phase rotation and scaling factor
    OutType* rot;
    rot = (OutType*) malloc(sizeof(OutType) * nStand*nChan);
    compute_phase_rotator(nStand, nChan, SampleRate, nChan/2, nChan, freq, fifo, frac, rot);
    
    #ifdef _OPENMP
        #pragma omp parallel default(shared) private(in, i, j, k, l, secStart, cleanFactor)
    #endif
    {
        in = (Complex32*) fftwf_malloc(sizeof(Complex32) * nChan*PFB_NTAP);
        
        #ifdef _OPENMP
            #pragma omp for schedule(OMP_SCHEDULER)
        #endif
        for(ij=0; ij<nStand*nFFT; ij++) {
            i = ij / nFFT;
            j = ij % nFFT;
            
            cleanFactor = 1.0;
            secStart = *(fifo + i) + nSamps*i + nChan*j/Overlap;
            
            for(k=0; k<nChan*PFB_NTAP; k++) {
                if( secStart - nChan*(PFB_NTAP-1) + k < nSamps*i ) {
                    in[k] = 0.0;
                } else {
                    in[k] = Complex32(*(data + 2*secStart - 2*nChan*(PFB_NTAP-1) + 2*k + 0), \
                                      *(data + 2*secStart - 2*nChan*(PFB_NTAP-1) + 2*k + 1));
                }
                
                if( Clip && abs(in[k]) >= Clip ) {
                    cleanFactor = 0.0;
                }
                
                in[k] *= *(pfb + k);
            }
            
            fftwf_execute_dft(p, \
                              reinterpret_cast<fftwf_complex*>(in), \
                              reinterpret_cast<fftwf_complex*>(in));
            
            for(l=1; l<PFB_NTAP; l++) {
                for(k=0; k<nChan; k++) {
                    in[k] += in[k+l*nChan];
                }
            }
            
            for(k=0; k<nChan/2; k++) {
                *(fdomain + nChan*nFFT*i + nFFT*k + j)  = (float) cleanFactor * in[k+nChan/2+nChan%2];
                *(fdomain + nChan*nFFT*i + nFFT*k + j) *= *(rot + nChan*i + k);
            }
            for(k=nChan/2; k<nChan; k++) {
                *(fdomain + nChan*nFFT*i + nFFT*k + j)  = (float) cleanFactor * in[k-nChan/2];
                *(fdomain + nChan*nFFT*i + nFFT*k + j) *= *(rot + nChan*i + k);
            }
            
            *(valid + nFFT*i + j) = (unsigned char) cleanFactor;
        }
        
        fftwf_free(in);
    }
    free(rot);
    
    free(pfb);
    fftwf_destroy_plan(p);
    fftwf_free(inP);
    
    Py_END_ALLOW_THREADS
}


static PyObject *FEngine(PyObject *self, PyObject *args, PyObject *kwds) {
    PyObject *signals, *freqs, *delays, *window=Py_None, *signalsF;
    PyArrayObject *data=NULL, *freq=NULL, *delay=NULL, *dataF=NULL, *validF=NULL, *windowData=NULL;
    int isReal;
    int nChan = 64;
    int Overlap = 1;
    int Clip = 0;
    double SampleRate = 196.0e6;

    long nStand, nSamps, nFFT;
    
    char const* kwlist[] = {"signals", "freqs", "delays", "LFFT", "overlap", "sample_rate", "clip_level", "window", NULL};
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "OOO|iidiO:set_callback", const_cast<char **>(kwlist), &signals, &freqs, &delays, &nChan, &Overlap, &SampleRate, &Clip, &window)) {
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
    freq = (PyArrayObject *) PyArray_ContiguousFromObject(freqs, NPY_DOUBLE, 1, 1);
    delay = (PyArrayObject *) PyArray_ContiguousFromObject(delays, NPY_DOUBLE, 2, 2);
    if( data == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input signals array as a 2-D array");
        goto fail;
    }
    if( freq == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input freq array to 1-D double");
        goto fail;
    }
    if( delay == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input delay array to 2-D double");
        goto fail;
    }
    
    // Check data dimensions
    if(PyArray_DIM(data, 0) != PyArray_DIM(delay, 0)) {
        PyErr_Format(PyExc_RuntimeError, "signals and delays have different stand counts");
        goto fail;
    }
    
    if(nChan != PyArray_DIM(freq, 0)) {
        PyErr_Format(PyExc_RuntimeError, "freqs has a different channel count than nChan");
        goto fail;
    }
    
    if(PyArray_DIM(freq, 0) != PyArray_DIM(delay, 1)) {
        PyErr_Format(PyExc_RuntimeError, "freqs and delays have different channel counts");
        goto fail;
    }
    
    // Get the properties of the data
    nStand = (long) PyArray_DIM(data, 0);
    nSamps = (long) PyArray_DIM(data, 1);
    isReal = 1 - PyArray_ISCOMPLEX(data);
    
    // Calculate the windowing function
    if( windowFunc != Py_None ) {
        window = Py_BuildValue("(i)", (1+isReal)*nChan);
        window = PyObject_CallObject(windowFunc, window);
        windowData = (PyArrayObject *) PyArray_ContiguousFromObject(window, NPY_DOUBLE, 1, 1);
        Py_DECREF(window);
    }
    
    // Compute the integer sample offset and the fractional sample delay for each stand
    long *fifo, fifoMax;
    double *frac;
    fifo = (long *) malloc(nStand*sizeof(long));
    frac = (double *) malloc(nStand*nChan*sizeof(double));
    if( fifo == NULL || frac == NULL ) {
        PyErr_Format(PyExc_MemoryError, "Cannot create fifo/fractional delay arrays");
        goto fail;
    }
    fifoMax = compute_delay_components(nStand, nChan, SampleRate, \
                                       (double*) PyArray_DATA(delay), fifo, frac);
    
    // Find out how large the output array needs to be and initialize it
    nFFT = (nSamps - fifoMax) / (((1+isReal)*nChan)/Overlap) - ((1+isReal)*nChan)/(((1+isReal)*nChan)/Overlap) + 1;
    npy_intp dims[3];
    dims[0] = (npy_intp) nStand;
    dims[1] = (npy_intp) nChan;
    dims[2] = (npy_intp) nFFT;
    dataF = (PyArrayObject*) PyArray_ZEROS(3, dims, NPY_COMPLEX64, 0);
    if(dataF == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array");
        free(fifo);
        free(frac);
        goto fail;
    }
    
    // Create an array to store whether or not the FFT window is valid (1) or not (0)
    npy_intp dimsV[2];
    dimsV[0] = (npy_intp) nStand;
    dimsV[1] = (npy_intp) nFFT;
    validF = (PyArrayObject*) PyArray_ZEROS(2, dimsV, NPY_UINT8, 0);
    if(validF == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create valid index array");
        free(fifo);
        free(frac);
        goto fail;
    }
    
#define LAUNCH_FENGINE_REAL(IterType) \
        compute_fengine_real<IterType>(nStand, nSamps, nFFT, nChan, Overlap, Clip, SampleRate, \
                                       (IterType*) PyArray_DATA(data), \
                                       (double*) PyArray_DATA(freq), \
                                       (long*) fifo, (double*) frac, \
                                       (double*) PyArray_SAFE_DATA(windowData), \
                                       (Complex32*) PyArray_DATA(dataF), \
                                       (unsigned char*) PyArray_DATA(validF))
#define LAUNCH_FENGINE_COMPLEX(IterType) \
        compute_fengine_complex<IterType>(nStand, nSamps, nFFT, nChan, Overlap, Clip, SampleRate, \
                                          (IterType*) PyArray_DATA(data), \
                                          (double*) PyArray_DATA(freq), \
                                          (long*) fifo, (double*) frac, \
                                          (double*) PyArray_SAFE_DATA(windowData), \
                                          (Complex32*) PyArray_DATA(dataF), \
                                          (unsigned char*) PyArray_DATA(validF))
    
    switch( PyArray_TYPE(data) ){
        case( NPY_INT8       ): LAUNCH_FENGINE_REAL(int8_t);    break;
        case( NPY_INT16      ): LAUNCH_FENGINE_REAL(int16_t);   break;
        case( NPY_INT32      ): LAUNCH_FENGINE_REAL(int);       break;
        case( NPY_INT64      ): LAUNCH_FENGINE_REAL(long);      break;
        case( NPY_FLOAT32    ): LAUNCH_FENGINE_REAL(float);     break;
        case( NPY_FLOAT64    ): LAUNCH_FENGINE_REAL(double);    break;
        case( NPY_COMPLEX64  ): LAUNCH_FENGINE_COMPLEX(float);  break;
        case( NPY_COMPLEX128 ): LAUNCH_FENGINE_COMPLEX(double); break;
        default: PyErr_Format(PyExc_RuntimeError, "Unsupport input data type"); goto fail;
    }
        
#undef LAUNCH_FENGINE_REAL
#undef LAUNCH_FENGINE_COMPLEX
    
    free(frac);
    free(fifo);
    
    signalsF = Py_BuildValue("(OO)", PyArray_Return(dataF), PyArray_Return(validF));
    
    Py_XDECREF(data);
    Py_XDECREF(freq);
    Py_XDECREF(delay);
    Py_XDECREF(windowData);
    Py_XDECREF(dataF);
    Py_XDECREF(validF);
    
    return signalsF;
    
fail:
    Py_XDECREF(data);
    Py_XDECREF(freq);
    Py_XDECREF(delay);
    Py_XDECREF(windowData);
    Py_XDECREF(dataF);
    Py_XDECREF(validF);
    
    return NULL;
}

PyDoc_STRVAR(FEngine_doc, \
"Perform a series of overlapped Fourier transforms on real-valued data using\n\
OpenMP and windows.\n\
\n\
Input arguments are:\n\
 * signals: 2-D numpy (stands by samples) array of data to FFT\n\
 * frequency: 1-D numpy.double array of frequency values in Hz for the\n\
              FFT channels\n\
 * delays: 1-D numpy.double array of delays to apply to each stand\n\
\n\
Input keywords are:\n\
 * LFFT: number of FFT channels to make (default=64)\n\
 * overlap: number of overlapped FFTs to use (default=1)\n\
 * sample_rate: sample rate of the data (default=196e6)\n\
 * window: Callable Python function for generating the window or None for no\n\
           window\n\
 * clip_level: count value of 'bad' data.  FFT windows with instantaneous\n\
               powers greater than or equal to this value greater are zeroed.  \n\
               Setting the ClipLeve to zero disables time-domain blanking\n\
\n\
Outputs:\n\
 * fsignals: 3-D numpy.complex64 (stands by channels by FFT_set) of FFTd\n\
             data\n\
 * valid: 2-D numpy.uint8 (stands by FFT_set) of whether or not the FFT\n\
          set is valid (1) or not (0)\n\
");


static PyObject *PFBEngine(PyObject *self, PyObject *args, PyObject *kwds) {
    PyObject *signals, *freqs, *delays, *window=Py_None, *signalsF;
    PyArrayObject *data=NULL, *freq=NULL, *delay=NULL, *dataF=NULL, *validF=NULL, *windowData=NULL;
    int isReal;
    int nChan = 64;
    int Overlap = 1;
    int Clip = 0;
    double SampleRate = 196.0e6;

    long nStand, nSamps, nFFT;
    
    char const* kwlist[] = {"signals", "freqs", "delays", "LFFT", "overlap", "sample_rate", "clip_level", "window", NULL};
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "OOO|iidiO:set_callback", const_cast<char **>(kwlist), &signals, &freqs, &delays, &nChan, &Overlap, &SampleRate, &Clip, &window)) {
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
    freq = (PyArrayObject *) PyArray_ContiguousFromObject(freqs, NPY_DOUBLE, 1, 1);
    delay = (PyArrayObject *) PyArray_ContiguousFromObject(delays, NPY_DOUBLE, 2, 2);
    if( data == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input signals array as a 2-D array");
        goto fail;
    }
    if( freq == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input freq array to 1-D double");
        goto fail;
    }
    if( delay == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input delay array to 2-D double");
        goto fail;
    }
    
    // Check data dimensions
    if(PyArray_DIM(data, 0) != PyArray_DIM(delay, 0)) {
        PyErr_Format(PyExc_RuntimeError, "signals and delays have different stand counts");
        goto fail;
    }
    
    if(nChan != PyArray_DIM(freq, 0)) {
        PyErr_Format(PyExc_RuntimeError, "freqs has a different channel count than nChan");
        goto fail;
    }
    
    if(PyArray_DIM(freq, 0) != PyArray_DIM(delay, 1)) {
        PyErr_Format(PyExc_RuntimeError, "freqs and delays have different channel counts");
        goto fail;
    }
    
    // Get the properties of the data
    nStand = (long) PyArray_DIM(data, 0);
    nSamps = (long) PyArray_DIM(data, 1);
    isReal = 1 - PyArray_ISCOMPLEX(data);
    
    // Calculate the windowing function
    if( windowFunc != Py_None ) {
        window = Py_BuildValue("(i)", (1+isReal)*nChan);
        window = PyObject_CallObject(windowFunc, window);
        windowData = (PyArrayObject *) PyArray_ContiguousFromObject(window, NPY_DOUBLE, 1, 1);
        Py_DECREF(window);
    }
    
    // Compute the integer sample offset and the fractional sample delay for each stand
    long *fifo, fifoMax;
    double *frac;
    fifo = (long *) malloc(nStand*sizeof(long));
    frac = (double *) malloc(nStand*nChan*sizeof(double));
    if( fifo == NULL || frac == NULL ) {
        PyErr_Format(PyExc_MemoryError, "Cannot create fifo/fractional delay arrays");
        goto fail;
    }
    fifoMax = compute_delay_components(nStand, nChan, SampleRate, \
                                       (double*) PyArray_DATA(delay), fifo, frac);
    
    // Find out how large the output array needs to be and initialize it
    nFFT = (nSamps - fifoMax) / (((1+isReal)*nChan)/Overlap) - ((1+isReal)*nChan)/(((1+isReal)*nChan)/Overlap) + 1;
    npy_intp dims[3];
    dims[0] = (npy_intp) nStand;
    dims[1] = (npy_intp) nChan;
    dims[2] = (npy_intp) nFFT;
    dataF = (PyArrayObject*) PyArray_ZEROS(3, dims, NPY_COMPLEX64, 0);
    if(dataF == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array");
        free(fifo);
        free(frac);
        goto fail;
    }
    
    // Create an array to store whether or not the FFT window is valid (1) or not (0)
    npy_intp dimsV[2];
    dimsV[0] = (npy_intp) nStand;
    dimsV[1] = (npy_intp) nFFT;
    validF = (PyArrayObject*) PyArray_ZEROS(2, dimsV, NPY_UINT8, 0);
    if(validF == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create valid index array");
        free(fifo);
        free(frac);
        goto fail;
    }
    
#define LAUNCH_PFBENGINE_REAL(IterType) \
        compute_pfbengine_real<IterType>(nStand, nSamps, nFFT, nChan, Overlap, Clip, SampleRate, \
                                         (IterType*) PyArray_DATA(data), \
                                         (double*) PyArray_DATA(freq), \
                                         (long*) fifo, (double*) frac, \
                                         (double*) PyArray_SAFE_DATA(windowData), \
                                         (Complex32*) PyArray_DATA(dataF), \
                                         (unsigned char*) PyArray_DATA(validF))
#define LAUNCH_PFBENGINE_COMPLEX(IterType) \
        compute_pfbengine_complex<IterType>(nStand, nSamps, nFFT, nChan, Overlap, Clip, SampleRate, \
                                            (IterType*) PyArray_DATA(data), \
                                            (double*) PyArray_DATA(freq), \
                                            (long*) fifo, (double*) frac, \
                                            (double*) PyArray_SAFE_DATA(windowData), \
                                            (Complex32*) PyArray_DATA(dataF), \
                                            (unsigned char*) PyArray_DATA(validF))
    
    switch( PyArray_TYPE(data) ){
        case( NPY_INT8       ): LAUNCH_PFBENGINE_REAL(int8_t);    break;
        case( NPY_INT16      ): LAUNCH_PFBENGINE_REAL(int16_t);   break;
        case( NPY_INT32      ): LAUNCH_PFBENGINE_REAL(int);       break;
        case( NPY_INT64      ): LAUNCH_PFBENGINE_REAL(long);      break;
        case( NPY_FLOAT32    ): LAUNCH_PFBENGINE_REAL(float);     break;
        case( NPY_FLOAT64    ): LAUNCH_PFBENGINE_REAL(double);    break;
        case( NPY_COMPLEX64  ): LAUNCH_PFBENGINE_COMPLEX(float);  break;
        case( NPY_COMPLEX128 ): LAUNCH_PFBENGINE_COMPLEX(double); break;
        default: PyErr_Format(PyExc_RuntimeError, "Unsupport input data type"); goto fail;
    }
        
#undef LAUNCH_PFBENGINE_REAL
#undef LAUNCH_PFBENGINE_COMPLEX
    
    free(frac);
    free(fifo);
    
    signalsF = Py_BuildValue("(OO)", PyArray_Return(dataF), PyArray_Return(validF));
    
    Py_XDECREF(data);
    Py_XDECREF(freq);
    Py_XDECREF(delay);
    Py_XDECREF(windowData);
    Py_XDECREF(dataF);
    Py_XDECREF(validF);
    
    return signalsF;
    
fail:
    Py_XDECREF(data);
    Py_XDECREF(freq);
    Py_XDECREF(delay);
    Py_XDECREF(windowData);
    Py_XDECREF(dataF);
    Py_XDECREF(validF);
    
    return NULL;
}

PyDoc_STRVAR(PFBEngine_doc, \
"Perform a series of overlapped polyphase filter bank transforms (4-tap plus a\n\
Hanning window) on real-valued data using\n\
OpenMP and windows.\n\
\n\
Input arguments are:\n\
 * signals: 2-D numpy (stands by samples) array of data to FFT\n\
 * frequency: 1-D numpy.double array of frequency values in Hz for the\n\
              FFT channels\n\
 * delays: 1-D numpy.double array of delays to apply to each stand\n\
\n\
Input keywords are:\n\
 * LFFT: number of FFT channels to make (default=64)\n\
 * overlap: number of overlapped FFTs to use (default=1)\n\
 * sample_rate: sample rate of the data (default=196e6)\n\
 * window: Callable Python function for generating the window or None for no\n\
           window\n\
 * clip_level: count value of 'bad' data.  FFT windows with instantaneous\n\
               powers greater than or equal to this value greater are zeroed.  \n\
               Setting the ClipLeve to zero disables time-domain blanking\n\
\n\
Outputs:\n\
 * fsignals: 3-D numpy.complex64 (stands by channels by FFT_set) of FFTd\n\
             data\n\
 * valid: 2-D numpy.uint8 (stands by FFT_set) of whether or not the FFT\n\
          set is valid (1) or not (0)\n\
");


/*
Cross-Multiplication And Accumulation Function ("X Engines")
    1. XEngine2 - XMAC two collections of signals
*/

template<typename InType, typename OutType>
void compute_xengine_two(long nStand,
                         long nChan,
                         long nFFT,
                         long nBL,
                         InType const* data1,
                         InType const* data2,
                         unsigned char const* valid1,
                         unsigned char const* valid2,
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
    OutType tempVis;
    
    // Time-domain blanking control
    long nActVis;
    
    #ifdef _OPENMP
        #pragma omp parallel default(shared) private(c, f, nActVis, tempVis)
    #endif
    {
        #ifdef _OPENMP
            #pragma omp for schedule(OMP_SCHEDULER)
        #endif
        for(bl=0; bl<nBL; bl++) {
            nActVis = 0;
            for(f=0; f<nFFT; f++) {
                nActVis += (long) (*(valid1 + mapper[bl][0]*nFFT + f) & *(valid2 + mapper[bl][1]*nFFT + f));
            }
            
            for(c=0; c<nChan; c++) {
                blas_dotc_sub(nFFT, (data2 + mapper[bl][1]*nChan*nFFT + c*nFFT), 1, (data1 + mapper[bl][0]*nChan*nFFT + c*nFFT), 1, &tempVis);
                *(dataA + bl*nChan + c) = tempVis / (float) nActVis;
            }
        }
    }
    
    Py_END_ALLOW_THREADS
    
}


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
    OutType tempVis;
    
    // Time-domain blanking control
    long nActVisPureX, nActVisPureY, nActVisCross0, nActVisCross1;
    
    #ifdef _OPENMP
        #pragma omp parallel default(shared) private(c, f, nActVisPureX, nActVisPureY, nActVisCross0, nActVisCross1, tempVis)
    #endif
    {
        #ifdef _OPENMP
            #pragma omp for schedule(OMP_SCHEDULER)
        #endif
        for(bl=0; bl<nBL; bl++) {
            nActVisPureX = 0;
            nActVisPureY = 0;
            nActVisCross0 = 0;
            nActVisCross1 = 0;
            for(f=0; f<nFFT; f++) {
                nActVisPureX  += (long) (*(validX + mapper[bl][0]*nFFT + f) & *(validX + mapper[bl][1]*nFFT + f));
                nActVisPureY  += (long) (*(validY + mapper[bl][0]*nFFT + f) & *(validY + mapper[bl][1]*nFFT + f));
                nActVisCross0 += (long) (*(validX + mapper[bl][0]*nFFT + f) & *(validY + mapper[bl][1]*nFFT + f));
                nActVisCross1 += (long) (*(validY + mapper[bl][0]*nFFT + f) & *(validX + mapper[bl][1]*nFFT + f));
            }
            
            for(c=0; c<nChan; c++) {
                // XX
                blas_dotc_sub(nFFT, (dataX + mapper[bl][1]*nChan*nFFT + c*nFFT), 1, (dataX + mapper[bl][0]*nChan*nFFT + c*nFFT), 1, &tempVis);
                *(dataA + 0*nBL*nChan + bl*nChan + c) = tempVis / (float) nActVisPureX;
                
                // XY
                blas_dotc_sub(nFFT, (dataY + mapper[bl][1]*nChan*nFFT + c*nFFT), 1, (dataX + mapper[bl][0]*nChan*nFFT + c*nFFT), 1, &tempVis);
                *(dataA + 1*nBL*nChan + bl*nChan + c) = tempVis / (float) nActVisCross0;
                
                // YX
                blas_dotc_sub(nFFT, (dataX + mapper[bl][1]*nChan*nFFT + c*nFFT), 1, (dataY + mapper[bl][0]*nChan*nFFT + c*nFFT), 1, &tempVis);
                *(dataA + 2*nBL*nChan + bl*nChan + c) = tempVis / (float) nActVisCross1;
                
                // YY
                blas_dotc_sub(nFFT, (dataY + mapper[bl][1]*nChan*nFFT + c*nFFT), 1, (dataY + mapper[bl][0]*nChan*nFFT + c*nFFT), 1, &tempVis);
                *(dataA + 3*nBL*nChan + bl*nChan + c) = tempVis / (float) nActVisPureY;
            }
        }
    }
    
    Py_END_ALLOW_THREADS
    
}


static PyObject *XEngine2(PyObject *self, PyObject *args) {
    PyObject *signals1, *signals2, *sigValid1, *sigValid2, *output;
    PyArrayObject *data1=NULL, *data2=NULL, *valid1=NULL, *valid2=NULL, *vis=NULL;
    long nStand, nChan, nFFT, nBL;	

    if(!PyArg_ParseTuple(args, "OOOO", &signals1, &signals2, &sigValid1, &sigValid2)) {
        PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
        goto fail;
    }

    // Bring the data into C and make it usable
    data1 = (PyArrayObject *) PyArray_ContiguousFromObject(signals1, 
                                                        PyArray_TYPE((PyArrayObject *) signals1), 
                                                        3, 3);
    data2 = (PyArrayObject *) PyArray_ContiguousFromObject(signals2, 
                                                        PyArray_TYPE((PyArrayObject *) signals1), 
                                                        3, 3);
    valid1 = (PyArrayObject *) PyArray_ContiguousFromObject(sigValid1, NPY_UINT8, 2, 2);
    valid2 = (PyArrayObject *) PyArray_ContiguousFromObject(sigValid2, NPY_UINT8, 2, 2);
    if( data1 == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input signals1 array as a 3-D array");
        goto fail;
    }
    if( data2 == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input signals2 array as a 3-D array");
        goto fail;
    }
    if( valid1 == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input sigValid1 array to 2-D uint8");
        goto fail;
    }
    if( valid2 == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input sigValid2 array to 2-D uint8");
        goto fail;
    }
    
    // Get channel count and number of FFTs stored
    nStand = (long) PyArray_DIM(data1, 0);
    nChan = (long) PyArray_DIM(data1, 1);
    nFFT = (long) PyArray_DIM(data1, 2);
    nBL = (nStand+1)*nStand/2;
    
    // Create the output visibility array and fill with zeros
    npy_intp dims[2];
    dims[0] = (npy_intp) nBL;
    dims[1] = (npy_intp) nChan;
    vis = (PyArrayObject*) PyArray_ZEROS(2, dims, NPY_COMPLEX64, 0);
    if(vis == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array");
        goto fail;
    }
    
    #define LAUNCH_XENGINE_TWO(IterType) \
        compute_xengine_two<IterType>(nStand, nChan, nFFT, nBL, \
                                      (IterType *) PyArray_DATA(data1), \
                                      (IterType *) PyArray_DATA(data2), \
                                      (unsigned char *) PyArray_DATA(valid1), \
                                      (unsigned char *) PyArray_DATA(valid2), \
                                      (Complex32 *) PyArray_DATA(vis))
    
    switch( PyArray_TYPE(data1) ){
        case( NPY_COMPLEX64  ): LAUNCH_XENGINE_TWO(Complex32); break;
        case( NPY_COMPLEX128 ): LAUNCH_XENGINE_TWO(Complex64); break;
        default: PyErr_Format(PyExc_RuntimeError, "Unsupport input data type"); goto fail;
    }
    
    #undef LAUNCH_XENGINE_TWO
    
    output = Py_BuildValue("O", PyArray_Return(vis));
    
    Py_XDECREF(data1);
    Py_XDECREF(data2);
    Py_XDECREF(valid1);
    Py_XDECREF(valid2);
    Py_XDECREF(vis);
    
    return output;
    
fail:
    Py_XDECREF(data1);
    Py_XDECREF(data2);
    Py_XDECREF(valid1);
    Py_XDECREF(valid2);
    Py_XDECREF(vis);
    
    return NULL;
}

PyDoc_STRVAR(XEngine2_doc, \
"Perform all XMACs for a data stream out of the F engine using OpenMP.\n\
\n\
.. versionchanged:: 0.5\n\
\tThe second signal is not longer input as a conjugated array.  Rather\n\
\tthe conjucation is performed as part of the cross-correlation.\n\
\n\
Input arguments are:\n\
 * fsignals1: 3-D numpy.complex64 (stand by channels by FFT_set) array of FFTd\n\
              data from an F engine.\n\
 * fsignals2: 3-D numpy.complex64 (stand by channels by FFT_set) array of\n\
              FFTd data from an F engine.\n\
 * sigValid1: 1-D numpy.uint8 (FFT_set) array of whether or not the FFT_set is\n\
              valid (1) or not (0) for the first signal.\n\
 * sigValid2: 1-D numpy.uint8 (FFT_set) array of whether or not the FFT_set is\n\
              valid (1) or not (0) for the second signal.\n\
\n\
Ouputs:\n\
 * visibility: 2-D numpy.complex64 (baseline by channel) array of cross-\n\
               correlated and averaged visibility data.\n\
");


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
creates the four linear polarization products\n\
\n\
.. versionadded:: 1.1.2\n\
\n\
Input arguments are:\n\
 * fsignals1: 3-D numpy.complex64 (stand by channels by FFT_set) array of FFTd\n\
              data from an F engine.\n\
 * fsignals2: 3-D numpy.complex64 (stand by channels by FFT_set) array of FFTd\n\
              data from an F engine.\n\
 * sigValid1: 1-D numpy.uint8 (FFT_set) array of whether or not the FFT_set is\n\
              valid (1) or not (0) for the first signal.\n\
 * sigValid2: 1-D numpy.uint8 (FFT_set) array of whether or not the FFT_set is\n\
              valid (1) or not (0) for the second signal.\n\
\n\
Ouputs:\n\
 * visibility: 3-D numpy.complex64 (Stokes parameter (XX,XY,YX,YY) by baseline by\n\
               channel) array of cross-correlated and averaged visibility data.\n\
");


/*
Module Setup - Function Definitions and Documentation
*/

static PyMethodDef CorrelatorMethods[] = {
    {"FEngine",   (PyCFunction) FEngine,   METH_VARARGS|METH_KEYWORDS, FEngine_doc  },  
    {"PFBEngine", (PyCFunction) PFBEngine, METH_VARARGS|METH_KEYWORDS, PFBEngine_doc}, 
    {"XEngine2",  (PyCFunction) XEngine2,  METH_VARARGS,               XEngine2_doc }, 
    {"XEngine3",  (PyCFunction) XEngine3,  METH_VARARGS,               XEngine3_doc }, 
    {NULL,        NULL,                    0,                          NULL         }
};

PyDoc_STRVAR(correlator_doc, \
"C-based F and X engines for the LWA software FX correlator.  These function\n\
are meant to provide an alternative to the lsl.correlator.fx.correlate function and \n\
provide a much-needed speed boost to cross-correlation.\n\
\n\
The function defined in this module are:\n\
 * FEngine - F-engine for computing a series of overlapped Fourier transforms with\n\
             delay corrections for a real-valued (TBW) or complex valued (TBN or DRX)\n\
             signals from a collection of stands all at once.\n\
 * PFBEngine - Similar to FEngine, but using a 4-tap + Hanning windowed polyphase\n\
               filter bank.\n\
 * XEngine2 - Cross multiply a single polarization stream created by FEngine or\n\
              PFBEngine\n\
 * XEngine3 - Similar to XEngine2, but works with all linear polarization products at\n\
              once.\n\
\n\
See the inidividual functions for more details.");


/*
Module Setup - Initialization
*/

MOD_INIT(_core) {
    char filename[256];
    PyObject *m, *all, *pModule, *pDataPath=NULL;
    
    Py_Initialize();
    
    // Module definitions and functions
    MOD_DEF(m, "_core", CorrelatorMethods, correlator_doc);
    if( m == NULL ) {
        return MOD_ERROR_VAL;
    }
    import_array();
    
    // Version and revision information
    PyModule_AddObject(m, "__version__", PyString_FromString("0.8"));
    
    // Function listings
    all = PyList_New(0);
    PyList_Append(all, PyString_FromString("FEngine"));
    PyList_Append(all, PyString_FromString("PFBEngine"));
    PyList_Append(all, PyString_FromString("XEngine2"));
    PyList_Append(all, PyString_FromString("XEngine3"));
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

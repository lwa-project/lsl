#include "Python.h"
#include <math.h>
#include <stdio.h>
#include <complex.h>
#include <cblas.h>
#include <fftw3.h>
#include <stdlib.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

#include "numpy/arrayobject.h"

#define TPI (2*NPY_PI*_Complex_I)


/*
 Load in FFTW wisdom.  Based on the read_wisdom function in PRESTO.
*/

void read_wisdom(char *filename, PyObject *m) {
	int status = 0;
	FILE *wisdomfile;
	
	wisdomfile = fopen(filename, "r");
	if( wisdomfile != NULL ) {
		status = fftwf_import_wisdom_from_file(wisdomfile);
		PyModule_AddObject(m, "useWisdom", PyBool_FromLong(status));
		fclose(wisdomfile);
	} else {
		PyModule_AddObject(m, "useWisdom", PyBool_FromLong(status));
	}
}


/*
  Holder for window function callback
*/

static PyObject *windowFunc = NULL;


/*
  Function to compute the interger and fractional delays for a set of inputs
*/

long computeDelayComponents(PyArrayObject *delays, double SampleRate, long *fifo, double *frac) {
	long i, j;
	long nStand, nChan, fifoMax;
	double minDelay;
	
	// Get the properties of the delays
	nStand = (long) PyArray_DIM(delays, 0);
	nChan = (long) PyArray_DIM(delays, 1);
	
	// Set up a way to access the data
	double *a;
	a = (double *) PyArray_DATA(delays);
	
	// Find the minimum delay
	/*
	minDelay = 1e9;
	for(i=0; i<nStand; i++) {
		for(j=0; j<nChan; j++) {
			if( *(a + nChan*i + j) < minDelay ) {
				minDelay = *(a + nChan*i + j);
			}
		}
	}
	*/
	minDelay = 0.0;
	
	// Compute the FIFO and fractional delays
	fifoMax = 0.0;
	for(i=0; i<nStand; i++) {
		*(fifo + i) = lround( (*(a + nChan*i + nChan/2) - minDelay) * SampleRate );
		if( *(fifo + i) > fifoMax) {
			fifoMax = *(fifo + i);
		}
		
		for(j=0; j<nChan; j++) {
			*(frac + nChan*i + j) = (*(a + nChan*i + j) - minDelay) - (double) *(fifo + i)/SampleRate;
		}
	}
	
	return fifoMax;
}


/*
  FFT Functions ("F-engines")
    1. FEngineR2 - FFT a real-valued collection of signals
    2. FEngineR3 - window the data and FFT a real-valued collection of signals
    3. FEngineC2 - FFT a complex-valued collection of signals
    4. FEngineC3 - window the data and FFT a complex-valued collection of signals
*/


static PyObject *FEngineR2(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *freqs, *delays, *signalsF;
	PyArrayObject *data=NULL, *freq=NULL, *delay=NULL, *dataF=NULL, *validF=NULL;
	int nChan = 64;
	int Overlap = 1;
	int Clip = 0;
	double SampleRate = 196.0e6;

	long ij, i, j, k, nStand, nSamps, nFFT;
	
	static char *kwlist[] = {"signals", "freqs", "delays", "LFFT", "Overlap", "SampleRate", "ClipLevel", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OOO|iidi", kwlist, &signals, &freqs, &delays, &nChan, &Overlap, &SampleRate, &Clip)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		goto fail;
	}

	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_INT16, 2, 2);
	freq = (PyArrayObject *) PyArray_ContiguousFromObject(freqs, NPY_DOUBLE, 1, 1);
	delay = (PyArrayObject *) PyArray_ContiguousFromObject(delays, NPY_DOUBLE, 2, 2);
	
	// Check data dimensions
	if(PyArray_DIM(data, 0) != PyArray_DIM(delay, 0)) {
		PyErr_Format(PyExc_TypeError, "signals and delays have different stand counts");
		goto fail;
	}
	
	if(nChan != PyArray_DIM(freq, 0)) {
		PyErr_Format(PyExc_RuntimeError, "freqs has a different channel count than nChan");
		goto fail;
	}
	
	if(PyArray_DIM(freq, 0) != PyArray_DIM(delay, 1)) {
		PyErr_Format(PyExc_TypeError, "freqs and delays have different channel counts");
		goto fail;
	}
	
	// Get the properties of the data
	nStand = (long) PyArray_DIM(data, 0);
	nSamps = (long) PyArray_DIM(data, 1);
	
	// Compute the integer sample offset and the fractional sample delay for each stand
	long *fifo, fifoMax;
	double *frac;
	fifo = (long *) malloc(nStand*sizeof(long));
	frac = (double *) malloc(nStand*nChan*sizeof(double));
	if( fifo == NULL || frac == NULL ) {
		PyErr_Format(PyExc_MemoryError, "Cannot create fifo/fractional delay arrays");
		goto fail;
	}
	fifoMax = computeDelayComponents(delay, SampleRate, fifo, frac);

	// Find out how large the output array needs to be and initialize it
	nFFT = (nSamps - fifoMax) / ((2*nChan)/Overlap) - (2*nChan)/((2*nChan)/Overlap) + 1;
	npy_intp dims[3];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChan;
	dims[2] = (npy_intp) nFFT;
	dataF = (PyArrayObject*) PyArray_ZEROS(3, dims, NPY_COMPLEX64, 0);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		free(fifo)
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
		free(fifo)
		free(frac);
		goto fail;
	}
	
	Py_BEGIN_ALLOW_THREADS
	
	// Create the FFTW plan                          
	float complex *inP, *in;                          
	inP = (float complex *) fftwf_malloc(sizeof(float complex) * 2*nChan);
	fftwf_plan p;
	p = fftwf_plan_dft_1d(2*nChan, inP, inP, FFTW_FORWARD, FFTW_ESTIMATE);
	
	// Data indexing and access
	long secStart;
	short int *a;
	float complex *b;
	double *c;
	unsigned char *d;
	a = (short int *) PyArray_DATA(data);
	b = (float complex *) PyArray_DATA(dataF);
	c = (double *) PyArray_DATA(freq);
	d = (unsigned char *) PyArray_DATA(validF);
	
	// Time-domain blanking control
	double cleanFactor;
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(in, i, j, k, secStart, cleanFactor)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(dynamic)
		#endif
		for(ij=0; ij<nStand*nFFT; ij++) {
			i = ij / nFFT;
			j = ij % nFFT;
			
			in = (float complex *) fftwf_malloc(sizeof(float complex) * 2*nChan);
			
			cleanFactor = 1.0;
			secStart = *(fifo + i) + nSamps*i + 2*nChan*j/Overlap;
			
			for(k=0; k<2*nChan; k++) {
				in[k] = (float complex) *(a + secStart + k);
				
				if( Clip && cabsf(in[k]) >= Clip ) {
					cleanFactor = 0.0;
				}
			}
			
			fftwf_execute_dft(p, in, in);
			
			for(k=0; k<nChan; k++) {
				*(b + nChan*nFFT*i + nFFT*k + j)  = cleanFactor*in[k];
				*(b + nChan*nFFT*i + nFFT*k + j) *= cexp(TPI * *(c + k) * *(frac + nChan*i + k));
				*(b + nChan*nFFT*i + nFFT*k + j) /= sqrt(2*nChan);
			}
			
			*(d + nFFT*i + j) = (unsigned char) cleanFactor;
			
			fftwf_free(in);
		}
	}
	fftwf_destroy_plan(p);
	fftwf_free(inP);
	free(frac);
	free(fifo);
	
	Py_END_ALLOW_THREADS
	
	signalsF = Py_BuildValue("(OO)", PyArray_Return(dataF), PyArray_Return(validF));
	
	Py_XDECREF(data);
	Py_XDECREF(freq);
	Py_XDECREF(delay);
	Py_XDECREF(dataF);
	Py_XDECREF(validF);
	
	return signalsF;
	
fail:
	Py_XDECREF(data);
	Py_XDECREF(freq);
	Py_XDECREF(delay);
	Py_XDECREF(dataF);
	Py_XDECREF(validF);
	
	return NULL;
}

PyDoc_STRVAR(FEngineR2_doc, \
"Perform a series of overlapped Fourier transforms on real-valued data using\n\
OpenMP.\n\
\n\
Input arguments are:\n\
 * signals: 2-D numpy.int16 (stands by samples) array of data to FFT\n\
 * frequency: 1-D numpy.double array of frequency values in Hz for the\n\
   FFT channels\n\
 * delays: 1-D numpy.double array of delays to apply to each stand\n\
\n\
Input keywords are:\n\
 * LFFT: number of FFT channels to make (default=64)\n\
 * Overlap: number of overlapped FFTs to use (default=1)\n\
 * SampleRate: sample rate of the data (default=196e6)\n\
 * ClipLevel: count value of 'bad' data.  FFT windows with instantaneous powers\n\
   greater than or equal to this value greater are zeroed.  Setting the ClipLevel\n\
   to zero disables time-domain blanking\n\
\n\
Outputs:\n\
 * fsignals: 3-D numpy.complex64 (stands by channels by FFT_set) of FFTd\n\
   data\n\
 * valid: 2-D numpy.uint8 (stands by FFT_set) of whether or not the FFT\n\
   set is valid (1) or not (0)\n\
");


static PyObject *FEngineR3(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *freqs, *delays, *window, *signalsF;
	PyArrayObject *data=NULL, *freq=NULL, *delay=NULL, *dataF=NULL, *validF=NULL, *windowData=NULL;
	int nChan = 64;
	int Overlap = 1;
	int Clip = 0;
	double SampleRate = 196.0e6;

	long ij, i, j, k, nStand, nSamps, nFFT;
	
	static char *kwlist[] = {"signals", "freqs", "delays", "LFFT", "Overlap", "SampleRate", "ClipLevel", "window", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OOO|iidiO:set_callback", kwlist, &signals, &freqs, &delays, &nChan, &Overlap, &SampleRate, &Clip, &window)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		goto fail;
	} else {
		if(!PyCallable_Check(window)) {
			PyErr_Format(PyExc_TypeError, "window must be a callable function");
			goto fail;
		}
		Py_XINCREF(window);
		Py_XDECREF(windowFunc);
		windowFunc = window;
	}

	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_INT16, 2, 2);
	freq = (PyArrayObject *) PyArray_ContiguousFromObject(freqs, NPY_DOUBLE, 1, 1);
	delay = (PyArrayObject *) PyArray_ContiguousFromObject(delays, NPY_DOUBLE, 2, 2);
	
	// Calculate the windowing function
	window = Py_BuildValue("(i)", 2*nChan);
	window = PyObject_CallObject(windowFunc, window);
	windowData = (PyArrayObject *) PyArray_ContiguousFromObject(window, NPY_DOUBLE, 1, 1);
	Py_DECREF(window);
	
	// Check data dimensions
	if(PyArray_DIM(data, 0) != PyArray_DIM(delay, 0)) {
		PyErr_Format(PyExc_TypeError, "signals and delays have different stand counts");
		goto fail;
	}
	
	if(nChan != PyArray_DIM(freq, 0)) {
		PyErr_Format(PyExc_RuntimeError, "freqs has a different channel count than nChan");
		goto fail;
	}
	
	if(PyArray_DIM(freq, 0) != PyArray_DIM(delay, 1)) {
		PyErr_Format(PyExc_TypeError, "freqs and delays have different channel counts");
		goto fail;
	}
	
	// Get the properties of the data
	nStand = (long) PyArray_DIM(data, 0);
	nSamps = (long) PyArray_DIM(data, 1);
	
	// Compute the integer sample offset and the fractional sample delay for each stand
	long *fifo, fifoMax;
	double *frac;
	fifo = (long *) malloc(nStand*sizeof(long));
	frac = (double *) malloc(nStand*nChan*sizeof(double));
	if( fifo == NULL || frac == NULL ) {
		PyErr_Format(PyExc_MemoryError, "Cannot create fifo/fractional delay arrays");
		goto fail;
	}
	fifoMax = computeDelayComponents(delay, SampleRate, fifo, frac);
	
	// Find out how large the output array needs to be and initialize it
	nFFT = (nSamps - fifoMax) / ((2*nChan)/Overlap) - (2*nChan)/((2*nChan)/Overlap) + 1;
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
	
	Py_BEGIN_ALLOW_THREADS
	
	// Create the FFTW plan                          
	float complex *inP, *in;                          
	inP = (float complex *) fftwf_malloc(sizeof(float complex) * 2*nChan);
	fftwf_plan p;
	p = fftwf_plan_dft_1d(2*nChan, inP, inP, FFTW_FORWARD, FFTW_ESTIMATE);
	
	// Data indexing and access
	long secStart;
	short int *a;
	float complex *b;
	double *c, *e;
	unsigned char *d;
	a = (short int *) PyArray_DATA(data);
	b = (float complex *) PyArray_DATA(dataF);
	c = (double *) PyArray_DATA(freq);
	d = (unsigned char *) PyArray_DATA(validF);
	e = (double *) PyArray_DATA(windowData);
	
	// Time-domain blanking control
	double cleanFactor;
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(in, i, j, k, secStart, cleanFactor)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(dynamic)
		#endif
		for(ij=0; ij<nStand*nFFT; ij++) {
			i = ij / nFFT;
			j = ij % nFFT;
			
			in = (float complex *) fftwf_malloc(sizeof(float complex) * 2*nChan);
			
			cleanFactor = 1.0;
			secStart = *(fifo + i) + nSamps*i + 2*nChan*j/Overlap;
			
			for(k=0; k<2*nChan; k++) {
				in[k] = (float complex) *(a + secStart + k);
				
				if( Clip && cabsf(in[k]) >= Clip ) {
					cleanFactor = 0.0;
				}
				
				in[k] *= *(e + k);
			}
			
			fftwf_execute_dft(p, in, in);
			
			for(k=0; k<nChan; k++) {
				*(b + nChan*nFFT*i + nFFT*k + j)  = cleanFactor*in[k];
				*(b + nChan*nFFT*i + nFFT*k + j) *= cexp(TPI * *(c + k) * *(frac + nChan*i + k));
				*(b + nChan*nFFT*i + nFFT*k + j) /= sqrt(2*nChan);
			}
			
			*(d + nFFT*i + j) = (unsigned char) cleanFactor;
			
			fftwf_free(in);
		}
	}
	fftwf_destroy_plan(p);
	fftwf_free(inP);
	free(frac);
	free(fifo);
	
	Py_END_ALLOW_THREADS
	
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

PyDoc_STRVAR(FEngineR3_doc, \
"Perform a series of overlapped Fourier transforms on real-valued data using\n\
OpenMP and windows.\n\
\n\
Input arguments are:\n\
 * signals: 2-D numpy.int16 (stands by samples) array of data to FFT\n\
 * frequency: 1-D numpy.double array of frequency values in Hz for the\n\
   FFT channels\n\
 * delays: 1-D numpy.double array of delays to apply to each stand\n\
\n\
Input keywords are:\n\
 * LFFT: number of FFT channels to make (default=64)\n\
 * Overlap: number of overlapped FFTs to use (default=1)\n\
 * SampleRate: sample rate of the data (default=196e6)\n\
 * window: Callable Python function for generating the window\n\
 * ClipLevel: count value of 'bad' data.  FFT windows with instantaneous powers\n\
   greater than or equal to this value greater are zeroed.  Setting the ClipLevel\n\
   to zero disables time-domain blanking\n\
\n\
Outputs:\n\
 * fsignals: 3-D numpy.complex64 (stands by channels by FFT_set) of FFTd\n\
   data\n\
 * valid: 2-D numpy.uint8 (stands by FFT_set) of whether or not the FFT\n\
   set is valid (1) or not (0)\n\
");


static PyObject *FEngineC2(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *freqs, *delays, *signalsF;
	PyArrayObject *data=NULL, *freq=NULL, *delay=NULL, *dataF=NULL, *validF=NULL;
	int nChan = 64;
	int Overlap = 1;
	int Clip = 0;
	double SampleRate = 1.0e5;
	
	long ij, i, j, k, nStand, nSamps, nFFT;
	
	static char *kwlist[] = {"signals", "freqs", "delays", "LFFT", "Overlap", "SampleRate", "ClipLevel", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OOO|iidi", kwlist, &signals, &freqs, &delays, &nChan, &Overlap, &SampleRate, &Clip)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		goto fail;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_COMPLEX64, 2, 2);
	freq = (PyArrayObject *) PyArray_ContiguousFromObject(freqs, NPY_DOUBLE, 1, 1);
	delay = (PyArrayObject *) PyArray_ContiguousFromObject(delays, NPY_DOUBLE, 2, 2);
	
	// Check data dimensions
	if(PyArray_DIM(data, 0) != PyArray_DIM(delay, 0)) {
		PyErr_Format(PyExc_TypeError, "signals and delays have different stand counts");
		goto fail;
	}
	
	if(nChan != PyArray_DIM(freq, 0)) {
		PyErr_Format(PyExc_RuntimeError, "freqs has a different channel count than nChan");
		goto fail;
	}
	
	if(PyArray_DIM(freq, 0) != PyArray_DIM(delay, 1)) {
		PyErr_Format(PyExc_TypeError, "freqs and delays have different channel counts");
		goto fail;
	}

	// Get the properties of the data
	nStand = (long) PyArray_DIM(data, 0);
	nSamps = (long) PyArray_DIM(data, 1);
	
	// Compute the integer sample offset and the fractional sample delay for each stand
	long *fifo, fifoMax;
	double *frac;
	fifo = (long *) malloc(nStand*sizeof(long));
	frac = (double *) malloc(nStand*nChan*sizeof(double));
	if( fifo == NULL || frac == NULL ) {
		PyErr_Format(PyExc_MemoryError, "Cannot create fifo/fractional delay arrays");
		goto fail;
	}
	fifoMax = computeDelayComponents(delay, SampleRate, fifo, frac);

	// Find out how large the output array needs to be and initialize it
	nFFT = (nSamps - fifoMax) / (nChan/Overlap) - nChan/(nChan/Overlap) + 1;
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
	
	Py_BEGIN_ALLOW_THREADS
	
	// Create the FFTW plan
	float complex *inP, *in;
	inP = (float complex*) fftwf_malloc(sizeof(float complex) * nChan);
	fftwf_plan p;
	p = fftwf_plan_dft_1d(nChan, inP, inP, FFTW_FORWARD, FFTW_ESTIMATE);
	
	// Data indexing and access
	long secStart;
	float complex *a, *b;
	double *c;
	unsigned char *d;
	a = (float complex *) PyArray_DATA(data);
	b = (float complex *) PyArray_DATA(dataF);
	c = (double *) PyArray_DATA(freq);
	d = (unsigned char *) PyArray_DATA(validF);
	
	// Time-domain blanking control
	double cleanFactor;
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(in, i, j, k, secStart, cleanFactor)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(dynamic)
		#endif
		for(ij=0; ij<nStand*nFFT; ij++) {
			i = ij / nFFT;
			j = ij % nFFT;
			
			in = (float complex*) fftwf_malloc(sizeof(float complex) * nChan);
			
			cleanFactor = 1.0;
			secStart = *(fifo + i) + nSamps*i + nChan*j/Overlap;
			
			for(k=0; k<nChan; k++) {
				in[k] = *(a + secStart + k);
				
				if( Clip && cabsf(in[k]) >= Clip ) {
					cleanFactor = 0.0;
				}
			}
			
			fftwf_execute_dft(p, in, in);
			
			for(k=0; k<nChan/2; k++) {
				*(b + nChan*nFFT*i + nFFT*k + j)  = cleanFactor*in[k+nChan/2+nChan%2];
				*(b + nChan*nFFT*i + nFFT*k + j) *= cexp(TPI * *(c + k) * *(frac + nChan*i + k));
				*(b + nChan*nFFT*i + nFFT*k + j) *= cexp(TPI * *(c + nChan/2) / SampleRate * *(fifo + i));
				*(b + nChan*nFFT*i + nFFT*k + j) /= sqrt(nChan);
			}
			for(k=nChan/2; k<nChan; k++) {
				*(b + nChan*nFFT*i + nFFT*k + j)  = cleanFactor*in[k-nChan/2];
				*(b + nChan*nFFT*i + nFFT*k + j) *= cexp(TPI * *(c + k) * *(frac + nChan*i + k));
				*(b + nChan*nFFT*i + nFFT*k + j) *= cexp(TPI * *(c + nChan/2) / SampleRate * *(fifo + i));
				*(b + nChan*nFFT*i + nFFT*k + j) /= sqrt(nChan);
			}
			
			*(d + nFFT*i + j) = (unsigned char) cleanFactor;
			
			fftwf_free(in);
		}
	}
	fftwf_destroy_plan(p);
	fftwf_free(inP);
	free(frac);
	free(fifo);
	
	Py_END_ALLOW_THREADS
	
	signalsF = Py_BuildValue("(OO)", PyArray_Return(dataF), PyArray_Return(validF));
	
	Py_XDECREF(data);
	Py_XDECREF(freq);
	Py_XDECREF(delay);
	Py_XDECREF(dataF);
	Py_XDECREF(validF);
	
	return signalsF;

fail:
	Py_XDECREF(data);
	Py_XDECREF(freq);
	Py_XDECREF(delay);
	Py_XDECREF(dataF);
	Py_XDECREF(validF);
	
	return NULL;
}

PyDoc_STRVAR(FEngineC2_doc, \
"Perform a series of overlapped Fourier transforms on complex-valued data\n\
using OpenMP.\n\
\n\
Input arguments are:\n\
 * signals: 2-D numpy.complex64 (stands by samples) array of data to FFT\n\
 * frequency: 1-D numpy.double array of frequency values in Hz for the\n\
   FFT channels\n\
 * delays: 1-D numpy.double array of delays to apply to each stand\n\
\n\
Input keywords are:\n\
 * LFFT: number of FFT channels to make (default=64)\n\
 * Overlap: number of overlapped FFTs to use (default=1)\n\
 * SampleRate: sample rate of the data (default=100e3)\n\
 * ClipLevel: count value of 'bad' data.  FFT windows with instantaneous powers\n\
   greater than or equal to this value greater are zeroed.  Setting the ClipLevel\n\
   to zero disables time-domain blanking\n\
\n\
Outputs:\n\
 * fsignals: 3-D numpy.complex64 (stands by channels by FFT_set) of FFTd\n\
   data\n\
 * valid: 2-D numpy.uint8 (stands by FFT_set) of whether or not the FFT\n\
   set is valid (1) or not (0)\n\
");


static PyObject *FEngineC3(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *freqs, *delays, *window, *signalsF;
	PyArrayObject *data=NULL, *freq=NULL, *delay=NULL, *dataF=NULL, *validF=NULL, *windowData=NULL;
	int nChan = 64;
	int Overlap = 1;
	int Clip = 0;
	double SampleRate = 1.0e5;

	long ij, i, j, k, nStand, nSamps, nFFT;

	static char *kwlist[] = {"signals", "freqs", "delays", "LFFT", "Overlap", "SampleRate", "ClipLevel", "window", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OOO|iidiO:set_callback", kwlist, &signals, &freqs, &delays, &nChan, &Overlap, &SampleRate, &Clip, &window)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		goto fail;
	} else {
		if(!PyCallable_Check(window)) {
			PyErr_Format(PyExc_TypeError, "window must be a callable function");
			goto fail;
		}
		Py_XINCREF(window);
		Py_XDECREF(windowFunc);
		windowFunc = window;
	}

	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_COMPLEX64, 2, 2);
	freq = (PyArrayObject *) PyArray_ContiguousFromObject(freqs, NPY_DOUBLE, 1, 1);
	delay = (PyArrayObject *) PyArray_ContiguousFromObject(delays, NPY_DOUBLE, 2, 2);
	
	// Calculate the windowing function
	window = Py_BuildValue("(i)", nChan);
	window = PyObject_CallObject(windowFunc, window);
	windowData = (PyArrayObject *) PyArray_ContiguousFromObject(window, NPY_DOUBLE, 1, 1);
	Py_DECREF(window);
	
	// Check data dimensions
	if(PyArray_DIM(data, 0) != PyArray_DIM(delay, 0)) {
		PyErr_Format(PyExc_TypeError, "signals and delays have different stand counts");
		goto fail;
	}
	
	if(nChan != PyArray_DIM(freq, 0)) {
		PyErr_Format(PyExc_RuntimeError, "freqs has a different channel count than nChan");
		goto fail;
	}
	
	if(PyArray_DIM(freq, 0) != PyArray_DIM(delay, 1)) {
		PyErr_Format(PyExc_TypeError, "freqs and delays have different channel counts");
		goto fail;
	}

	// Get the properties of the data
	nStand = (long) PyArray_DIM(data, 0);
	nSamps = (long) PyArray_DIM(data, 1);
	
	// Compute the integer sample offset and the fractional sample delay for each stand
	long *fifo, fifoMax;
	double *frac;
	fifo = (long *) malloc(nStand*sizeof(long));
	frac = (double *) malloc(nStand*nChan*sizeof(double));
	if( fifo == NULL || frac == NULL ) {
		PyErr_Format(PyExc_MemoryError, "Cannot create fifo/fractional delay arrays");
		goto fail;
	}
	fifoMax = computeDelayComponents(delay, SampleRate, fifo, frac);

	// Find out how large the output array needs to be and initialize it
	nFFT = (nSamps - fifoMax) / (nChan/Overlap) - nChan/(nChan/Overlap) + 1;
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
	
	Py_BEGIN_ALLOW_THREADS
	
	// Create the FFTW plan
	float complex *inP, *in;
	inP = (float complex*) fftwf_malloc(sizeof(float complex) * nChan);
	fftwf_plan p;
	p = fftwf_plan_dft_1d(nChan, inP, inP, FFTW_FORWARD, FFTW_ESTIMATE);
	
	// Data indexing and access
	long secStart;
	float complex *a, *b;
	double *c, *e;
	unsigned char *d;
	a = (float complex *) PyArray_DATA(data);
	b = (float complex *) PyArray_DATA(dataF);
	c = (double *) PyArray_DATA(freq);
	d = (unsigned char *) PyArray_DATA(validF);
	e = (double *) PyArray_DATA(windowData);
	
	// Time-domain blanking control
	double cleanFactor;
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(in, i, j, k, secStart, cleanFactor)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(dynamic)
		#endif
		for(ij=0; ij<nStand*nFFT; ij++) {
			i = ij / nFFT;
			j = ij % nFFT;
			
			in = (float complex*) fftwf_malloc(sizeof(float complex) * nChan);
			
			cleanFactor = 1.0;
			secStart = *(fifo + i) + nSamps*i + nChan*j/Overlap;
			
			for(k=0; k<nChan; k++) {
				in[k] = *(a + secStart + k);
				
				if( Clip && cabsf(in[k]) >= Clip ) {
					cleanFactor = 0.0;
				}
				
				in[k] *= *(e + k);
			}
			
			fftwf_execute_dft(p, in, in);
			
			for(k=0; k<nChan/2; k++) {
				*(b + nChan*nFFT*i + nFFT*k + j)  = cleanFactor*in[k+nChan/2+nChan%2];
				*(b + nChan*nFFT*i + nFFT*k + j) *= cexp(TPI * *(c + k) * *(frac + nChan*i + k));
				*(b + nChan*nFFT*i + nFFT*k + j) *= cexp(TPI * *(c + nChan/2) / SampleRate * *(fifo + i));
				*(b + nChan*nFFT*i + nFFT*k + j) /= sqrt(nChan);
			}
			for(k=nChan/2; k<nChan; k++) {
				*(b + nChan*nFFT*i + nFFT*k + j)  = cleanFactor*in[k-nChan/2];
				*(b + nChan*nFFT*i + nFFT*k + j) *= cexp(TPI * *(c + k) * *(frac + nChan*i + k));
				*(b + nChan*nFFT*i + nFFT*k + j) *= cexp(TPI * *(c + nChan/2) / SampleRate * *(fifo + i));
				*(b + nChan*nFFT*i + nFFT*k + j) /= sqrt(nChan);
			}
			
			*(d + nFFT*i + j) = (unsigned char) cleanFactor;
			
			fftwf_free(in);
		}
	}
	fftwf_destroy_plan(p);
	fftwf_free(inP);
	free(frac);
	free(fifo);
	
	Py_END_ALLOW_THREADS
	
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

PyDoc_STRVAR(FEngineC3_doc, \
"Perform a series of overlapped Fourier transforms on complex-valued data\n\
using OpenMP and allow for windowing of the data.\n\
\n\
Input arguments are:\n\
 * signals: 2-D numpy.complex64 (stands by samples) array of data to FFT\n\
 * frequency: 1-D numpy.double array of frequency values in Hz for the\n\
   FFT channels\n\
 * delays: 1-D numpy.double array of delays to apply to each stand\n\
\n\
Input keywords are:\n\
 * LFFT: number of FFT channels to make (default=64)\n\
 * Overlap: number of overlapped FFTs to use (default=1)\n\
 * SampleRate: sample rate of the data (default=100e3)\n\
 * window: Callable Python function for generating the window\n\
 * ClipLevel: count value of 'bad' data.  FFT windows with instantaneous powers\n\
   greater than or equal to this value greater are zeroed.  Setting the ClipLevel\n\
   to zero disables time-domain blanking\n\
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

static PyObject *XEngine2(PyObject *self, PyObject *args) {
	PyObject *signals1, *signals2, *sigValid1, *sigValid2, *output;
	PyArrayObject *data1=NULL, *data2=NULL, *valid1=NULL, *valid2=NULL, *vis=NULL;
	long nStand, nChan, nFFT, nBL;	

	if(!PyArg_ParseTuple(args, "OOOO", &signals1, &signals2, &sigValid1, &sigValid2)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		goto fail;
	}

	// Bring the data into C and make it usable
	data1 = (PyArrayObject *) PyArray_ContiguousFromObject(signals1, NPY_COMPLEX64, 3, 3);
	data2 = (PyArrayObject *) PyArray_ContiguousFromObject(signals2, NPY_COMPLEX64, 3, 3);
	valid1 = (PyArrayObject *) PyArray_ContiguousFromObject(sigValid1, NPY_UINT8, 2, 2);
	valid2 = (PyArrayObject *) PyArray_ContiguousFromObject(sigValid2, NPY_UINT8, 2, 2);

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
	float complex tempVis;
	float complex *a, *b, *v;
	a = (float complex *) PyArray_DATA(data1);
	b = (float complex *) PyArray_DATA(data2);
	v = (float complex *) PyArray_DATA(vis);
	
	// Time-domain blanking control
	long nActVis;
	unsigned char *u1, *u2;
	u1 = (unsigned char *) PyArray_DATA(valid1);
	u2 = (unsigned char *) PyArray_DATA(valid2);
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(c, f, nActVis, tempVis)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(dynamic)
		#endif
		for(bl=0; bl<nBL; bl++) {
			nActVis = 0;
			for(f=0; f<nFFT; f++) {
				nActVis += (long) (*(u1 + mapper[bl][0]*nFFT + f) & *(u2 + mapper[bl][1]*nFFT + f));
			}
			
			for(c=0; c<nChan; c++) {
				cblas_cdotc_sub(nFFT, (b + mapper[bl][1]*nChan*nFFT + c*nFFT), 1, (a + mapper[bl][0]*nChan*nFFT + c*nFFT), 1, &tempVis);
				*(v + bl*nChan + c) = tempVis / nActVis;
			}
		}
	}
	
	Py_END_ALLOW_THREADS
	
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
	dataX = (PyArrayObject *) PyArray_ContiguousFromObject(signalsX, NPY_COMPLEX64, 3, 3);
	dataY = (PyArrayObject *) PyArray_ContiguousFromObject(signalsY, NPY_COMPLEX64, 3, 3);
	validX = (PyArrayObject *) PyArray_ContiguousFromObject(sigValidX, NPY_UINT8, 2, 2);
	validY = (PyArrayObject *) PyArray_ContiguousFromObject(sigValidY, NPY_UINT8, 2, 2);

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
	float complex tempVis;
	float complex *a, *b, *v;
	a = (float complex *) PyArray_DATA(dataX);
	b = (float complex *) PyArray_DATA(dataY);
	v = (float complex *) PyArray_DATA(vis);
	
	// Time-domain blanking control
	long nActVisPureX, nActVisPureY, nActVisCross;
	unsigned char *u1, *u2;
	u1 = (unsigned char *) PyArray_DATA(validX);
	u2 = (unsigned char *) PyArray_DATA(validY);
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(c, f, nActVisPureX, nActVisPureY, nActVisCross, tempVis)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(dynamic)
		#endif
		for(bl=0; bl<nBL; bl++) {
			nActVisPureX = 0;
			nActVisPureY = 0;
			nActVisCross = 0;
			for(f=0; f<nFFT; f++) {
				nActVisPureX += (long) (*(u1 + mapper[bl][0]*nFFT + f) & *(u1 + mapper[bl][1]*nFFT + f));
				nActVisPureY += (long) (*(u2 + mapper[bl][0]*nFFT + f) & *(u2 + mapper[bl][1]*nFFT + f));
				nActVisCross += (long) (*(u1 + mapper[bl][0]*nFFT + f) & *(u2 + mapper[bl][1]*nFFT + f));
			}
			
			for(c=0; c<nChan; c++) {
				// XX
				cblas_cdotc_sub(nFFT, (a + mapper[bl][1]*nChan*nFFT + c*nFFT), 1, (a + mapper[bl][0]*nChan*nFFT + c*nFFT), 1, &tempVis);
				*(v + 0*nBL*nChan + bl*nChan + c) = tempVis / nActVisPureX;
				
				// XY
				cblas_cdotc_sub(nFFT, (b + mapper[bl][1]*nChan*nFFT + c*nFFT), 1, (a + mapper[bl][0]*nChan*nFFT + c*nFFT), 1, &tempVis);
				*(v + 1*nBL*nChan + bl*nChan + c) = tempVis / nActVisCross;
				
				// YX
				*(v + 2*nBL*nChan + bl*nChan + c) = conjf(*(v + 1*nBL*nChan + bl*nChan + c));
				
				// YY
				cblas_cdotc_sub(nFFT, (b + mapper[bl][1]*nChan*nFFT + c*nFFT), 1, (b + mapper[bl][0]*nChan*nFFT + c*nFFT), 1, &tempVis);
				*(v + 3*nBL*nChan + bl*nChan + c) = tempVis / nActVisPureY;
			}
		}
	}
	
	Py_END_ALLOW_THREADS
	
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
 * fsignals1: 3-D numpy.cdouble (stand by channels by FFT_set) array of FFTd\n\
   data from an F engine.\n\
 * fsignals2: 3-D numpy.cdouble (stand by channels by FFT_set) array of FFTd\n\
   data from an F engine.\n\
 * sigValid1: 1-D numpy.uint8 (FFT_set) array of whether or not the FFT_set is\n\
   valid (1) or not (0) for the first signal.\n\
 * sigValid2: 1-D numpy.uint8 (FFT_set) array of whether or not the FFT_set is\n\
   valid (1) or not (0) for the second signal.\n\
\n\
Ouputs:\n\
  * visibility: 3-D numpy.cdouble (Stokes parameter (XX,XY,YX,YY) by baseline by\n\
  channel) array of cross-correlated and averaged visibility data.\n\
");


/*
  Module Setup - Function Definitions and Documentation
*/

static PyMethodDef CorrelatorMethods[] = {
	{"FEngineR2", (PyCFunction) FEngineR2, METH_VARARGS|METH_KEYWORDS, FEngineR2_doc}, 
	{"FEngineR3", (PyCFunction) FEngineR3, METH_VARARGS|METH_KEYWORDS, FEngineR3_doc}, 
	{"FEngineC2", (PyCFunction) FEngineC2, METH_VARARGS|METH_KEYWORDS, FEngineC2_doc}, 
	{"FEngineC3", (PyCFunction) FEngineC3, METH_VARARGS|METH_KEYWORDS, FEngineC3_doc}, 
	{"XEngine2",  (PyCFunction) XEngine2,  METH_VARARGS,               XEngine2_doc }, 
	{"XEngine3",  (PyCFunction) XEngine3,  METH_VARARGS,               XEngine3_doc }, 
	{NULL,        NULL,      0,                          NULL         }
};

PyDoc_STRVAR(correlator_doc, \
"C-based F and X engines for the LWA software FX correlator.  These function\n\
are meant to provide an alternative to the lsl.correlator.fx.correlate function and \n\
provide a much-needed speed boost to cross-correlation.\n\
\n\
The function defined in this module are:\n\
  * FEngineR2 -F-engine for computing a series of overlapped Fourier transforms with\n\
    delay corrections for a real-valued (TBW) signal from a collection of stands all at\n\
    once.\n\
  * FEngineR3 - Similar to FEngineR2, but allows for a window function to be applied\n\
    to the data.\n\
  * FEngineC2 - F-engine for computing a series of overlapped Fourier transforms with\n\
    delay corrections for a complex-valued (TBN) signal from a collection of stands all at\n\
    once.\n\
  * FEngineC3 - Similar to FEngineC2, but allows for a window function to be applied\n\
    to the data.\n\
  * XEngine2 - Similar to XEngine, but works with a collection of stands all at\n\
    once.\n\
\n\
See the inidividual functions for more details.");


/*
  Module Setup - Initialization
*/

PyMODINIT_FUNC init_core(void) {
	char filename[256];
	PyObject *m, *pModule, *pDataPath;

	// Module definitions and functions
	m = Py_InitModule3("_core", CorrelatorMethods, correlator_doc);
	import_array();
	
	// Version and revision information
	PyModule_AddObject(m, "__version__", PyString_FromString("0.6"));
	PyModule_AddObject(m, "__revision__", PyString_FromString("$Rev$"));
	
	// LSL FFTW Wisdom
	pModule = PyImport_ImportModule("lsl.common.paths");
	if( pModule != NULL ) {
		pDataPath = PyObject_GetAttrString(pModule, "data");
		sprintf(filename, "%s/fftwf_wisdom.txt", PyString_AsString(pDataPath));
		read_wisdom(filename, m);
	} else {
		PyErr_Warn(PyExc_RuntimeWarning, "Cannot load the LSL FFTWF wisdom");
	}
}


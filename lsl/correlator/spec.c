#include "Python.h"
#include <math.h>
#include <stdio.h>
#ifdef _MKL
	#include "mkl_cblas.h"
	#include "fftw3.h"
#else
	#include <cblas.h>
	#include <fftw3.h>
#endif
#include <stdlib.h>
#include <complex.h>

#ifdef _OPENMP
	#include <omp.h>
	#ifdef _MKL
		#include "fftw3_mkl.h"
	#endif
#endif

#include "numpy/arrayobject.h"

#define PI 3.1415926535898
#define imaginary _Complex_I

/*
  Static declaration of the number of taps to use
*/

#define nTaps 64


/*
  Holder for window function callback
*/

static PyObject *windowFunc = NULL;


/*
  Sinc function for use by the polyphase filter bank
*/

double sinc(double x) {
	if(x == 0.0) {
		return 1.0;
	} else {
		return sin(x*PI)/(x*PI);
	}
}


/*
  FFT Functions
    1. FPSDR2 - FFT a real-valued collection of signals
    2. FPSDR3 - window the data and FFT a real-valued collection of signals
    3. FPSDC2 - FFT a complex-valued collection of signals
    4. FPSDC3 - window the data and FFT a complex-valued collection of signals
*/

static PyObject *FPSDR2(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *signalsF;
	PyArrayObject *data, *dataF;
	int nChan = 64;
	int Overlap = 1;
	int Clip = 0;

	long i, j, k, nStand, nSamps, nFFT;
	
	static char *kwlist[] = {"signals", "LFFT", "Overlap", "ClipLevel", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|iii", kwlist, &signals, &nChan, &Overlap, &Clip)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_INT16, 2, 2);
	
	// Get the properties of the data
	nStand = (long) data->dimensions[0];
	nSamps = (long) data->dimensions[1];
	
	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / 2 / nChan * Overlap - Overlap + 1;
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChan;
	dataF = (PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_DOUBLE);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		return NULL;
	}
	PyArray_FILLWBYTE(dataF, 0);
	
	// Create the FFTW plan                          
	fftw_complex *inP, *in;                          
	inP = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * 2*nChan);
	fftw_plan p;
	p = fftw_plan_dft_1d(2*nChan, inP, inP, FFTW_FORWARD, FFTW_MEASURE);
	
	// Integer delay, FFT, and fractional delay
	long secStart, fftIndex;
	short int *a;
	double *b;
	a = (short int *) data->data;
	b = (double *) dataF->data;
	
	// Time-domain blanking control
	double cleanFactor;
	long nActFFT[nStand];
	
	#ifdef _OPENMP
		#ifdef _MKL
			fftw3_mkl.number_of_user_threads = omp_get_num_threads();
		#endif
		#pragma omp parallel default(shared) private(in, secStart, i, j, k, fftIndex, cleanFactor)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			nActFFT[i] = 0;
			in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * 2*nChan);
			
			for(j=0; j<nFFT; j++) {
				cleanFactor = 1.0;
				secStart = (long) (nSamps * i) + 2*nChan*j/Overlap;
				
				for(k=0; k<2*nChan; k++) {
					in[k][0] = (double) *(a + secStart + k);
					in[k][1] = 0.0;
					
					if( Clip && (in[k][0] >= Clip || in[k][0] <= -Clip) ) {
						cleanFactor = 0.0;
					}
				}
				
				fftw_execute_dft(p, in, in);
				
				for(k=0; k<nChan; k++) {
					fftIndex = k;
					*(b + nChan*i + k) += cleanFactor*in[fftIndex][0]*in[fftIndex][0];
					*(b + nChan*i + k) += cleanFactor*in[fftIndex][1]*in[fftIndex][1];
				}
				
				nActFFT[i] += (long) cleanFactor;
			}
			fftw_free(in);
		}
	}
	fftw_destroy_plan(p);
	fftw_free(inP);
	
	// cblas_dscal(nChan*nStand, 1.0/(2*nChan*nFFT), b, 1);
	for(i=0; i<nStand; i++) {
		cblas_dscal(nChan, 1.0/(2*nChan*nActFFT[i]), (b + i*nChan), 1);
	}

	Py_XDECREF(data);

	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(FPSDR2_doc, \
"Perform a series of Fourier transforms on real-valued data to get the PSD.\n\
\n\
Input arguments are:\n\
 * signals: 2-D numpy.int16 (stands by samples) array of data to FFT\n\
\n\
Input keywords are:\n\
 * LFFT: number of FFT channels to make (default=64)\n\
 * Overlap: number of overlapped FFTs to use (default=1)\n\
 * ClipLevel: count value of 'bad' data.  FFT windows with instantaneous powers\n\
   greater than or equal to this value greater are zeroed.  Setting the ClipLevel\n\
   to zero disables time-domain blanking\n\
\n\
Outputs:\n\
 * psd: 2-D numpy.double (stands by channels) of PSD data\n\
");


static PyObject *FPSDR3(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *signalsF, *window;
	PyArrayObject *data, *dataF, *windowData;
	int nChan = 64;
	int Overlap = 1;
	int Clip = 0;
	
	long i, j, k, nStand, nSamps, nFFT;
	
	static char *kwlist[] = {"signals", "LFFT", "Overlap", "ClipLevel", "window", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|iiiO:set_callback", kwlist, &signals, &nChan, &Overlap, &Clip, &window)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	} else {
		if(!PyCallable_Check(window)) {
			PyErr_Format(PyExc_TypeError, "window must be a callable function");
		}
		Py_XINCREF(window);
		Py_XDECREF(windowFunc);
		windowFunc = window;
	}

	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_INT16, 2, 2);
	
	// Calculate the windowing function
	window = Py_BuildValue("(i)", 2*nChan);
	window = PyObject_CallObject(windowFunc, window);
	windowData = (PyArrayObject *) PyArray_ContiguousFromObject(window, NPY_DOUBLE, 1, 1);
	Py_DECREF(window);
	
	// Get the properties of the data
	nStand = (long) data->dimensions[0];
	nSamps = (long) data->dimensions[1];

	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / 2 / nChan * Overlap - Overlap + 1;
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChan;
	dataF = (PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_DOUBLE);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		Py_XDECREF(windowData);
		return NULL;
	}
	PyArray_FILLWBYTE(dataF, 0);
	
	// Create the FFTW plan                          
	fftw_complex *inP, *in;                          
	inP = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * 2*nChan);
	fftw_plan p;
	p = fftw_plan_dft_1d(2*nChan, inP, inP, FFTW_FORWARD, FFTW_ESTIMATE);
	
	// Integer delay, FFT, and fractional delay
	long secStart, fftIndex;
	short int *a;
	double *b, *c;
	a = (short int *) data->data;
	b = (double *) dataF->data;
	c = (double *) windowData->data;
	
	// Time-domain blanking control
	double cleanFactor;
	long nActFFT[nStand];
	
	#ifdef _OPENMP
		#ifdef _MKL
			fftw3_mkl.number_of_user_threads = omp_get_num_threads();
		#endif
		#pragma omp parallel default(shared) private(in, secStart, i, j, k, fftIndex, cleanFactor)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			nActFFT[i] = 0;
			in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * 2*nChan);
			
			for(j=0; j<nFFT; j++) {
				cleanFactor = 1.0;
				secStart = nSamps * i + 2*nChan*j/Overlap;
				
				for(k=0; k<2*nChan; k++) {
					in[k][0] = (double) *(a + secStart + k) * *(c + k);
					in[k][1] = 0.0;
					
					if( Clip && (*(a + secStart + k) >= Clip || *(a + secStart + k) <= -Clip) ) {
						cleanFactor = 0.0;
					}
				}
				
				fftw_execute_dft(p, in, in);
				
				for(k=0; k<nChan; k++) {
					fftIndex = k;
					*(b + nChan*i + k) += cleanFactor*in[fftIndex][0]*in[fftIndex][0];
					*(b + nChan*i + k) += cleanFactor*in[fftIndex][1]*in[fftIndex][1];
				}
				
				nActFFT[i] += (long) cleanFactor;
			}
			fftw_free(in);
		}
	}
	fftw_destroy_plan(p);
	fftw_free(inP);
	
	// cblas_dscal(nChan*nStand, 1.0/(2*nChan*nFFT), b, 1);
	for(i=0; i<nStand; i++) {
		cblas_dscal(nChan, 1.0/(2*nChan*nActFFT[i]), (b + i*nChan), 1);
	}

	Py_XDECREF(data);
	Py_XDECREF(windowData);

	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(FPSDR3_doc, \
"Perform a series of Fourier transforms with windows on real-valued data to\n\
get the PSD.\n\
\n\
Input arguments are:\n\
 * signals: 2-D numpy.int16 (stands by samples) array of data to FFT\n\
\n\
Input keywords are:\n\
 * LFFT: number of FFT channels to make (default=64)\n\
 * Overlap: number of overlapped FFTs to use (default=1)\n\
 * window: Callable Python function for generating the window\n\
 * ClipLevel: count value of 'bad' data.  FFT windows with instantaneous powers\n\
   greater than or equal to this value greater are zeroed.  Setting the ClipLevel\n\
   to zero disables time-domain blanking\n\
\n\
Outputs:\n\
 * psd: 2-D numpy.double (stands by channels) of PSD data\n\
");


static PyObject *FPSDC2(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *signalsF;
	PyArrayObject *data, *dataF;
	int nChan = 64;
	int Overlap = 1;
	int Clip = 0;

	long i, j, k, nStand, nSamps, nFFT;

	static char *kwlist[] = {"signals", "LFFT", "Overlap", "ClipLevel", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|iii", kwlist, &signals, &nChan, &Overlap, &Clip)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_COMPLEX64, 2, 2);

	// Get the properties of the data
	nStand = (long) data->dimensions[0];
	nSamps = (long) data->dimensions[1];
	
	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / nChan * Overlap - Overlap + 1;
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChan;
	dataF = (PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_DOUBLE);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		return NULL;
	}
	PyArray_FILLWBYTE(dataF, 0);

	// Create the FFTW plan
	fftw_complex *inP, *in;
	inP = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nChan);
	fftw_plan p;
	p = fftw_plan_dft_1d(nChan, inP, inP, FFTW_FORWARD, FFTW_MEASURE);

	// Integer delay, FFT, and fractional delay
	long secStart;
	float complex *a;
	double *b;
	a = (float complex *) data->data;
	b = (double *) dataF->data;
	
	// Time-domain blanking control
	double cleanFactor;
	long nActFFT[nStand];
	
	#ifdef _OPENMP
		#ifdef _MKL
			fftw3_mkl.number_of_user_threads = omp_get_num_threads();
		#endif
		#pragma omp parallel default(shared) private(in, secStart, i, j, k, cleanFactor)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			nActFFT[i] = 0;
			in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nChan);
			
			for(j=0; j<nFFT; j++) {
				cleanFactor = 1.0;
				secStart = nSamps * i + nChan*j/Overlap;
				
				for(k=0; k<nChan; k++) {
					in[k][0] = creal(*(a + secStart + k));
					in[k][1] = cimag(*(a + secStart + k));
					
					if( Clip && cabs(*(a + secStart + k)) >= Clip ) {
						cleanFactor = 0.0;
					}
				}
				
				fftw_execute_dft(p, in, in);
				
				for(k=0; k<nChan; k++) {
					*(b + nChan*i + k) += cleanFactor*in[k][0]*in[k][0];
					*(b + nChan*i + k) += cleanFactor*in[k][1]*in[k][1];
				}
				
				nActFFT[i] += (long) cleanFactor;
			}
			fftw_free(in);
		}
	}
	fftw_destroy_plan(p);
	fftw_free(inP);

	// Shift and scale FFTs
	double *temp, *temp2;
	temp2 = (double *) malloc(sizeof(double)*nChan/2);
	for(i=0; i<nStand; i++) {
		temp = b + nChan*i;
		memcpy(temp2, temp, sizeof(double)*nChan/2);
		memmove(temp, temp+nChan/2, sizeof(double)*nChan/2);
		memcpy(temp+nChan/2, temp2, sizeof(double)*nChan/2);
		
		cblas_dscal(nChan, 1.0/(nChan*nActFFT[i]), (b + i*nChan), 1);
	}
	free(temp2);

	Py_XDECREF(data);
	
	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(FPSDC2_doc, \
"Perform a series of Fourier transforms on complex-valued data to get the\n\
PSD.\n\
\n\
Input arguments are:\n\
 * signals: 2-D numpy.complex64 (stands by samples) array of data to FFT\n\
\n\
Input keywords are:\n\
 * LFFT: number of FFT channels to make (default=64)\n\
 * Overlap: number of overlapped FFTs to use (default=1)\n\
 * ClipLevel: count value of 'bad' data.  FFT windows with instantaneous powers\n\
   greater than or equal to this value greater are zeroed.  Setting the ClipLevel\n\
   to zero disables time-domain blanking\n\
\n\
Outputs:\n\
 * psd: 2-D numpy.double (stands by channels) of PSD data\n\
");


static PyObject *FPSDC3(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *signalsF, *window;
	PyArrayObject *data, *dataF, *windowData;
	int nChan = 64;
	int Overlap = 1;
	int Clip = 0;

	long i, j, k, nStand, nSamps, nFFT;

	static char *kwlist[] = {"signals", "LFFT", "Overlap", "ClipLevel", "window", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|iiiO:set_callback", kwlist, &signals, &nChan, &Overlap, &Clip, &window)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	} else {
		if(!PyCallable_Check(window)) {
			PyErr_Format(PyExc_TypeError, "window must be a callable function");
		}
		Py_XINCREF(window);
		Py_XDECREF(windowFunc);
		windowFunc = window;
	}

	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_COMPLEX64, 2, 2);
	
	// Calculate the windowing function
	window = Py_BuildValue("(i)", nChan);
	window = PyObject_CallObject(windowFunc, window);
	windowData = (PyArrayObject *) PyArray_ContiguousFromObject(window, NPY_DOUBLE, 1, 1);
	Py_DECREF(window);

	// Get the properties of the data
	nStand = (long) data->dimensions[0];
	nSamps = (long) data->dimensions[1];

	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / nChan * Overlap - Overlap + 1;
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChan;
	dataF = (PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_DOUBLE);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		Py_XDECREF(windowData);
		return NULL;
	}
	PyArray_FILLWBYTE(dataF, 0);

	// Create the FFTW plan
	fftw_complex *inP, *in;
	inP = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nChan);
	fftw_plan p;
	p = fftw_plan_dft_1d(nChan, inP, inP, FFTW_FORWARD, FFTW_MEASURE);

	// Integer delay, FFT, and fractional delay
	long secStart;
	float complex *a;
	double *b, *c;
	a = (float complex *) data->data;
	b = (double *) dataF->data;
	c = (double *) windowData->data;
	
	// Time-domain blanking control
	double cleanFactor;
	long nActFFT[nStand];
	
	#ifdef _OPENMP
		#ifdef _MKL
			fftw3_mkl.number_of_user_threads = omp_get_num_threads();
		#endif
		#pragma omp parallel default(shared) private(in, secStart, i, j, k, cleanFactor)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			nActFFT[i] = 0;
			in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nChan);
			
			for(j=0; j<nFFT; j++) {
				cleanFactor = 1.0;
				secStart = nSamps * i + nChan*j/Overlap;
				
				for(k=0; k<nChan; k++) {
					in[k][0] = creal(*(a + secStart + k)) * *(c + k);
					in[k][1] = cimag(*(a + secStart + k)) * *(c + k);
					
					if( Clip && cabs(*(a + secStart + k)) >= Clip ) {
						cleanFactor = 0.0;
					}
				}
				
				fftw_execute_dft(p, in, in);
				
				for(k=0; k<nChan; k++) {
					*(b + nChan*i + k) += cleanFactor*in[k][0]*in[k][0];
					*(b + nChan*i + k) += cleanFactor*in[k][1]*in[k][1];
				}
				
				nActFFT[i] += (long) cleanFactor;
			}
			fftw_free(in);
		}
	}
	fftw_destroy_plan(p);
	fftw_free(inP);
	
	// Shift and scale FFTs
	double *temp, *temp2;
	temp2 = (double *) malloc(sizeof(double)*nChan/2);
	for(i=0; i<nStand; i++) {
		temp = b + nChan*i;
		memcpy(temp2, temp, sizeof(double)*nChan/2);
		memmove(temp, temp+nChan/2, sizeof(double)*nChan/2);
		memcpy(temp+nChan/2, temp2, sizeof(double)*nChan/2);
		
		cblas_dscal(nChan, 1.0/(nChan*nActFFT[i]), (b + i*nChan), 1);
	}
	free(temp2);
	
	Py_XDECREF(data);
	Py_XDECREF(windowData);

	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(FPSDC3_doc, \
"Perform a series of Fourier transforms with windows on complex-valued data\n\
to get the PSD.\n\
\n\
Input arguments are:\n\
 * signals: 2-D numpy.complex64 (stands by samples) array of data to FFT\n\
\n\
Input keywords are:\n\
 * LFFT: number of FFT channels to make (default=64)\n\
 * Overlap: number of overlapped FFTs to use (default=1)\n\
 * window: Callable Python function for generating the window\n\
 * ClipLevel: count value of 'bad' data.  FFT windows with instantaneous powers\n\
   greater than or equal to this value greater are zeroed.  Setting the ClipLevel\n\
   to zero disables time-domain blanking\n\
\n\
Outputs:\n\
 * psd: 2-D numpy.double (stands by channels) of PSD data\n\
");


/*
  Polyphase Filterbank (PFB) Functions
    1. PPSDR2 - PFB a real-valued collection of signals
    2. PPSDR3 - window the data and PFB a real-valued collection of signals
    3. PPSDC2 - PFB a complex-valued collection of signals
    4. PPSDC3 - window the data and PFB a complex-valued collection of signals
*/

static PyObject *PPSDR2(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *signalsF;
	PyArrayObject *data, *dataF;
	int nChan = 64;
	int Overlap = 1;
	int Clip = 0;

	long i, j, k, m, nStand, nSamps, nFFT;
	
	static char *kwlist[] = {"signals", "LFFT", "Overlap", "ClipLevel", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|iii", kwlist, &signals, &nChan, &Overlap, &Clip)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_INT16, 2, 2);
	
	// Get the properties of the data
	nStand = (long) data->dimensions[0];
	nSamps = (long) data->dimensions[1];

	// Compute the filterbank window for the correct number of taps
	double fbWindow[2*nChan*nTaps];
	double complex tempFB[nChan-1];
	double tempB[nChan-1];
	for(i=0; i<2*nChan*nTaps; i++) {
		fbWindow[i] = sinc((double) (i - nTaps*nChan + 0.5)/2/nChan);
	}
	
	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / 2 / nChan * Overlap - Overlap + 1;
	nFFT = (nFFT / nTaps) * nTaps;
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChan;
	dataF = (PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_DOUBLE);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		return NULL;
	}
	PyArray_FILLWBYTE(dataF, 0);
	
	// Create the FFTW plan                          
	fftw_complex *inP, *in;                          
	inP = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * 2*nChan);
	fftw_plan p;
	p = fftw_plan_dft_1d(2*nChan, inP, inP, FFTW_FORWARD, FFTW_MEASURE);
	
	// Integer delay, FFT, and fractional delay
	long secStart, fftIndex;
	short int *a;
	double *b;
	a = (short int *) data->data;
	b = (double *) dataF->data;
	
	// Time-domain blanking control
	double cleanFactor;
	long nActFFT[nStand];
	
	#ifdef _OPENMP
		#ifdef _MKL
			fftw3_mkl.number_of_user_threads = omp_get_num_threads();
		#endif
		#pragma omp parallel default(shared) private(in, secStart, tempFB, tempB, i, j, k, m, fftIndex, cleanFactor)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			nActFFT[i] = 0;
			in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * 2*nChan);
			
			for(j=0; j<nFFT; j+=nTaps) {
				cleanFactor = 1.0;
				cblas_zdscal(nChan, 0.0, tempFB, 1);
				cblas_dscal(nChan, 0.0, tempB, 1);
				
				for(m=0; m<nTaps; m++) {
					secStart = nSamps * i + 2*nChan*(j+m)/Overlap;
					
					for(k=0; k<2*nChan; k++) {
						in[k][0] = (double) *(a + secStart + k) * fbWindow[2*nChan*m + k];
						in[k][1] = 0.0;
						
						if( Clip && (*(a + secStart + k) >= Clip || *(a + secStart + k) <= -Clip) ) {
							cleanFactor = 0.0;
						}
					}
				
					fftw_execute_dft(p, in, in);
				
					for(k=0; k<nChan; k++) {
						fftIndex = k;
						tempFB[k] += cleanFactor*in[fftIndex][0] + imaginary *  cleanFactor*in[fftIndex][1];
					}
				}
				
				#ifdef _MKL
					vzAbs(nChan, tempFB, tempB);
					vdSqr(nChan, tempB, tempB);
					vdAdd(nChan, (b+nChan*i), tempB, (b+nChan*i));
				#else
					for(k=0; k<nChan; k++) {
						*(b + nChan*i + k) += pow(cabs(tempFB[k]), 2);
					}
				#endif
				
				nActFFT[i] += (long) cleanFactor;
			}
			fftw_free(in);
		}
	}
	fftw_destroy_plan(p);
	fftw_free(inP);
	
	// cblas_dscal(nChan*nStand, ((float) nTaps)/(2*nChan*nFFT), b, 1);
	for(i=0; i<nStand; i++) {
		cblas_dscal(nChan, ((float) nTaps)/(2*nChan*nActFFT[i]), (b + i*nChan), 1);
	}

	Py_XDECREF(data);

	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(PPSDR2_doc, \
"Perform a series of filter bank transforms on real-valued data to get the\n\
PSD.\n\
\n\
Input arguments are:\n\
 * signals: 2-D numpy.int16 (stands by samples) array of data to FFT\n\
\n\
Input keywords are:\n\
 * LFFT: number of FFT channels to make (default=64)\n\
 * Overlap: number of overlapped FFTs to use (default=1)\n\
 * ClipLevel: count value of 'bad' data.  FFT windows with instantaneous powers\n\
   greater than or equal to this value greater are zeroed.  Setting the ClipLevel\n\
   to zero disables time-domain blanking\n\
\n\
Outputs:\n\
 * psd: 2-D numpy.double (stands by channels) of PSD data\n\
");


static PyObject *PPSDR3(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *signalsF, *window;
	PyArrayObject *data, *dataF, *windowData;
	int nChan = 64;
	int Overlap = 1;
	int Clip = 0;

	long i, j, k, m, nStand, nSamps, nFFT;
	
	static char *kwlist[] = {"signals", "LFFT", "Overlap", "ClipLevel", "window", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|iiiO:set_callback", kwlist, &signals, &nChan, &Overlap, &Clip, &window)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	} else {
		if(!PyCallable_Check(window)) {
			PyErr_Format(PyExc_TypeError, "window must be a callable function");
		}
		Py_XINCREF(window);
		Py_XDECREF(windowFunc);
		windowFunc = window;
	}

	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_INT16, 2, 2);
	
	// Calculate the windowing function
	window = Py_BuildValue("(i)", 2*nTaps*nChan);
	window = PyObject_CallObject(windowFunc, window);
	windowData = (PyArrayObject *) PyArray_ContiguousFromObject(window, NPY_DOUBLE, 1, 1);
	Py_DECREF(window);
	
	// Get the properties of the data
	nStand = (long) data->dimensions[0];
	nSamps = (long) data->dimensions[1];

	// Compute the filterbank window for the correct number of taps
	npy_intp *qLoc;
	double fbWindow[2*nChan*nTaps];
	double complex tempFB[nChan-1];
	double tempB[nChan-1];
	qLoc = PyDimMem_NEW(1);
	for(i=0; i<2*nChan*nTaps; i++) {
		qLoc[0] = (npy_intp) i;
		fbWindow[i] = sinc((double) (i - nTaps*nChan + 0.5)/2/nChan);
		fbWindow[i] *= *(double *) PyArray_GetPtr(windowData, qLoc);
	}
	PyDimMem_FREE(qLoc);

	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / 2 / nChan * Overlap - Overlap + 1;
	nFFT = (nFFT / nTaps) * nTaps;
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChan;
	dataF = (PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_DOUBLE);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		Py_XDECREF(windowData);
		return NULL;
	}
	PyArray_FILLWBYTE(dataF, 0);
	
	// Create the FFTW plan                          
	fftw_complex *inP, *in;                          
	inP = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * 2*nChan);
	fftw_plan p;
	p = fftw_plan_dft_1d(2*nChan, inP, inP, FFTW_FORWARD, FFTW_MEASURE);
	
	// Integer delay, FFT, and fractional delay
	long secStart, fftIndex;
	short int *a;
	double *b;
	a = (short int *) data->data;
	b = (double *) dataF->data;
	
	// Time-domain blanking control
	double cleanFactor;
	long nActFFT[nStand];

	#ifdef _OPENMP
		#ifdef _MKL
			fftw3_mkl.number_of_user_threads = omp_get_num_threads();
		#endif
		#pragma omp parallel default(shared) private(in, secStart, tempFB, tempB, i, j, k, m, fftIndex, cleanFactor)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			nActFFT[i] = 0;
			in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * 2*nChan);
			
			for(j=0; j<nFFT; j+=nTaps) {
				cleanFactor = 1.0;
				cblas_zdscal(nChan, 0.0, tempFB, 1);
				cblas_dscal(nChan, 0.0, tempB, 1);
				
				for(m=0; m<nTaps; m++) {
					secStart = nSamps * i + 2*nChan*(j+m)/Overlap;
					
					for(k=0; k<2*nChan; k++) {
						in[k][0] = (double) *(a + secStart + k) * fbWindow[2*nChan*m + k];
						in[k][1] = 0.0;
						
						if( Clip && (*(a + secStart + k) >= Clip || *(a + secStart + k) <= -Clip) ) {
							cleanFactor = 0.0;
						}
					}
				
					fftw_execute_dft(p, in, in);
				
					for(k=0; k<nChan; k++) {
						fftIndex = k;
						tempFB[k] += cleanFactor*in[fftIndex][0] + imaginary * cleanFactor*in[fftIndex][1];
					}
				}
				
				#ifdef _MKL
					vzAbs(nChan, tempFB, tempB);
					vdSqr(nChan, tempB, tempB);
					vdAdd(nChan, (b+nChan*i), tempB, (b+nChan*i));
				#else
					for(k=0; k<nChan; k++) {
						*(b + nChan*i + k) += pow(cabs(tempFB[k]), 2);
					}
				#endif
				
				nActFFT[i] += (long) cleanFactor;
			}
			fftw_free(in);
		}
	}
	fftw_destroy_plan(p);
	fftw_free(inP);
	
	// cblas_dscal(nChan*nStand, ((float) nTaps)/(2*nChan*nFFT), b, 1);
	for(i=0; i<nStand; i++) {
		cblas_dscal(nChan, ((float) nTaps)/(2*nChan*nActFFT[i]), (b + i*nChan), 1);
	}

	Py_XDECREF(data);
	Py_XDECREF(windowData);

	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(PPSDR3_doc, \
"Perform a series of filter bank transforms with windows on real-valued data\n\
to get the PSD.\n\
\n\
Input arguments are:\n\
 * signals: 2-D numpy.int16 (stands by samples) array of data to FFT\n\
\n\
Input keywords are:\n\
 * LFFT: number of FFT channels to make (default=64)\n\
 * Overlap: number of overlapped FFTs to use (default=1)\n\
 * window: Callable Python function for generating the window\n\
 * ClipLevel: count value of 'bad' data.  FFT windows with instantaneous powers\n\
   greater than or equal to this value greater are zeroed.  Setting the ClipLevel\n\
   to zero disables time-domain blanking\n\
\n\
Outputs:\n\
 * psd: 2-D numpy.double (stands by channels) of PSD data\n\
");


static PyObject *PPSDC2(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *signalsF;
	PyArrayObject *data, *dataF;
	int nChan = 64;
	int Overlap = 1;
	int Clip = 0;

	long i, j, k, m, nStand, nSamps, nFFT;

	static char *kwlist[] = {"signals", "LFFT", "Overlap", "ClipLevel", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|iii", kwlist, &signals, &nChan, &Overlap, &Clip)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_COMPLEX64, 2, 2);

	// Get the properties of the data
	nStand = (long) data->dimensions[0];
	nSamps = (long) data->dimensions[1];

	// Compute the filterbank window for the correct number of taps
	double fbWindow[nChan*nTaps];
	double complex tempFB[nChan-1];
	double tempB[nChan-1];
	for(i=0; i<nChan*nTaps; i++) {
		fbWindow[i] = sinc((double) (i - nTaps*nChan/2 + 0.5)/nChan);
	}

	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / nChan * Overlap - Overlap + 1;
	nFFT = (nFFT / nTaps) * nTaps;
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChan;
	dataF = (PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_DOUBLE);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		return NULL;
	}
	PyArray_FILLWBYTE(dataF, 0);

	// Create the FFTW plan
	fftw_complex *inP, *in;
	inP = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nChan);
	fftw_plan p;
	p = fftw_plan_dft_1d(nChan, inP, inP, FFTW_FORWARD, FFTW_ESTIMATE);

	// Integer delay, FFT, and fractional delay
	long secStart, fftIndex;
	float complex *a;
	double *b;
	a = (float complex *) data->data;
	b = (double *) dataF->data;
	
	// Time-domain blanking control
	double cleanFactor;
	long nActFFT[nStand];
	
	#ifdef _OPENMP
		#ifdef _MKL
			fftw3_mkl.number_of_user_threads = omp_get_num_threads();
		#endif
		#pragma omp parallel default(shared) private(in, secStart, tempFB, tempB, i, j, k, m, fftIndex, cleanFactor)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			nActFFT[i] = 0;
			in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nChan);
			
			for(j=0; j<nFFT; j+=nTaps) {
				cleanFactor = 1.0;
				cblas_zdscal(nChan, 0.0, tempFB, 1);
				cblas_dscal(nChan, 0.0, tempB, 1);

				for(m=0; m<nTaps; m++) {
					secStart = nSamps * i + nChan*(j+m)/Overlap;
					
					for(k=0; k<nChan; k++) {
						in[k][0] = creal(*(a + secStart + k)) * fbWindow[nChan*m + k];
						in[k][1] = cimag(*(a + secStart + k)) * fbWindow[nChan*m + k];
						
						if( Clip && cabs(*(a + secStart + k)) >= Clip ) {
							cleanFactor = 0.0;
						}
					}
				
					fftw_execute_dft(p, in, in);
				
					for(k=0; k<nChan; k++) {
						fftIndex = (k + nChan/2) % nChan;
						tempFB[k] += cleanFactor*in[fftIndex][0] + imaginary * cleanFactor*in[fftIndex][1];
					}
				}
				
				#ifdef _MKL
					vzAbs(nChan, tempFB, tempB);
					vdSqr(nChan, tempB, tempB);
					vdAdd(nChan, (b+nChan*i), tempB, (b+nChan*i));
				#else
					for(k=0; k<nChan; k++) {
						*(b + nChan*i + k) += pow(cabs(tempFB[k]), 2);
					}
				#endif
				
				nActFFT[i] += (long) cleanFactor;
			}
			fftw_free(in);
		}
	}
	fftw_destroy_plan(p);
	fftw_free(inP);

	// cblas_dscal(nChan*nStand, ((float) nTaps)/(nChan*nFFT), b, 1);
	for(i=0; i<nStand; i++) {
		cblas_dscal(nChan, ((float) nTaps)/(nChan*nActFFT[i]), (b + i*nChan), 1);
	}
	
	Py_XDECREF(data);

	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(PPSDC2_doc, \
"Perform a series of filter bank transforms on complex-valued data to get\n\
the PSD.\n\
\n\
Input arguments are:\n\
 * signals: 2-D numpy.complex64 (stands by samples) array of data to FFT\n\
\n\
Input keywords are:\n\
 * LFFT: number of FFT channels to make (default=64)\n\
 * Overlap: number of overlapped FFTs to use (default=1)\n\
 * ClipLevel: count value of 'bad' data.  FFT windows with instantaneous powers\n\
   greater than or equal to this value greater are zeroed.  Setting the ClipLevel\n\
   to zero disables time-domain blanking\n\
\n\
Outputs:\n\
 * psd: 2-D numpy.double (stands by channels) of PSD data\n\
");


static PyObject *PPSDC3(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *signalsF, *window;
	PyArrayObject *data, *dataF, *windowData;
	int nChan = 64;
	int Overlap = 1;
	int Clip = 0;

	long i, j, k, m, nStand, nSamps, nFFT;

	static char *kwlist[] = {"signals", "LFFT", "Overlap", "ClipLevel", "window", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|iiiO:set_callback", kwlist, &signals, &nChan, &Overlap, &Clip, &window)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	} else {
		if(!PyCallable_Check(window)) {
			PyErr_Format(PyExc_TypeError, "window must be a callable function");
		}
		Py_XINCREF(window);
		Py_XDECREF(windowFunc);
		windowFunc = window;
	}

	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_COMPLEX64, 2, 2);
	
	// Calculate the windowing function
	window = Py_BuildValue("(i)", nTaps*nChan);
	window = PyObject_CallObject(windowFunc, window);
	windowData = (PyArrayObject *) PyArray_ContiguousFromObject(window, NPY_DOUBLE, 1, 1);
	Py_DECREF(window);

	// Get the properties of the data
	nStand = (long) data->dimensions[0];
	nSamps = (long) data->dimensions[1];

	// Compute the filterbank window for the correct number of taps
	double fbWindow[nChan*nTaps];
	double complex tempFB[nChan-1];
	double tempB[nChan-1];
	double *c;
	c = (double *) windowData->data;
	for(i=0; i<nChan*nTaps; i++) {
		fbWindow[i] = sinc((double) (i - nTaps*nChan/2 + 0.5)/nChan);
		fbWindow[i] *= *(c + i);
	}

	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / nChan * Overlap - Overlap + 1;
	nFFT = (nFFT / nTaps) * nTaps;
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChan;
	dataF = (PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_DOUBLE);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		Py_XDECREF(windowData);
		return NULL;
	}
	PyArray_FILLWBYTE(dataF, 0);

	// Create the FFTW plan
	fftw_complex *inP, *in;
	inP = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nChan);
	fftw_plan p;
	p = fftw_plan_dft_1d(nChan, inP, inP, FFTW_FORWARD, FFTW_MEASURE);

	// Integer delay, FFT, and fractional delay
	long secStart, fftIndex;
	float complex *a;
	double *b;
	a = (float complex *) data->data;
	b = (double *) dataF->data;
	
	// Time-domain blanking control
	double cleanFactor;
	long nActFFT[nStand];
	
	#ifdef _OPENMP
		#ifdef _MKL
			fftw3_mkl.number_of_user_threads = omp_get_num_threads();
		#endif
		#pragma omp parallel default(shared) private(in, secStart, tempFB, tempB, i, j, k, m, fftIndex, cleanFactor)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			nActFFT[i] = 0;
			in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nChan);
			
			for(j=0; j<nFFT; j+=nTaps) {
				cleanFactor = 1.0;
				cblas_zdscal(nChan, 0.0, tempFB, 1);
				cblas_dscal(nChan, 0.0, tempB, 1);

				for(m=0; m<nTaps; m++) {
					secStart = nSamps * i + nChan*(j+m)/Overlap;
					
					for(k=0; k<nChan; k++) {
						in[k][0] = creal(*(a + secStart + k)) * fbWindow[nChan*m + k];
						in[k][1] = cimag(*(a + secStart + k)) * fbWindow[nChan*m + k];
						
						if( Clip && cabs(*(a + secStart + k)) >= Clip ) {
							cleanFactor = 0.0;
						}
					}
				
					fftw_execute_dft(p, in, in);
				
					for(k=0; k<nChan; k++) {
						fftIndex = (k + nChan/2) % nChan;
						tempFB[k] += cleanFactor*in[fftIndex][0] + imaginary * cleanFactor*in[fftIndex][1];
					}
				}
				
				#ifdef _MKL
					vzAbs(nChan, tempFB, tempB);
					vdSqr(nChan, tempB, tempB);
					vdAdd(nChan, (b+nChan*i), tempB, (b+nChan*i));
				#else
					for(k=0; k<nChan; k++) {
						*(b + nChan*i + k) += pow(cabs(tempFB[k]), 2);
					}
				#endif
				
				nActFFT[i] += (long) cleanFactor;
			}
			fftw_free(in);
		}
	}
	fftw_destroy_plan(p);
	fftw_free(inP);

	// cblas_dscal(nChan*nStand, ((float) nTaps)/(nChan*nFFT), b, 1);
	for(i=0; i<nStand; i++) {
		cblas_dscal(nChan, ((float) nTaps)/(nChan*nActFFT[i]), (b + i*nChan), 1);
	}
	
	Py_XDECREF(data);
	Py_XDECREF(windowData);

	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(PPSDC3_doc, \
"Perform a series of filter bank transforms with windows on complex-valued\n\
data to get the PSD.\n\
\n\
Input arguments are:\n\
 * signals: 2-D numpy.complex64 (stands by samples) array of data to FFT\n\
\n\
Input keywords are:\n\
 * LFFT: number of FFT channels to make (default=64)\n\
 * Overlap: number of overlapped FFTs to use (default=1)\n\
 * window: Callable Python function for generating the window\n\
 * ClipLevel: count value of 'bad' data.  FFT windows with instantaneous powers\n\
   greater than or equal to this value greater are zeroed.  Setting the ClipLevel\n\
   to zero disables time-domain blanking\n\
\n\
Outputs:\n\
 * psd: 2-D numpy.double (stands by channels) of PSD data\n\
");


/*
  Module Setup - Function Definitions and Documentation
*/

static PyMethodDef SpecMethods[] = {
	{"FPSDR2",  (PyCFunction) FPSDR2,  METH_VARARGS|METH_KEYWORDS, FPSDR2_doc}, 
	{"FPSDR3",  (PyCFunction) FPSDR3,  METH_VARARGS|METH_KEYWORDS, FPSDR3_doc}, 
	{"FPSDC2",  (PyCFunction) FPSDC2,  METH_VARARGS|METH_KEYWORDS, FPSDC2_doc}, 
	{"FPSDC3",  (PyCFunction) FPSDC3,  METH_VARARGS|METH_KEYWORDS, FPSDC3_doc}, 
	{"PPSDR2",  (PyCFunction) PPSDR2,  METH_VARARGS|METH_KEYWORDS, PPSDR2_doc}, 
	{"PPSDR3",  (PyCFunction) PPSDR3,  METH_VARARGS|METH_KEYWORDS, PPSDR3_doc}, 
	{"PPSDC2",  (PyCFunction) PPSDC2,  METH_VARARGS|METH_KEYWORDS, PPSDC2_doc}, 
	{"PPSDC3",  (PyCFunction) PPSDC3,  METH_VARARGS|METH_KEYWORDS, PPSDC3_doc}, 
	{NULL,      NULL,    0,                          NULL      }
};

PyDoc_STRVAR(spec_doc, \
"Extension to take timeseries data and convert it to the frequency domain.\n\
\n\
The functions defined in this module are:\n\
  * FPSDR2 -  FFT and integrate function for computing a series of overlapped\n\
    Fourier transforms for a real-valued (TBW) signal from a collection of\n\
    stands all at once.\n\
  * FPSDR3 - Similar to FPSDR2, but allows for a window function to be applied\n\
    to the data.\n\
  * FPSDC2 - FFT and integrate function for computing a series of overlapped\n\
    Fourier transforms for a complex-valued (TBN and DRX) signal from a \n\
    collection of stands/beams all at once.\n\
  * FPSDC3 - Similar to FPSDC2, but allows for a window function to be applied\n\
    to the data.\n\
  * PPSDR2 - Polyphase filterband (PFB) and integration function for computing\n\
    a series of overlapped PFB transforms for a real-valued (TBW) signal from\n\
    a collection of stands all at once.\n\
  * PPSDR3 - Similar to PPSDR2, but allows for a window function to be applied\n\
    to the data.\n\
  * PPSDC2 - PFB and integrate function for computing a series of overlapped\n\
    PFB transforms for a complex-valued (TBN and DRX) signal from a collection\n\
    of stands/beams all at once.\n\
  * PPSDC3 - Similar to PPSDC2, but allows for a window function to be applied\n\
    to the data.\n\
\n\
See the inidividual functions for more details.");


/*
  Module Setup - Initialization
*/

PyMODINIT_FUNC init_spec(void) {
	PyObject *m;

	// Module definitions and functions
	m = Py_InitModule3("_spec", SpecMethods, spec_doc);
	import_array();
	
	// Version and revision information
	PyModule_AddObject(m, "__version__", PyString_FromString("0.3"));
	PyModule_AddObject(m, "__revision__", PyString_FromString("$Rev$"));
	
}

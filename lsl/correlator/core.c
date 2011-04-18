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
  FFT Functions ("F-engines")
    1. FEngineR2 - FFT a real-valued collection of signals
    2. FEngineR3 - window the data and FFT a real-valued collection of signals
    3. FEngineC2 - FFT a complex-valued collection of signals
    4. FEngineC3 - window the data and FFT a complex-valued collection of signals
*/


static PyObject *FEngineR2(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *freq, *delays, *signalsF;
	PyArrayObject *data, *fq, *times, *dataF;
	int nChan = 64;
	int Overlap = 1;
	double SampleRate = 196.0e6;

	long i, j, k, nStand, nSamps, nFFT;
	
	static char *kwlist[] = {"signals", "freq", "delays", "LFFT", "Overlap", "SampleRate", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OOO|iid", kwlist, &signals, &freq, &delays, &nChan, &Overlap, &SampleRate)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it useable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_INT16, 2, 2);
	fq = (PyArrayObject *) PyArray_ContiguousFromObject(freq, NPY_DOUBLE, 1, 1);
	times = (PyArrayObject *) PyArray_ContiguousFromObject(delays, NPY_DOUBLE, 2, 2);
	
	// Check data dimensions
	if(data->dimensions[0] != times->dimensions[0]) {
		PyErr_Format(PyExc_TypeError, "signals and delays have different stand counts");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		return NULL;
	}
	
	if(nChan != (fq->dimensions[0]+1)) {
		PyErr_Format(PyExc_RuntimeError, "freq has a different channel count than nChan");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		return NULL;
	}
	
	if(fq->dimensions[0] != times->dimensions[1]) {
		PyErr_Format(PyExc_TypeError, "freq and delays have different channel counts");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		return NULL;
	}
	
	// Get the properties of the data
	nStand = (long) data->dimensions[0];
	nSamps = (long) data->dimensions[1];
	
	// Compute the interger sample offset and the fractional sample delay for each stand
	npy_intp *tLoc;
	long start[nStand];
	long startMax = 0;
	double frac[nStand][nChan];
	tLoc = PyDimMem_NEW(2);
	for(i=0; i<nStand; i++) {
		tLoc[0] = (npy_intp) i;
		tLoc[1] = (npy_intp) (nChan / 2);
		start[i] = (long) round(*(double *) PyArray_GetPtr(times, tLoc) * SampleRate);
		if(start[i] > startMax) {
			startMax = start[i];
		}

		for(j=0; j<nChan; j++) {
			tLoc[1] = (npy_intp) j;
			frac[i][j] = *(double *) PyArray_GetPtr(times, tLoc) - (double) start[i]/SampleRate;
		}
	}
	PyDimMem_FREE(tLoc);

	// Find out how large the output array needs to be and initialize it
	nFFT = (nSamps - startMax) / 2 / nChan * Overlap - Overlap + 1;
	npy_intp dims[3];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) (nChan - 1);
	dims[2] = (npy_intp) nFFT;
	dataF = (PyArrayObject*) PyArray_SimpleNew(3, dims, NPY_CDOUBLE);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		return NULL;
	}
	
	// Create the FFTW plan                          
	fftw_complex *inP, *in;                          
	inP = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * 2*nChan);
	fftw_plan p;
	p = fftw_plan_dft_1d(2*nChan, inP, inP, FFTW_FORWARD, FFTW_MEASURE);
	
	// Integer delay, FFT, and fractional delay
	long secStart, fftIndex;
	npy_intp *dLoc, *fLoc, *qLoc;
	
	#ifdef _OPENMP
		#ifdef _MKL
			fftw3_mkl.number_of_user_threads = omp_get_num_threads();
		#endif
		#pragma omp parallel default(shared) private(dLoc, fLoc, qLoc, in, secStart, j, k, fftIndex)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			dLoc = PyDimMem_NEW(2);
			fLoc = PyDimMem_NEW(3);
			qLoc = PyDimMem_NEW(1);
			
			in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * 2*nChan);
			
			dLoc[0] = (npy_intp) i;
			fLoc[0] = (npy_intp) i;
			
			for(j=0; j<nFFT; j++) {
				secStart = start[i] + ((long) (2*nChan*((float) j)/Overlap));
				
				for(k=0; k<2*nChan; k++) {
					dLoc[1] = (npy_intp) (secStart + k);
					in[k][0] = *(double *) PyArray_GetPtr(data, dLoc);
					in[k][1] = 0.0;
				}	
			
				fftw_execute_dft(p, in, in);
			
				fLoc[2] = (npy_intp) j;
				for(k=0; k<(nChan-1); k++) {
					fLoc[1] = (npy_intp) k;
					qLoc[0] = (npy_intp) k;
					fftIndex = k + 1;
					*(double complex *) PyArray_GetPtr(dataF, fLoc) = (in[fftIndex][0] + imaginary*in[fftIndex][1]);
					*(double complex *) PyArray_GetPtr(dataF, fLoc) *= cexp(-2*imaginary*PI* *(double *) PyArray_GetPtr(fq, qLoc) * frac[i][k]);
				}
			}
			
			PyDimMem_FREE(dLoc);
			PyDimMem_FREE(fLoc);
			PyDimMem_FREE(qLoc);

			fftw_free(in);
		}
	}
	fftw_destroy_plan(p);
	fftw_free(inP);

	Py_XDECREF(data);
	Py_XDECREF(fq);
	Py_XDECREF(times);

	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(FEngineR2_doc, \
"Perform a series of overlaped Fourier transforms on real-valued data using\n\
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
\n\
Outputs:\n\
 * fsignals: 3-D numpy.cdouble (stands by channels by FFT_set) of FFTd\n\
   data\n\
");


static PyObject *FEngineR3(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *freq, *delays, *signalsF, *window;
	PyArrayObject *data, *fq, *times, *dataF, *windowData;
	int nChan = 64;
	int Overlap = 1;
	double SampleRate = 196.0e6;

	long i, j, k, nStand, nSamps, nFFT;
	
	static char *kwlist[] = {"signals", "freq", "delays", "LFFT", "Overlap", "SampleRate", "window", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OOO|iidO:set_callback", kwlist, &signals, &freq, &delays, &nChan, &Overlap, &SampleRate, &window)) {
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

	// Bring the data into C and make it useable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_INT16, 2, 2);
	fq = (PyArrayObject *) PyArray_ContiguousFromObject(freq, NPY_DOUBLE, 1, 1);
	times = (PyArrayObject *) PyArray_ContiguousFromObject(delays, NPY_DOUBLE, 2, 2);
	
	// Calculate the windowing function
	window = Py_BuildValue("(i)", 2*nChan);
	window = PyObject_CallObject(windowFunc, window);
	windowData = (PyArrayObject *) PyArray_ContiguousFromObject(window, NPY_DOUBLE, 1, 1);
	Py_DECREF(window);
	
	// Check data dimensions
	if(data->dimensions[0] != times->dimensions[0]) {
		PyErr_Format(PyExc_TypeError, "signals and delays have different stand counts");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		Py_XDECREF(windowData);
		return NULL;
	}
	
	if(nChan != (fq->dimensions[0]+1)) {
		PyErr_Format(PyExc_RuntimeError, "freq has a different channel count than nChan");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		Py_XDECREF(windowData);
		return NULL;
	}
	
	if(fq->dimensions[0] != times->dimensions[1]) {
		PyErr_Format(PyExc_TypeError, "freq and delays have different channel counts");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		Py_XDECREF(windowData);
		return NULL;
	}
	
	// Get the properties of the data
	nStand = (long) data->dimensions[0];
	nSamps = (long) data->dimensions[1];
	
	// Compute the interger sample offset and the fractional sample delay for each stand
	npy_intp *tLoc;
	long start[nStand];
	long startMax = 0;
	double frac[nStand][nChan];
	tLoc = PyDimMem_NEW(2);
	for(i=0; i<nStand; i++) {
		tLoc[0] = (npy_intp) i;
		tLoc[1] = (npy_intp) (nChan / 2);
		start[i] = (long) round(*(double *) PyArray_GetPtr(times, tLoc) * SampleRate);
		if(start[i] > startMax) {
			startMax = start[i];
		}

		for(j=0; j<nChan; j++) {
			tLoc[1] = (npy_intp) j;
			frac[i][j] = *(double *) PyArray_GetPtr(times, tLoc) - (double) start[i]/SampleRate;
		}
	}
	PyDimMem_FREE(tLoc);

	// Find out how large the output array needs to be and initialize it
	nFFT = (nSamps - startMax) / 2 / nChan * Overlap - Overlap + 1;
	npy_intp dims[3];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) (nChan - 1);
	dims[2] = (npy_intp) nFFT;
	dataF = (PyArrayObject*) PyArray_SimpleNew(3, dims, NPY_CDOUBLE);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		Py_XDECREF(windowData);
		return NULL;
	}
	
	// Create the FFTW plan                          
	fftw_complex *inP, *in;                          
	inP = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * 2*nChan);
	fftw_plan p;
	p = fftw_plan_dft_1d(2*nChan, inP, inP, FFTW_FORWARD, FFTW_MEASURE);
	
	// Integer delay, FFT, and fractional delay
	long secStart, fftIndex;
	npy_intp *dLoc, *fLoc, *qLoc;
	
	#ifdef _OPENMP
		#ifdef _MKL
			fftw3_mkl.number_of_user_threads = omp_get_num_threads();
		#endif
		#pragma omp parallel default(shared) private(dLoc, fLoc, qLoc, in, secStart, j, k, fftIndex)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			dLoc = PyDimMem_NEW(2);
			fLoc = PyDimMem_NEW(3);
			qLoc = PyDimMem_NEW(1);
			
			in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * 2*nChan);
			
			dLoc[0] = (npy_intp) i;
			fLoc[0] = (npy_intp) i;
			
			for(j=0; j<nFFT; j++) {
				secStart = start[i] + ((long) (2*nChan*((float) j)/Overlap));
				
				for(k=0; k<2*nChan; k++) {
					dLoc[1] = (npy_intp) (secStart + k);
					qLoc[0] = (npy_intp) k;
					in[k][0] = *(double *) PyArray_GetPtr(data, dLoc) * *(double *) PyArray_GetPtr(windowData, qLoc);
					in[k][1] = 0.0;
				}	
			
				fftw_execute_dft(p, in, in);
			
				fLoc[2] = (npy_intp) j;
				for(k=0; k<(nChan-1); k++) {
					fLoc[1] = (npy_intp) k;
					qLoc[0] = (npy_intp) k;
					fftIndex = k + 1;
					*(double complex *) PyArray_GetPtr(dataF, fLoc) = (in[fftIndex][0] + imaginary*in[fftIndex][1]);
					*(double complex *) PyArray_GetPtr(dataF, fLoc) *= cexp(-2*imaginary*PI* *(double *) PyArray_GetPtr(fq, qLoc) * frac[i][k]);
				}
			}
			
			PyDimMem_FREE(dLoc);
			PyDimMem_FREE(fLoc);
			PyDimMem_FREE(qLoc);

			fftw_free(in);
		}
	}
	fftw_destroy_plan(p);
	fftw_free(inP);

	Py_XDECREF(data);
	Py_XDECREF(fq);
	Py_XDECREF(times);
	Py_XDECREF(windowData);

	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(FEngineR3_doc, \
"Perform a series of overlaped Fourier transforms on real-valued data using\n\
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
\n\
Outputs:\n\
 * fsignals: 3-D numpy.cdouble (stands by channels by FFT_set) of FFTd\n\
   data\n\
");


static PyObject *FEngineC2(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *freq, *delays, *signalsF;
	PyArrayObject *data, *fq, *times, *dataF;
	int nChan = 64;
	int Overlap = 1;
	double SampleRate = 1.0e5;

	long i, j, k, nStand, nSamps, nFFT;

	static char *kwlist[] = {"signals", "freq", "delays", "LFFT", "Overlap", "SampleRate", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OOO|iid", kwlist, &signals, &freq, &delays, &nChan, &Overlap, &SampleRate)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it useable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_COMPLEX64, 2, 2);
	fq = (PyArrayObject *) PyArray_ContiguousFromObject(freq, NPY_DOUBLE, 1, 1);
	times = (PyArrayObject *) PyArray_ContiguousFromObject(delays, NPY_DOUBLE, 2, 2);
	
	// Check data dimensions
	if(data->dimensions[0] != times->dimensions[0]) {
		PyErr_Format(PyExc_TypeError, "signals and delays have different stand counts");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		return NULL;
	}
	
	if(nChan != (fq->dimensions[0]+1)) {
		PyErr_Format(PyExc_RuntimeError, "freq has a different channel count than nChan");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		return NULL;
	}
	
	if(fq->dimensions[0] != times->dimensions[1]) {
		PyErr_Format(PyExc_TypeError, "freq and delays have different channel counts");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		return NULL;
	}

	// Get the properties of the data
	nStand = (long) data->dimensions[0];
	nSamps = (long) data->dimensions[1];
	
	// Compute the interger sample offset and the fractional sample delay for each stand
	npy_intp *tLoc;
	long start[nStand];
	long startMax = 0;
	double frac[nStand][nChan];
	tLoc = PyDimMem_NEW(2);
	for(i=0; i<nStand; i++) {
		tLoc[0] = (npy_intp) i;
		tLoc[1] = (npy_intp) (nChan / 2);
		start[i] = (long) round(*(double *) PyArray_GetPtr(times, tLoc) * SampleRate);
		if(start[i] > startMax) {
			startMax = start[i];
		}

		for(j=0; j<nChan; j++) {
			tLoc[1] = (npy_intp) j;
			frac[i][j] = *(double *) PyArray_GetPtr(times, tLoc) - (double) start[i]/SampleRate;
		}
	}
	PyDimMem_FREE(tLoc);

	// Find out how large the output array needs to be and initialize it
	nFFT = (nSamps - startMax) / nChan * Overlap - Overlap + 1;
	npy_intp dims[3];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) (nChan - 1);
	dims[2] = (npy_intp) nFFT;
	dataF = (PyArrayObject*) PyArray_SimpleNew(3, dims, NPY_CDOUBLE);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		return NULL;
	}

	// Create the FFTW plan
	fftw_complex *inP, *in;
	inP = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nChan);
	fftw_plan p;
	p = fftw_plan_dft_1d(nChan, inP, inP, FFTW_FORWARD, FFTW_MEASURE);

	// Integer delay, FFT, and fractional delay
	long secStart, fftIndex;
	npy_intp *dLoc, *fLoc, *qLoc;
	
	#ifdef _OPENMP
		#ifdef _MKL
			fftw3_mkl.number_of_user_threads = omp_get_num_threads();
		#endif
		#pragma omp parallel default(shared) private(dLoc, fLoc, qLoc, in, secStart, j, k, fftIndex)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			dLoc = PyDimMem_NEW(2);
			fLoc = PyDimMem_NEW(3);
			qLoc = PyDimMem_NEW(1);
			
			in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nChan);
			
			dLoc[0] = (npy_intp) i;
			fLoc[0] = (npy_intp) i;
			
			for(j=0; j<nFFT; j++) {
				secStart = start[i] + ((long) (nChan*((float) j)/Overlap));
				
				for(k=0; k<nChan; k++) {
					dLoc[1] = (npy_intp) (secStart + k);
					in[k][0] = creal(*(double complex *) PyArray_GetPtr(data, dLoc));
					in[k][1] = cimag(*(double complex *) PyArray_GetPtr(data, dLoc));
				}	
			
				fftw_execute_dft(p, in, in);
			
				fLoc[2] = (npy_intp) j;
				for(k=0; k<(nChan-1); k++) {
					fLoc[1] = (npy_intp) k;
					qLoc[0] = (npy_intp) k;
					fftIndex = ((k+1) + nChan/2) % nChan;
					*(double complex *) PyArray_GetPtr(dataF, fLoc) = (in[fftIndex][0] + imaginary*in[fftIndex][1]);
					*(double complex *) PyArray_GetPtr(dataF, fLoc) *= cexp(-2*imaginary*PI* *(double *) PyArray_GetPtr(fq, qLoc) * frac[i][k]);
				}
			}
			
			PyDimMem_FREE(dLoc);
			PyDimMem_FREE(fLoc);
			PyDimMem_FREE(qLoc);

			fftw_free(in);
		}
	}
	fftw_destroy_plan(p);
	fftw_free(inP);

	Py_XDECREF(data);
	Py_XDECREF(fq);
	Py_XDECREF(times);

	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(FEngineC2_doc, \
"Perform a series of overlaped Fourier transforms on complex-valued data\n\
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
\n\
Outputs:\n\
 * fsignals: 3-D numpy.cdouble (stands by channels by FFT_set) of FFTd\n\
   data\n\
");


static PyObject *FEngineC3(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *freq, *delays, *signalsF, *window;
	PyArrayObject *data, *fq, *times, *dataF, *windowData;
	int nChan = 64;
	int Overlap = 1;
	double SampleRate = 1.0e5;

	long i, j, k, nStand, nSamps, nFFT;

	static char *kwlist[] = {"signals", "freq", "delays", "LFFT", "Overlap", "SampleRate", "window", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OOO|iidO:set_callback", kwlist, &signals, &freq, &delays, &nChan, &Overlap, &SampleRate, &window)) {
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

	// Bring the data into C and make it useable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_COMPLEX64, 2, 2);
	fq = (PyArrayObject *) PyArray_ContiguousFromObject(freq, NPY_DOUBLE, 1, 1);
	times = (PyArrayObject *) PyArray_ContiguousFromObject(delays, NPY_DOUBLE, 2, 2);
	
	// Calculate the windowing function
	window = Py_BuildValue("(i)", nChan);
	window = PyObject_CallObject(windowFunc, window);
	windowData = (PyArrayObject *) PyArray_ContiguousFromObject(window, NPY_DOUBLE, 1, 1);
	Py_DECREF(window);
	
	// Check data dimensions
	if(data->dimensions[0] != times->dimensions[0]) {
		PyErr_Format(PyExc_TypeError, "signals and delays have different stand counts");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		Py_XDECREF(windowData);
		return NULL;
	}
	
	if(nChan != (fq->dimensions[0]+1)) {
		PyErr_Format(PyExc_RuntimeError, "freq has a different channel count than nChan");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		Py_XDECREF(windowData);
		return NULL;
	}
	
	if(fq->dimensions[0] != times->dimensions[1]) {
		PyErr_Format(PyExc_TypeError, "freq and delays have different channel counts");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		Py_XDECREF(windowData);
		return NULL;
	}

	// Get the properties of the data
	nStand = (long) data->dimensions[0];
	nSamps = (long) data->dimensions[1];
	
	// Compute the interger sample offset and the fractional sample delay for each stand
	npy_intp *tLoc;
	long start[nStand];
	long startMax = 0;
	double frac[nStand][nChan];
	tLoc = PyDimMem_NEW(2);
	for(i=0; i<nStand; i++) {
		tLoc[0] = (npy_intp) i;
		tLoc[1] = (npy_intp) (nChan / 2);
		start[i] = (long) round(*(double *) PyArray_GetPtr(times, tLoc) * SampleRate);
		if(start[i] > startMax) {
			startMax = start[i];
		}

		for(j=0; j<nChan; j++) {
			tLoc[1] = (npy_intp) j;
			frac[i][j] = *(double *) PyArray_GetPtr(times, tLoc) - (double) start[i]/SampleRate;
		}
	}
	PyDimMem_FREE(tLoc);

	// Find out how large the output array needs to be and initialize it
	nFFT = (nSamps - startMax) / nChan * Overlap - Overlap + 1;
	npy_intp dims[3];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) (nChan - 1);
	dims[2] = (npy_intp) nFFT;
	dataF = (PyArrayObject*) PyArray_SimpleNew(3, dims, NPY_CDOUBLE);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		Py_XDECREF(windowData);
		return NULL;
	}

	// Create the FFTW plan
	fftw_complex *inP, *in;
	inP = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nChan);
	fftw_plan p;
	p = fftw_plan_dft_1d(nChan, inP, inP, FFTW_FORWARD, FFTW_MEASURE);

	// Integer delay, FFT, and fractional delay
	long secStart, fftIndex;
	npy_intp *dLoc, *fLoc, *qLoc;
	
	#ifdef _OPENMP
		#ifdef _MKL
			fftw3_mkl.number_of_user_threads = omp_get_num_threads();
		#endif
		#pragma omp parallel default(shared) private(dLoc, fLoc, qLoc, in, secStart, j, k, fftIndex)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			dLoc = PyDimMem_NEW(2);
			fLoc = PyDimMem_NEW(3);
			qLoc = PyDimMem_NEW(1);
			
			in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nChan);
			
			dLoc[0] = (npy_intp) i;
			fLoc[0] = (npy_intp) i;
			
			for(j=0; j<nFFT; j++) {
				secStart = start[i] + ((long) (nChan*((float) j)/Overlap));
				
				for(k=0; k<nChan; k++) {
					dLoc[1] = (npy_intp) (secStart + k);
					qLoc[0] = (npy_intp) k;
					in[k][0] = creal(*(double complex *) PyArray_GetPtr(data, dLoc) * *(double *) PyArray_GetPtr(windowData, qLoc));
					in[k][1] = cimag(*(double complex *) PyArray_GetPtr(data, dLoc) * *(double *) PyArray_GetPtr(windowData, qLoc));
				}	
			
				fftw_execute_dft(p, in, in);
			
				fLoc[2] = (npy_intp) j;
				for(k=0; k<(nChan-1); k++) {
					fLoc[1] = (npy_intp) k;
					qLoc[0] = (npy_intp) k;
					fftIndex = ((k+1) + nChan/2) % nChan;
					*(double complex *) PyArray_GetPtr(dataF, fLoc) = (in[fftIndex][0] + imaginary*in[fftIndex][1]);
					*(double complex *) PyArray_GetPtr(dataF, fLoc) *= cexp(-2*imaginary*PI* *(double *) PyArray_GetPtr(fq, qLoc) * frac[i][k]);
				}
			}
			
			PyDimMem_FREE(dLoc);
			PyDimMem_FREE(fLoc);
			PyDimMem_FREE(qLoc);

			fftw_free(in);
		}
	}
	fftw_destroy_plan(p);
	fftw_free(inP);

	Py_XDECREF(data);
	Py_XDECREF(fq);
	Py_XDECREF(times);
	Py_XDECREF(windowData);

	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(FEngineC3_doc, \
"Perform a series of overlaped Fourier transforms on complex-valued data\n\
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
\n\
Outputs:\n\
 * fsignals: 3-D numpy.cdouble (stands by channels by FFT_set) of FFTd\n\
   data\n\
");


/*
  Cross-Multiplication And Accumulation Function ("X Engines")
    1. XEngine  - XMAC two signals
    2. XEngine2 -  XMAC two collections of signals
*/

static PyObject *XEngine(PyObject *self, PyObject *args) {
	PyObject *signal1, *signal2, *output;
	PyArrayObject *data1, *data2, *vis;
	long nChan, nFFT1, nFFT2, nFFT;

	if(!PyArg_ParseTuple(args, "OO", &signal1, &signal2)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it useable
	data1 = (PyArrayObject *) PyArray_ContiguousFromObject(signal1, NPY_CDOUBLE, 2, 2);
	data2 = (PyArrayObject *) PyArray_ContiguousFromObject(signal2, NPY_CDOUBLE, 2, 2);

	// Check data dimensions
	if(data1->dimensions[0] != data2->dimensions[0]) {
		PyErr_Format(PyExc_TypeError, "signal1 and signal2 have different channel counts");                                                           
		Py_XDECREF(data1);
		Py_XDECREF(data2);
		return NULL;                       
	}

	// Get channel count and number of FFTs stored
	nChan = (long) data1->dimensions[0];
	nFFT1 = (long) data1->dimensions[1];
	nFFT2 = (long) data2->dimensions[1];
	if(nFFT1 < nFFT2) {
		nFFT = nFFT1;
	} else {
		nFFT = nFFT2;
	}

	// Create the output visibility array
	npy_intp dims[1];
	dims[0] = (npy_intp) nChan;
	vis = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
	if(vis == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data1);
		Py_XDECREF(data2);
		return NULL;
	}

	// Cross-multiplication and accumulation
	long c, f;
	npy_intp *vLoc, *dLoc;
	vLoc = PyDimMem_NEW(1);
	dLoc = PyDimMem_NEW(2);
	for(c=0; c<nChan; c++) {
		vLoc[0] = (npy_intp) c;
		dLoc[0] = (npy_intp) c;
		*(double complex *) PyArray_GetPtr(vis, vLoc) = 0.0;
		for(f=0; f<nFFT; f++) {
			dLoc[1] = (npy_intp) f;
			*(double complex *) PyArray_GetPtr(vis, vLoc) += *(double complex *) PyArray_GetPtr(data1, dLoc) * *(double complex *) PyArray_GetPtr(data2, dLoc) / nFFT;
		}
	}
	PyDimMem_FREE(vLoc);
	PyDimMem_FREE(dLoc);
	
	Py_XDECREF(data1);
	Py_XDECREF(data2);

	output = Py_BuildValue("O", PyArray_Return(vis));
	Py_XDECREF(vis);

	return output;
}

PyDoc_STRVAR(XEngine_doc, \
"Perform XMAC on two data streams out of the F engine.\n\
\n\
Input arguments are:\n\
 * fsignals1: 3-D numpy.cdouble (stand by channels by FFT_set) array of FFTd\n\
   data from an F engine.\n\
 * fsignals2: 3-D numpy.cdouble (stand by channels by FFT_set) array of\n\
   conjudated FFTd data from an F engine.\n\
\n\
Ouputs:\n\
  * visibility: 3-D numpy.cdouble (baseline by channel) array of cross-\n\
    correlated and average visibility data.\n\
\n\
.. note::\n\
\tThis function is *slower* than a pure numpy version of the same function.\n\
");


static PyObject *XEngine2(PyObject *self, PyObject *args) {
	PyObject *signals, *signalsC, *output;
	PyArrayObject *data, *dataC, *vis;
	long nStand, nChan, nFFT, nBL;	

	if(!PyArg_ParseTuple(args, "OO", &signals, &signalsC)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it useable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_CDOUBLE, 3, 3);
	dataC = (PyArrayObject *) PyArray_ContiguousFromObject(signalsC, NPY_CDOUBLE, 3, 3);

	// Get channel count and number of FFTs stored
	nStand = (long) data->dimensions[0];
	nChan = (long) data->dimensions[1];
	nFFT = (long) data->dimensions[2];
	nBL = (nStand+1)*nStand/2;
	
	// Create the output visibility array and fill with zeros
	npy_intp dims[2];
	dims[0] = (npy_intp) nBL;
	dims[1] = (npy_intp) nChan;
	vis = (PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_CDOUBLE);
	if(vis == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		return NULL;
	}
	PyArray_FILLWBYTE(vis, 0);
	
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
	npy_intp *vLoc, *dLoc1, *dLoc2;
	double complex tempVis, dataBL1[nFFT], dataBL2[nFFT];
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(vLoc, dLoc1, dLoc2, c, f, dataBL1, dataBL2, tempVis)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(bl=0; bl<nBL; bl++) {
			vLoc = PyDimMem_NEW(2);
			dLoc1 = PyDimMem_NEW(3);
			dLoc2 = PyDimMem_NEW(3);
			
			vLoc[0] = (npy_intp) bl;
			dLoc1[0] = (npy_intp) mapper[bl][0];
			dLoc2[0] = (npy_intp) mapper[bl][1];
			
			for(c=0; c<nChan; c++) {
				vLoc[1] = (npy_intp) c;
				dLoc1[1] = (npy_intp) c;
				dLoc2[1] = (npy_intp) c;
				
				for(f=0; f<nFFT; f++) {
					dLoc1[2] = (npy_intp) f;
					dLoc2[2] = (npy_intp) f;
					
					dataBL1[f] = *(double complex *) PyArray_GetPtr(data, dLoc1);
					dataBL2[f] = *(double complex *) PyArray_GetPtr(dataC, dLoc2);
				}
			
				cblas_zdotu_sub(nFFT, dataBL1, 1, dataBL2, 1, &tempVis);
				*(double complex *) PyArray_GetPtr(vis, vLoc) = tempVis / nFFT;
			}
			
			PyDimMem_FREE(vLoc);
			PyDimMem_FREE(dLoc1);
			PyDimMem_FREE(dLoc2);
		}
	}
	Py_XDECREF(data);
	Py_XDECREF(dataC);

	output = Py_BuildValue("O", PyArray_Return(vis));
	Py_XDECREF(vis);

	return output;
}

PyDoc_STRVAR(XEngine2_doc, \
"Perform all XMACs for a data stream out of the F engine using OpenMP.\n\
\n\
Input arguments are:\n\
 * fsignals1: 3-D numpy.cdouble (stand by channels by FFT_set) array of FFTd\n\
   data from an F engine.\n\
 * fsignals2: 3-D numpy.cdouble (stand by channels by FFT_set) array of\n\
   conjudated FFTd data from an F engine.\n\
\n\
Ouputs:\n\
  * visibility: 3-D numpy.cdouble (baseline by channel) array of cross-\n\
    correlated and average visibility data.\n\
");
    
    
/* 
  Polyphase Filterbank Functions ("P Engines")
    1. PEngineR2 - Filterbank analog to FEngineR2 using 4 taps
    2. PEngineR3 - Filterbank analog to FEngineR3 using 4 taps
    3. PEngineC2 - Filterbank analog to FEnginedC2 using 4 taps
*/

static PyObject *PEngineR2(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *freq, *delays, *signalsF;
	PyArrayObject *data, *fq, *times, *dataF;
	int nChan = 64;
	int Overlap = 1;
	double SampleRate = 196.0e6;

	long i, j, k, nStand, nSamps, nFFT;
	
	static char *kwlist[] = {"signals", "freq", "delays", "LFFT", "Overlap", "SampleRate", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OOO|iid", kwlist, &signals, &freq, &delays, &nChan, &Overlap, &SampleRate)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it useable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, PyArray_DOUBLE, 2, 2);
	fq = (PyArrayObject *) PyArray_ContiguousFromObject(freq, PyArray_DOUBLE, 1, 1);
	times = (PyArrayObject *) PyArray_ContiguousFromObject(delays, PyArray_DOUBLE, 2, 2);
	
	// Check data dimensions
	if(data->dimensions[0] != times->dimensions[0]) {
		PyErr_Format(PyExc_TypeError, "signals and delays have different stand counts");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		return NULL;
	}
	
	if(nChan != (fq->dimensions[0]+1)) {
		PyErr_Format(PyExc_RuntimeError, "freq has a different channel count than nChan");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		return NULL;
	}
	
	if(fq->dimensions[0] != times->dimensions[1]) {
		PyErr_Format(PyExc_TypeError, "freq and delays have different channel counts");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		return NULL;
	}
	
	// Get the properties of the data
	nStand = (long) data->dimensions[0];
	nSamps = (long) data->dimensions[1];
	
	// Compute the interger sample offset and the fractional sample delay for each stand
	npy_intp *tLoc;
	long start[nStand];
	long startMax = 0;
	double frac[nStand][nChan];
	tLoc = PyDimMem_NEW(2);
	for(i=0; i<nStand; i++) {
		tLoc[0] = (npy_intp) i;
		tLoc[1] = (npy_intp) (nChan / 2);
		start[i] = (long) round(*(double *) PyArray_GetPtr(times, tLoc) * SampleRate);
		if(start[i] > startMax) {
			startMax = start[i];
		}

		for(j=0; j<nChan; j++) {
			tLoc[1] = (npy_intp) j;
			frac[i][j] = *(double *) PyArray_GetPtr(times, tLoc) - (double) start[i]/SampleRate;
		}
	}
	PyDimMem_FREE(tLoc);

	// Compute the filterbank window for the correct numer of taps
	double fbWindow[2*nChan*nTaps];
	for(i=0; i<2*nChan*nTaps; i++) {
		fbWindow[i] = sinc((double) (i - nTaps*nChan + 0.5)/2/nChan);
	}

	// Find out how large the output array needs to be and initialize it
	nFFT = (nSamps - startMax) / 2 / nChan * Overlap - Overlap + 1;
	nFFT = (nFFT / nTaps) * nTaps;
	npy_intp dims[3];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) (nChan - 1);
	dims[2] = (npy_intp) (nFFT / nTaps);
	dataF = (PyArrayObject*) PyArray_SimpleNew(3, dims, PyArray_CDOUBLE);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
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
	npy_intp *dLoc, *fLoc, *qLoc;
	
	#ifdef _OPENMP
		#ifdef _MKL
			fftw3_mkl.number_of_user_threads = omp_get_num_threads();
		#endif
		#pragma omp parallel default(shared) private(dLoc, fLoc, qLoc, in, secStart, j, k, fftIndex)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			dLoc = PyDimMem_NEW(2);
			fLoc = PyDimMem_NEW(3);
			qLoc = PyDimMem_NEW(1);
			
			in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * 2*nChan);
			
			dLoc[0] = (npy_intp) i;
			fLoc[0] = (npy_intp) i;

			for(j=0; j<nFFT; j++) {
				secStart = start[i] + ((long) (2*nChan*((float) j)/Overlap));
				for(k=0; k<2*nChan; k++) {
					dLoc[1] = (npy_intp) (secStart + k);
					in[k][0] = *(double *) PyArray_GetPtr(data, dLoc) * fbWindow[2*nChan*(j % nTaps) + k];
					in[k][1] = 0.0;
				}	
			
				fftw_execute_dft(p, in, in);
			
				fLoc[2] = (npy_intp) (j / nTaps);
				for(k=0; k<(nChan-1); k++) {
					fLoc[1] = (npy_intp) k;
					qLoc[0] = (npy_intp) k;
					fftIndex = k + 1;
					*(double complex *) PyArray_GetPtr(dataF, fLoc) += (in[fftIndex][0] + imaginary*in[fftIndex][1]) * cexp(-2*imaginary*PI* *(double *) PyArray_GetPtr(fq, qLoc) * frac[i][k]);
				}
			}
			
			PyDimMem_FREE(dLoc);
			PyDimMem_FREE(fLoc);
			PyDimMem_FREE(qLoc);

			fftw_free(in);
		}
	}
	fftw_destroy_plan(p);
	fftw_free(inP);

	Py_XDECREF(data);
	Py_XDECREF(fq);
	Py_XDECREF(times);

	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(PEngineR2_doc, "Perform a series of overlaped filter bank tranforms on real-valued data using OpenMP.");


static PyObject *PEngineR3(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *freq, *delays, *signalsF, *window;
	PyArrayObject *data, *fq, *times, *dataF, *windowData;
	int nChan = 64;
	int Overlap = 1;
	double SampleRate = 196.0e6;

	long i, j, k, nStand, nSamps, nFFT;
	
	static char *kwlist[] = {"signals", "freq", "delays", "LFFT", "Overlap", "SampleRate", "window", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OOO|iidO:set_callback", kwlist, &signals, &freq, &delays, &nChan, &Overlap, &SampleRate, &window)) {
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

	// Bring the data into C and make it useable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, PyArray_DOUBLE, 2, 2);
	fq = (PyArrayObject *) PyArray_ContiguousFromObject(freq, PyArray_DOUBLE, 1, 1);
	times = (PyArrayObject *) PyArray_ContiguousFromObject(delays, PyArray_DOUBLE, 2, 2);
	
	// Calculate the windowing function
	window = Py_BuildValue("(i)", 2*nTaps*nChan);
	window = PyObject_CallObject(windowFunc, window);
	windowData = (PyArrayObject *) PyArray_ContiguousFromObject(window, PyArray_DOUBLE, 1, 1);
	Py_DECREF(window);
	
	// Check data dimensions
	if(data->dimensions[0] != times->dimensions[0]) {
		PyErr_Format(PyExc_TypeError, "signals and delays have different stand counts");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		Py_XDECREF(windowData);
		return NULL;
	}
	
	if(nChan != (fq->dimensions[0]+1)) {
		PyErr_Format(PyExc_RuntimeError, "freq has a different channel count than nChan");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		Py_XDECREF(windowData);
		return NULL;
	}
	
	if(fq->dimensions[0] != times->dimensions[1]) {
		PyErr_Format(PyExc_TypeError, "freq and delays have different channel counts");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		Py_XDECREF(windowData);
		return NULL;
	}
	
	// Get the properties of the data
	nStand = (long) data->dimensions[0];
	nSamps = (long) data->dimensions[1];
	
	// Compute the interger sample offset and the fractional sample delay for each stand
	npy_intp *tLoc;
	long start[nStand];
	long startMax = 0;
	double frac[nStand][nChan];
	tLoc = PyDimMem_NEW(2);
	for(i=0; i<nStand; i++) {
		tLoc[0] = (npy_intp) i;
		tLoc[1] = (npy_intp) (nChan / 2);
		start[i] = (long) round(*(double *) PyArray_GetPtr(times, tLoc) * SampleRate);
		if(start[i] > startMax) {
			startMax = start[i];
		}

		for(j=0; j<nChan; j++) {
			tLoc[1] = (npy_intp) j;
			frac[i][j] = *(double *) PyArray_GetPtr(times, tLoc) - (double) start[i]/SampleRate;
		}
	}
	PyDimMem_FREE(tLoc);

	// Compute the filterbank window for the correct numer of taps
	npy_intp *qLoc;
	double fbWindow[2*nChan*nTaps];
	qLoc = PyDimMem_NEW(1);
	for(i=0; i<2*nChan*nTaps; i++) {
		qLoc[0] = (npy_intp) i;
		fbWindow[i] = sinc((double) (i - nTaps*nChan + 0.5)/2/nChan);
		fbWindow[i] *= *(double *) PyArray_GetPtr(windowData, qLoc);
	}
	PyDimMem_FREE(qLoc);

	// Find out how large the output array needs to be and initialize it
	nFFT = (nSamps - startMax) / 2 / nChan * Overlap - Overlap + 1;
	nFFT = (nFFT / nTaps) * nTaps;
	npy_intp dims[3];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) (nChan - 1);
	dims[2] = (npy_intp) (nFFT / nTaps);
	dataF = (PyArrayObject*) PyArray_SimpleNew(3, dims, PyArray_CDOUBLE);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		Py_XDECREF(windowData);
		return NULL;
	}
	
	// Create the FFTW plan                          
	fftw_complex *inP, *in;                          
	inP = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * 2*nChan);
	fftw_plan p;
	p = fftw_plan_dft_1d(2*nChan, inP, inP, FFTW_FORWARD, FFTW_MEASURE);
	
	// Integer delay, FFT, and fractional delay
	long secStart, fftIndex;
	npy_intp *dLoc, *fLoc;
	
	#ifdef _OPENMP
		#ifdef _MKL
			fftw3_mkl.number_of_user_threads = omp_get_num_threads();
		#endif
		#pragma omp parallel default(shared) private(dLoc, fLoc, qLoc, in, secStart, j, k, fftIndex)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			dLoc = PyDimMem_NEW(2);
			fLoc = PyDimMem_NEW(3);
			qLoc = PyDimMem_NEW(1);
			
			in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * 2*nChan);
			
			dLoc[0] = (npy_intp) i;
			fLoc[0] = (npy_intp) i;
			
			for(j=0; j<nFFT; j++) {
				secStart = start[i] + ((long) (2*nChan*((float) j)/Overlap));
				
				for(k=0; k<2*nChan; k++) {
					dLoc[1] = (npy_intp) (secStart + k);
					qLoc[0] = (npy_intp) k;
					in[k][0] = *(double *) PyArray_GetPtr(data, dLoc) * fbWindow[2*nChan*(j % nTaps) + k];
					in[k][1] = 0.0;
				}	
			
				fftw_execute_dft(p, in, in);
			
				fLoc[2] = (npy_intp) (j / nTaps);
				for(k=0; k<(nChan-1); k++) {
					fLoc[1] = (npy_intp) k;
					qLoc[0] = (npy_intp) k;
					fftIndex = k + 1;
					*(double complex *) PyArray_GetPtr(dataF, fLoc) += (in[fftIndex][0] + imaginary*in[fftIndex][1]) * cexp(-2*imaginary*PI* *(double *) PyArray_GetPtr(fq, qLoc) * frac[i][k]);
				}
			}
			
			PyDimMem_FREE(dLoc);
			PyDimMem_FREE(fLoc);
			PyDimMem_FREE(qLoc);

			fftw_free(in);
		}
	}
	fftw_destroy_plan(p);
	fftw_free(inP);

	Py_XDECREF(data);
	Py_XDECREF(fq);
	Py_XDECREF(times);
	Py_XDECREF(windowData);

	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(PEngineR3_doc, "Perform a series of overlaped filter bank tranforms on real-valued data using OpenMP and windows.");


static PyObject *PEngineC2(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *freq, *delays, *signalsF;
	PyArrayObject *data, *fq, *times, *dataF;
	int nChan = 64;
	int Overlap = 1;
	double SampleRate = 1.0e5;

	long i, j, k, m, nStand, nSamps, nFFT;

	static char *kwlist[] = {"signals", "freq", "delays", "LFFT", "Overlap", "SampleRate", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OOO|iid", kwlist, &signals, &freq, &delays, &nChan, &Overlap, &SampleRate)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it useable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, PyArray_CDOUBLE, 2, 2);
	fq = (PyArrayObject *) PyArray_ContiguousFromObject(freq, PyArray_DOUBLE, 1, 1);
	times = (PyArrayObject *) PyArray_ContiguousFromObject(delays, PyArray_DOUBLE, 2, 2);
	
	// Check data dimensions
	if(data->dimensions[0] != times->dimensions[0]) {
		PyErr_Format(PyExc_TypeError, "signals and delays have different stand counts");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		return NULL;
	}
	
	if(nChan != (fq->dimensions[0]+1)) {
		PyErr_Format(PyExc_RuntimeError, "freq has a different channel count than nChan");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		return NULL;
	}
	
	if(fq->dimensions[0] != times->dimensions[1]) {
		PyErr_Format(PyExc_TypeError, "freq and delays have different channel counts");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		return NULL;
	}

	// Get the properties of the data
	nStand = (long) data->dimensions[0];
	nSamps = (long) data->dimensions[1];
	
	// Compute the interger sample offset and the fractional sample delay for each stand
	npy_intp *tLoc;
	long start[nStand];
	long startMax = 0;
	double frac[nStand][nChan];
	tLoc = PyDimMem_NEW(2);
	for(i=0; i<nStand; i++) {
		tLoc[0] = (npy_intp) i;
		tLoc[1] = (npy_intp) (nChan / 2);
		start[i] = (long) round(*(double *) PyArray_GetPtr(times, tLoc) * SampleRate);
		if(start[i] > startMax) {
			startMax = start[i];
		}

		for(j=0; j<nChan; j++) {
			tLoc[1] = (npy_intp) j;
			frac[i][j] = *(double *) PyArray_GetPtr(times, tLoc) - (double) start[i]/SampleRate;
		}
	}
	PyDimMem_FREE(tLoc);

	// Compute the filterbank window for the correct numer of taps
	double fbWindow[nChan*nTaps];
	double complex tempFB[nChan-1];
	for(i=0; i<nChan*nTaps; i++) {
		fbWindow[i] = sinc((double) (i - nTaps*nChan/2 + 0.5)/nChan);
	}

	// Find out how large the output array needs to be and initialize it
	nFFT = (nSamps - startMax) / nChan * Overlap - Overlap + 1;
	nFFT = (nFFT / nTaps) * nTaps;
	npy_intp dims[3];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) (nChan - 1);
	dims[2] = (npy_intp) (nFFT / nTaps);
	dataF = (PyArrayObject*) PyArray_SimpleNew(3, dims, PyArray_CDOUBLE);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
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
	npy_intp *fLoc;
	double complex *a;
	double *f;
	a = (double complex *) data->data;
	f = (double *) fq->data;
	
	#ifdef _OPENMP
		#ifdef _MKL
			fftw3_mkl.number_of_user_threads = omp_get_num_threads();
		#endif
		#pragma omp parallel default(shared) private(fLoc, in, secStart, tempFB, i, j, k, m, fftIndex)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			fLoc = PyDimMem_NEW(3);
			
			in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nChan);
			
			fLoc[0] = (npy_intp) i;
			
			for(j=0; j<nFFT; j+=nTaps) {
				cblas_zdscal((nChan-1), 0.0, tempFB, 1);
				
				for(m=0; m<nTaps; m++) {
					secStart = start[i] + nSamps * i + nChan*j/Overlap;
					
					for(k=0; k<nChan; k++) {
						in[k][0] = creal(*(a + secStart + k)) * fbWindow[nChan*m + k];
						in[k][1] = cimag(*(a + secStart + k)) * fbWindow[nChan*m + k];
					}
				
					fftw_execute_dft(p, in, in);
				
					for(k=0; k<(nChan-1); k++) {
						fftIndex = ((k+1) + nChan/2) % nChan;
						tempFB[k] += in[fftIndex][0] + imaginary * in[fftIndex][1];
					}
				}
			
				fLoc[2] = (npy_intp) (j / nTaps);
				for(k=0; k<(nChan-1); k++) {
					fLoc[1] = (npy_intp) k;
					*(double complex *) PyArray_GetPtr(dataF, fLoc) = tempFB[k] * cexp(-2*imaginary*PI * *(f + k) * frac[i][k]);
				}
			}
			
			PyDimMem_FREE(fLoc);

			fftw_free(in);
		}
	}
	fftw_destroy_plan(p);
	fftw_free(inP);

	Py_XDECREF(data);
	Py_XDECREF(fq);
	Py_XDECREF(times);

	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(PEngineC2_doc, "Perform a series of overlaped filter bank transfomrs on complex-valued data using OpenMP.");


static PyObject *PEngineC3(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *freq, *delays, *signalsF, *window;
	PyArrayObject *data, *fq, *times, *dataF, *windowData;
	int nChan = 64;
	int Overlap = 1;
	double SampleRate = 1.0e5;

	long i, j, k, m, nStand, nSamps, nFFT;

	static char *kwlist[] = {"signals", "freq", "delays", "LFFT", "Overlap", "SampleRate", "window", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OOO|iidO:set_callback", kwlist, &signals, &freq, &delays, &nChan, &Overlap, &SampleRate, &window)) {
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

	// Bring the data into C and make it useable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, PyArray_CDOUBLE, 2, 2);
	fq = (PyArrayObject *) PyArray_ContiguousFromObject(freq, PyArray_DOUBLE, 1, 1);
	times = (PyArrayObject *) PyArray_ContiguousFromObject(delays, PyArray_DOUBLE, 2, 2);
	
	// Calculate the windowing function
	window = Py_BuildValue("(i)", nTaps*nChan);
	window = PyObject_CallObject(windowFunc, window);
	windowData = (PyArrayObject *) PyArray_ContiguousFromObject(window, PyArray_DOUBLE, 1, 1);
	Py_DECREF(window);

	// Check data dimensions
	if(data->dimensions[0] != times->dimensions[0]) {
		PyErr_Format(PyExc_TypeError, "signals and delays have different stand counts");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		Py_XDECREF(windowData);
		return NULL;
	}
	
	if(nChan != (fq->dimensions[0]+1)) {
		PyErr_Format(PyExc_RuntimeError, "freq has a different channel count than nChan");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		Py_XDECREF(windowData);
		return NULL;
	}
	
	if(fq->dimensions[0] != times->dimensions[1]) {
		PyErr_Format(PyExc_TypeError, "freq and delays have different channel counts");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		Py_XDECREF(windowData);
		return NULL;
	}

	// Get the properties of the data
	nStand = (long) data->dimensions[0];
	nSamps = (long) data->dimensions[1];
	
	// Compute the interger sample offset and the fractional sample delay for each stand
	npy_intp *tLoc;
	long start[nStand];
	long startMax = 0;
	double frac[nStand][nChan];
	tLoc = PyDimMem_NEW(2);
	for(i=0; i<nStand; i++) {
		tLoc[0] = (npy_intp) i;
		tLoc[1] = (npy_intp) (nChan / 2);
		start[i] = (long) round(*(double *) PyArray_GetPtr(times, tLoc) * SampleRate);
		if(start[i] > startMax) {
			startMax = start[i];
		}

		for(j=0; j<nChan; j++) {
			tLoc[1] = (npy_intp) j;
			frac[i][j] = *(double *) PyArray_GetPtr(times, tLoc) - (double) start[i]/SampleRate;
		}
	}
	PyDimMem_FREE(tLoc);

	// Compute the filterbank window for the correct numer of taps
	double fbWindow[nChan*nTaps];
	double complex tempFB[nChan-1];
	double *c;
	c = (double *) windowData->data;
	for(i=0; i<nChan*nTaps; i++) {
		fbWindow[i] = sinc((double) (i - nTaps*nChan/2 + 0.5)/nChan);
		fbWindow[i] *= *(c + i);
	}

	// Find out how large the output array needs to be and initialize it
	nFFT = (nSamps - startMax) / nChan * Overlap - Overlap + 1;
	nFFT = (nFFT / nTaps) * nTaps;
	npy_intp dims[3];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) (nChan - 1);
	dims[2] = (npy_intp) (nFFT / nTaps);
	dataF = (PyArrayObject*) PyArray_SimpleNew(3, dims, PyArray_CDOUBLE);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
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
	npy_intp *fLoc;
	double complex *a;
	double *f;
	a = (double complex *) data->data;
	f = (double *) fq->data;
	
	#ifdef _OPENMP
		#ifdef _MKL
			fftw3_mkl.number_of_user_threads = omp_get_num_threads();
		#endif
		#pragma omp parallel default(shared) private(fLoc, in, secStart, tempFB, i, j, k, m, fftIndex)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			fLoc = PyDimMem_NEW(3);
			
			in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nChan);
			
			fLoc[0] = (npy_intp) i;
			
			for(j=0; j<nFFT; j+=nTaps) {
				cblas_zdscal((nChan-1), 0.0, tempFB, 1);
				
				for(m=0; m<nTaps; m++) {
					secStart = start[i] + nSamps*i + nChan*j/Overlap;
					
					for(k=0; k<nChan; k++) {
						in[k][0] = creal(*(a + secStart + k)) * fbWindow[nChan*m + k];
						in[k][1] = cimag(*(a + secStart + k)) * fbWindow[nChan*m + k];
					}
				
					fftw_execute_dft(p, in, in);
				
					for(k=0; k<(nChan-1); k++) {
						fftIndex = ((k+1) + nChan/2) % nChan;
						tempFB[k] += in[fftIndex][0] + imaginary * in[fftIndex][1];
					}
				}
			
				fLoc[2] = (npy_intp) (j / nTaps);
				for(k=0; k<(nChan-1); k++) {
					fLoc[1] = (npy_intp) k;
					*(double complex *) PyArray_GetPtr(dataF, fLoc) = tempFB[k] * cexp(-2*imaginary*PI * *(f + k) * frac[i][k]);
				}
			}
			
			PyDimMem_FREE(fLoc);

			fftw_free(in);
		}
	}
	fftw_destroy_plan(p);
	fftw_free(inP);

	Py_XDECREF(data);
	Py_XDECREF(fq);
	Py_XDECREF(times);
	Py_XDECREF(windowData);

	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(PEngineC3_doc, "Perform a series of overlaped filter bank transfomrs on complex-valued data using OpenMP and a window.");


/*
  Module Setup - Function Definitions and Documentation
*/

static PyMethodDef CorrelatorMethods[] = {
	{"FEngineR2", FEngineR2, METH_KEYWORDS, FEngineR2_doc}, 
	{"FEngineR3", FEngineR3, METH_KEYWORDS, FEngineR3_doc}, 
	{"FEngineC2", FEngineC2, METH_KEYWORDS, FEngineC2_doc}, 
	{"FEngineC3", FEngineC3, METH_KEYWORDS, FEngineC3_doc}, 
	{"XEngine",   XEngine,   METH_VARARGS,  XEngine_doc}, 
	{"XEngine2",  XEngine2,  METH_VARARGS,  XEngine2_doc}, 
	{"PEngineR2", PEngineR2, METH_KEYWORDS, PEngineR2_doc}, 
	{"PEngineR3", PEngineR3, METH_KEYWORDS, PEngineR3_doc}, 
	{"PEngineC2", PEngineC2, METH_KEYWORDS, PEngineC2_doc}, 
	{"PEngineC3", PEngineC3, METH_KEYWORDS, PEngineC3_doc}, 
	{NULL, NULL, 0, NULL}
};

PyDoc_STRVAR(correlator_doc, \
"C-based F and X engines for the LWA software FX correlator.  These function\n\
are meant to provide an alternative to the lsl.correlator.fx.correlate function and \n\
provide a much-needed speed boost to cross-correlation.\n\
\n\
The function defined in this module are:\n\
  * FEngineR2 -F-engine for computing a series of overlapped Fourier transforms with\n\
    delay corrections for a real-valued (TBW) signal from a collection of stands all at\n\
    once\n\
  * FEngineR3 - Similar to FEngineR2, but allows for a window function to be applied\n\
    to the data.  The window function needs to be evaluated for the correct FFT length\n\
    before being passed to FEngineR3\n\
  * FEngineC2 - F-engine for computing a series of overlapped Fourier transforms with\n\
    delay corrections for a complex-valued (TBN) signal from a collection of stands all at\n\
    once\n\
  * FEngineC3 - Similar to FEngineC2, but allows for a window function to be applied\n\
    to the data.  The window function needs to be evaluated for the correct FFT length\n\
    before being passed to FEngineC3\n\
  * XEngine - Cross-multipy and accumulate cross-power spectra using the inputs for\n\
    two stands\n\
  * XEngine2 - Similar to XEngine, but works with a collection of stands all at\n\
    once\n\
  * PEngineR2 - F-engined based on a 4-tap uniform DFT filter bank for computing a series\n\
    of overlapped transforms with delay corrections for a real-valued (TBW) signals from\n\
    a collection of stands\n\
  * PEngineR3 - Similar to PEngineR2, but allows for a window function to be applied to\n\
    the data.  The window function needs to be evaluated to the correct FFT length and\n\
    number of taps before being passed to PEngineR3\n\
  * PEngineC2 - F-engined based on a 4-tap uniform DFT filter bank for computing a series\n\
    of overlapped transforms with delay corrections for a complex-valued (TBN) signals \n\
    from a collection of stands\n\
  * PEngineC3 - Similar to PEngineC3, but allows for a window function to be applied to\n\
    the data.  The window function needs to be evaluated to the correct FFT length and\n\
    number of taps before being passed to PEngineC3\n\
See the inidividual functions for more details.\n\
\n\
.. note::\n\
	The 'P' engines do not preserve phase and do not yield valid results for use in\n\
	cross-correlation.  Currently only the 'F' engines are know to be useful.\n\
");


/*
  Module Setup - Initialization
*/

PyMODINIT_FUNC init_core(void) {
	(void) Py_InitModule3("_core", CorrelatorMethods, correlator_doc);
	import_array();
}


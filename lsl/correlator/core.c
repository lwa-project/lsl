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
  FFT Functions ("F-engines")
    1. FEngineR2 - FFT a real-valued collection of signals
    2. FEngineR3 - window the data and FFT a real-valued collection of signals
    3. FEngineC2 - FFT a complex-valued collection of signals
    4. FEngineC3 - window the data and FFT a complex-valued collection of signals
*/


static PyObject *FEngineR2(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *freq, *delays, *signalsF;
	PyArrayObject *data, *fq, *times, *dataF, *validF;
	int nChan = 64;
	int Overlap = 1;
	int Clip = 0;
	double SampleRate = 196.0e6;

	long i, j, k, nStand, nSamps, nFFT;
	
	static char *kwlist[] = {"signals", "freq", "delays", "LFFT", "Overlap", "SampleRate", "ClipLevel", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OOO|iidi", kwlist, &signals, &freq, &delays, &nChan, &Overlap, &SampleRate, &Clip)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it usable
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
	
	if(nChan != fq->dimensions[0]) {
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
	
	// Compute the integer sample offset and the fractional sample delay for each stand
	npy_intp *tLoc;
	long start[nStand];
	long startMax = 0;
	double *frac;
	tLoc = PyDimMem_NEW(2);
	frac = (double*) malloc(nStand*nChan * sizeof(double));
	if(frac == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create fractional delay array");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		return NULL;
	}
	
	for(i=0; i<nStand; i++) {
		tLoc[0] = (npy_intp) i;
		tLoc[1] = (npy_intp) (nChan / 2);
		start[i] = (long) round(*(double *) PyArray_GetPtr(times, tLoc) * SampleRate);
		if(start[i] > startMax) {
			startMax = start[i];
		}

		for(j=0; j<nChan; j++) {
			tLoc[1] = (npy_intp) j;
			*(frac + nChan*i + j) = *(double *) PyArray_GetPtr(times, tLoc) - (double) start[i]/SampleRate;
		}
	}
	PyDimMem_FREE(tLoc);

	// Find out how large the output array needs to be and initialize it
	nFFT = (nSamps - startMax) / ((2*nChan)/Overlap) - (2*nChan)/((2*nChan)/Overlap) + 1;
	npy_intp dims[3];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChan;
	dims[2] = (npy_intp) nFFT;
	dataF = (PyArrayObject*) PyArray_SimpleNew(3, dims, NPY_COMPLEX64);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		free(frac);
		return NULL;
	}
	
	// Create an array to store whether or not the FFT window is valid (1) or not (0)
	npy_intp dimsV[2];
	dimsV[0] = (npy_intp) nStand;
	dimsV[1] = (npy_intp) nFFT;
	validF = (PyArrayObject*) PyArray_SimpleNew(2, dimsV, NPY_UINT8);
	if(validF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create valid index array");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		Py_XDECREF(dataF);
		free(frac);
		return NULL;
	}
	
	// Create the FFTW plan                          
	fftwf_complex *inP, *in;                          
	inP = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * 2*nChan);
	fftwf_plan p;
	p = fftwf_plan_dft_1d(2*nChan, inP, inP, FFTW_FORWARD, FFTW_ESTIMATE);
	
	// Integer delay, FFT, and fractional delay
	long secStart, fftIndex;
	npy_intp *dLoc, *fLoc, *qLoc, *vLoc;
	
	// Time-domain blanking control
	double cleanFactor;
	
	#ifdef _OPENMP
		#ifdef _MKL
			fftw3_mkl.number_of_user_threads = omp_get_num_threads();
		#endif
		#pragma omp parallel default(shared) private(dLoc, fLoc, qLoc, vLoc, in, secStart, j, k, fftIndex, cleanFactor)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			dLoc = PyDimMem_NEW(2);
			fLoc = PyDimMem_NEW(3);
			qLoc = PyDimMem_NEW(1);
			vLoc = PyDimMem_NEW(2);
			
			in = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * 2*nChan);
			
			dLoc[0] = (npy_intp) i;
			fLoc[0] = (npy_intp) i;
			vLoc[0] = (npy_intp) i;
			
			for(j=0; j<nFFT; j++) {
				cleanFactor = 1.0;
				secStart = start[i] + ((long) (2*nChan*((float) j)/Overlap));
				
				for(k=0; k<2*nChan; k++) {
					dLoc[1] = (npy_intp) (secStart + k);
					in[k][0] = *(short int *) PyArray_GetPtr(data, dLoc);
					in[k][1] = 0.0;
					
					if( Clip && (in[k][0] >= Clip || in[k][0] <= -Clip) ) {
						cleanFactor = 0.0;
					}
				}	
			
				fftwf_execute_dft(p, in, in);
			
				fLoc[2] = (npy_intp) j;
				vLoc[1] = (npy_intp) j;
				for(k=0; k<nChan; k++) {
					fLoc[1] = (npy_intp) k;
					qLoc[0] = (npy_intp) k;
					fftIndex = k;
					*(float complex *) PyArray_GetPtr(dataF, fLoc) = (cleanFactor*in[fftIndex][0] + imaginary*cleanFactor*in[fftIndex][1]);
					*(float complex *) PyArray_GetPtr(dataF, fLoc) *= cexp(2*imaginary*PI* *(double *) PyArray_GetPtr(fq, qLoc) * *(frac + nChan*i + k));
				}
				
				*(unsigned char *) PyArray_GetPtr(validF, vLoc) = (unsigned char) cleanFactor;
			}
			
			PyDimMem_FREE(dLoc);
			PyDimMem_FREE(fLoc);
			PyDimMem_FREE(qLoc);
			PyDimMem_FREE(vLoc);

			fftwf_free(in);
		}
	}
	fftwf_destroy_plan(p);
	fftwf_free(inP);
	free(frac);

	Py_XDECREF(data);
	Py_XDECREF(fq);
	Py_XDECREF(times);

	signalsF = Py_BuildValue("(OO)", PyArray_Return(dataF), PyArray_Return(validF));
	Py_XDECREF(dataF);
	Py_XDECREF(validF);

	return signalsF;
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
	PyObject *signals, *freq, *delays, *window, *signalsF;
	PyArrayObject *data, *fq, *times, *dataF, *validF, *windowData;
	int nChan = 64;
	int Overlap = 1;
	int Clip = 0;
	double SampleRate = 196.0e6;

	long i, j, k, nStand, nSamps, nFFT;
	
	static char *kwlist[] = {"signals", "freq", "delays", "LFFT", "Overlap", "SampleRate", "ClipLevel", "window", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OOO|iidiO:set_callback", kwlist, &signals, &freq, &delays, &nChan, &Overlap, &SampleRate, &Clip, &window)) {
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
	
	if(nChan != fq->dimensions[0]) {
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
	
	// Compute the integer sample offset and the fractional sample delay for each stand
	npy_intp *tLoc;
	long start[nStand];
	long startMax = 0;
	double *frac;
	tLoc = PyDimMem_NEW(2);
	frac = (double*) malloc(nStand*nChan * sizeof(double));
	if(frac == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create fractional delay array");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		return NULL;
	}
	
	for(i=0; i<nStand; i++) {
		tLoc[0] = (npy_intp) i;
		tLoc[1] = (npy_intp) (nChan / 2);
		start[i] = (long) round(*(double *) PyArray_GetPtr(times, tLoc) * SampleRate);
		if(start[i] > startMax) {
			startMax = start[i];
		}

		for(j=0; j<nChan; j++) {
			tLoc[1] = (npy_intp) j;
			*(frac + nChan*i + j) = *(double *) PyArray_GetPtr(times, tLoc) - (double) start[i]/SampleRate;
		}
	}
	PyDimMem_FREE(tLoc);

	// Find out how large the output array needs to be and initialize it
	nFFT = (nSamps - startMax) / ((2*nChan)/Overlap) - (2*nChan)/((2*nChan)/Overlap) + 1;
	npy_intp dims[3];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChan;
	dims[2] = (npy_intp) nFFT;
	dataF = (PyArrayObject*) PyArray_SimpleNew(3, dims, NPY_COMPLEX64);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		Py_XDECREF(windowData);
		free(frac);
		return NULL;
	}
	
	// Create an array to store whether or not the FFT window is valid (1) or not (0)
	npy_intp dimsV[2];
	dimsV[0] = (npy_intp) nStand;
	dimsV[1] = (npy_intp) nFFT;
	validF = (PyArrayObject*) PyArray_SimpleNew(2, dimsV, NPY_UINT8);
	if(validF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create valid index array");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		Py_XDECREF(dataF);
		Py_XDECREF(windowData);
		free(frac);
		return NULL;
	}
	
	// Create the FFTW plan                          
	fftwf_complex *inP, *in;                          
	inP = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * 2*nChan);
	fftwf_plan p;
	p = fftwf_plan_dft_1d(2*nChan, inP, inP, FFTW_FORWARD, FFTW_ESTIMATE);
	
	// Integer delay, FFT, and fractional delay
	long secStart, fftIndex;
	npy_intp *dLoc, *fLoc, *qLoc, *vLoc;
	
	// Time-domain blanking control
	double cleanFactor;
	
	#ifdef _OPENMP
		#ifdef _MKL
			fftw3_mkl.number_of_user_threads = omp_get_num_threads();
		#endif
		#pragma omp parallel default(shared) private(dLoc, fLoc, qLoc, vLoc, in, secStart, j, k, fftIndex, cleanFactor)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			dLoc = PyDimMem_NEW(2);
			fLoc = PyDimMem_NEW(3);
			qLoc = PyDimMem_NEW(1);
			vLoc = PyDimMem_NEW(2);
			
			in = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * 2*nChan);
			
			dLoc[0] = (npy_intp) i;
			fLoc[0] = (npy_intp) i;
			vLoc[0] = (npy_intp) i;
			
			for(j=0; j<nFFT; j++) {
				cleanFactor = 1.0;
				secStart = start[i] + ((long) (2*nChan*((float) j)/Overlap));
				
				for(k=0; k<2*nChan; k++) {
					dLoc[1] = (npy_intp) (secStart + k);
					qLoc[0] = (npy_intp) k;
					in[k][0] = *(short int *) PyArray_GetPtr(data, dLoc) * *(double *) PyArray_GetPtr(windowData, qLoc);
					in[k][1] = 0.0;
					
					if( Clip && (*(short int *) PyArray_GetPtr(data, dLoc) >= Clip || *(short int *) PyArray_GetPtr(data, dLoc) <= -Clip) ) {
						cleanFactor = 0.0;
					}
				}	
			
				fftwf_execute_dft(p, in, in);
			
				fLoc[2] = (npy_intp) j;
				vLoc[1] = (npy_intp) j;
				for(k=0; k<nChan; k++) {
					fLoc[1] = (npy_intp) k;
					qLoc[0] = (npy_intp) k;
					fftIndex = k;
					*(float complex *) PyArray_GetPtr(dataF, fLoc) = (cleanFactor*in[fftIndex][0] + imaginary*cleanFactor*in[fftIndex][1]);
					*(float complex *) PyArray_GetPtr(dataF, fLoc) *= cexp(2*imaginary*PI* *(double *) PyArray_GetPtr(fq, qLoc) * *(frac + nChan*i + k));
				}
				
				*(unsigned char *) PyArray_GetPtr(validF, vLoc) = (unsigned char) cleanFactor;
			}
			
			PyDimMem_FREE(dLoc);
			PyDimMem_FREE(fLoc);
			PyDimMem_FREE(qLoc);
			PyDimMem_FREE(vLoc);

			fftwf_free(in);
		}
	}
	fftwf_destroy_plan(p);
	fftwf_free(inP);
	free(frac);

	Py_XDECREF(data);
	Py_XDECREF(fq);
	Py_XDECREF(times);
	Py_XDECREF(windowData);

	signalsF = Py_BuildValue("(OO)", PyArray_Return(dataF), PyArray_Return(validF));
	Py_XDECREF(dataF);
	Py_XDECREF(validF);
	
	return signalsF;
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
	PyObject *signals, *freq, *delays, *signalsF;
	PyArrayObject *data, *fq, *times, *dataF, *validF;
	int nChan = 64;
	int Overlap = 1;
	int Clip = 0;
	double SampleRate = 1.0e5;

	long i, j, k, nStand, nSamps, nFFT;

	static char *kwlist[] = {"signals", "freq", "delays", "LFFT", "Overlap", "SampleRate", "ClipLevel", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OOO|iidi", kwlist, &signals, &freq, &delays, &nChan, &Overlap, &SampleRate, &Clip)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it usable
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
	
	if(nChan != fq->dimensions[0]) {
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
	
	// Compute the integer sample offset and the fractional sample delay for each stand
	npy_intp *tLoc;
	long start[nStand];
	long startMax = 0;
	double *frac;
	tLoc = PyDimMem_NEW(2);
	frac = (double*) malloc(nStand*nChan * sizeof(double));
	if(frac == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create fractional delay array");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		return NULL;
	}
	
	for(i=0; i<nStand; i++) {
		tLoc[0] = (npy_intp) i;
		tLoc[1] = (npy_intp) (nChan / 2);
		start[i] = (long) round(*(double *) PyArray_GetPtr(times, tLoc) * SampleRate);
		if(start[i] > startMax) {
			startMax = start[i];
		}

		for(j=0; j<nChan; j++) {
			tLoc[1] = (npy_intp) j;
			*(frac + nChan*i + j) = *(double *) PyArray_GetPtr(times, tLoc) - (double) start[i]/SampleRate;
		}
	}
	PyDimMem_FREE(tLoc);

	// Find out how large the output array needs to be and initialize it
	nFFT = (nSamps - startMax) / (nChan/Overlap) - nChan/(nChan/Overlap) + 1;
	npy_intp dims[3];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChan;
	dims[2] = (npy_intp) nFFT;
	dataF = (PyArrayObject*) PyArray_SimpleNew(3, dims, NPY_COMPLEX64);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		free(frac);
		return NULL;
	}
	
	// Create an array to store whether or not the FFT window is valid (1) or not (0)
	npy_intp dimsV[2];
	dimsV[0] = (npy_intp) nStand;
	dimsV[1] = (npy_intp) nFFT;
	validF = (PyArrayObject*) PyArray_SimpleNew(2, dimsV, NPY_UINT8);
	if(validF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create valid index array");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		Py_XDECREF(dataF);
		free(frac);
		return NULL;
	}

	// Create the FFTW plan
	fftwf_complex *inP, *in;
	inP = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nChan);
	fftwf_plan p;
	p = fftwf_plan_dft_1d(nChan, inP, inP, FFTW_FORWARD, FFTW_ESTIMATE);

	// Integer delay, FFT, and fractional delay
	long secStart, fftIndex;
	npy_intp *dLoc, *fLoc, *qLoc, *vLoc;
	
	// Time-domain blanking control
	double cleanFactor;
	
	#ifdef _OPENMP
		#ifdef _MKL
			fftw3_mkl.number_of_user_threads = omp_get_num_threads();
		#endif
		#pragma omp parallel default(shared) private(dLoc, fLoc, qLoc, vLoc, in, secStart, j, k, fftIndex, cleanFactor)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			dLoc = PyDimMem_NEW(2);
			fLoc = PyDimMem_NEW(3);
			qLoc = PyDimMem_NEW(1);
			vLoc = PyDimMem_NEW(2);
			
			in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nChan);
			
			dLoc[0] = (npy_intp) i;
			fLoc[0] = (npy_intp) i;
			vLoc[0] = (npy_intp) i;
			
			for(j=0; j<nFFT; j++) {
				cleanFactor = 1.0;
				secStart = start[i] + ((long) (nChan*((float) j)/Overlap));
				
				for(k=0; k<nChan; k++) {
					dLoc[1] = (npy_intp) (secStart + k);
					in[k][0] = creal(*(float complex *) PyArray_GetPtr(data, dLoc));
					in[k][1] = cimag(*(float complex *) PyArray_GetPtr(data, dLoc));
					
					if( Clip && cabs(*(float complex *) PyArray_GetPtr(data, dLoc)) >= Clip ) {
						cleanFactor = 0.0;
					}
				}	
			
				fftwf_execute_dft(p, in, in);
			
				fLoc[2] = (npy_intp) j;
				vLoc[1] = (npy_intp) j;
				for(k=0; k<nChan; k++) {
					fLoc[1] = (npy_intp) k;
					qLoc[0] = (npy_intp) k;
					fftIndex = (k + nChan/2) % nChan;
					*(float complex *) PyArray_GetPtr(dataF, fLoc) = (cleanFactor*in[fftIndex][0] + imaginary*cleanFactor*in[fftIndex][1]);
					*(float complex *) PyArray_GetPtr(dataF, fLoc) *= cexp(2*imaginary*PI* *(double *) PyArray_GetPtr(fq, qLoc) * *(frac + nChan*i + k));
				}
				
				*(unsigned char *) PyArray_GetPtr(validF, vLoc) = (unsigned char) cleanFactor;
			}
			
			PyDimMem_FREE(dLoc);
			PyDimMem_FREE(fLoc);
			PyDimMem_FREE(qLoc);
			PyDimMem_FREE(vLoc);

			fftwf_free(in);
		}
	}
	fftwf_destroy_plan(p);
	fftwf_free(inP);
	free(frac);

	Py_XDECREF(data);
	Py_XDECREF(fq);
	Py_XDECREF(times);

	signalsF = Py_BuildValue("(OO)", PyArray_Return(dataF), PyArray_Return(validF));
	Py_XDECREF(dataF);
	Py_XDECREF(validF);
	
	return signalsF;
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
	PyObject *signals, *freq, *delays, *window, *signalsF;
	PyArrayObject *data, *fq, *times, *dataF, *validF, *windowData;
	int nChan = 64;
	int Overlap = 1;
	int Clip = 0;
	double SampleRate = 1.0e5;

	long i, j, k, nStand, nSamps, nFFT;

	static char *kwlist[] = {"signals", "freq", "delays", "LFFT", "Overlap", "SampleRate", "ClipLevel", "window", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OOO|iidiO:set_callback", kwlist, &signals, &freq, &delays, &nChan, &Overlap, &SampleRate, &Clip, &window)) {
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
	
	if(nChan != fq->dimensions[0]) {
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
	
	// Compute the integer sample offset and the fractional sample delay for each stand
	npy_intp *tLoc;
	long start[nStand];
	long startMax = 0;
	double *frac;
	tLoc = PyDimMem_NEW(2);
	frac = (double*) malloc(nStand*nChan * sizeof(double));
	if(frac == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create fractional delay array");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		return NULL;
	}
	
	for(i=0; i<nStand; i++) {
		tLoc[0] = (npy_intp) i;
		tLoc[1] = (npy_intp) (nChan / 2);
		start[i] = (long) round(*(double *) PyArray_GetPtr(times, tLoc) * SampleRate);
		if(start[i] > startMax) {
			startMax = start[i];
		}

		for(j=0; j<nChan; j++) {
			tLoc[1] = (npy_intp) j;
			*(frac + nChan*i + j) = *(double *) PyArray_GetPtr(times, tLoc) - (double) start[i]/SampleRate;
		}
	}
	PyDimMem_FREE(tLoc);

	// Find out how large the output array needs to be and initialize it
	nFFT = (nSamps - startMax) / (nChan/Overlap) - nChan/(nChan/Overlap) + 1;
	npy_intp dims[3];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChan;
	dims[2] = (npy_intp) nFFT;
	dataF = (PyArrayObject*) PyArray_SimpleNew(3, dims, NPY_COMPLEX64);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		Py_XDECREF(windowData);
		free(frac);
		return NULL;
	}
	
	// Create an array to store whether or not the FFT window is valid (1) or not (0)
	npy_intp dimsV[2];
	dimsV[0] = (npy_intp) nStand;
	dimsV[1] = (npy_intp) nFFT;
	validF = (PyArrayObject*) PyArray_SimpleNew(2, dimsV, NPY_UINT8);
	if(validF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create valid index array");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		Py_XDECREF(windowData);
		Py_XDECREF(dataF);
		free(frac);
		return NULL;
	}

	// Create the FFTW plan
	fftwf_complex *inP, *in;
	inP = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nChan);
	fftwf_plan p;
	p = fftwf_plan_dft_1d(nChan, inP, inP, FFTW_FORWARD, FFTW_ESTIMATE);

	// Integer delay, FFT, and fractional delay
	long secStart, fftIndex;
	npy_intp *dLoc, *fLoc, *qLoc, *vLoc;;
	
	// Time-domain blanking control
	double cleanFactor;
	
	#ifdef _OPENMP
		#ifdef _MKL
			fftw3_mkl.number_of_user_threads = omp_get_num_threads();
		#endif
		#pragma omp parallel default(shared) private(dLoc, fLoc, qLoc, vLoc, in, secStart, j, k, fftIndex, cleanFactor)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			dLoc = PyDimMem_NEW(2);
			fLoc = PyDimMem_NEW(3);
			qLoc = PyDimMem_NEW(1);
			vLoc = PyDimMem_NEW(2);
			
			in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nChan);
			
			dLoc[0] = (npy_intp) i;
			fLoc[0] = (npy_intp) i;
			vLoc[0] = (npy_intp) i;
			
			for(j=0; j<nFFT; j++) {
				cleanFactor = 1.0;
				secStart = start[i] + ((long) (nChan*((float) j)/Overlap));
				
				for(k=0; k<nChan; k++) {
					dLoc[1] = (npy_intp) (secStart + k);
					qLoc[0] = (npy_intp) k;
					in[k][0] = creal(*(float complex *) PyArray_GetPtr(data, dLoc) * *(double *) PyArray_GetPtr(windowData, qLoc));
					in[k][1] = cimag(*(float complex *) PyArray_GetPtr(data, dLoc) * *(double *) PyArray_GetPtr(windowData, qLoc));
					
					if( Clip && cabs(*(float complex *) PyArray_GetPtr(data, dLoc)) >= Clip ) {
						cleanFactor = 0.0;
					}
				}	
			
				fftwf_execute_dft(p, in, in);
			
				fLoc[2] = (npy_intp) j;
				vLoc[1] = (npy_intp) j;
				for(k=0; k<nChan; k++) {
					fLoc[1] = (npy_intp) k;
					qLoc[0] = (npy_intp) k;
					fftIndex = (k + nChan/2) % nChan;
					*(float complex *) PyArray_GetPtr(dataF, fLoc) = (cleanFactor*in[fftIndex][0] + imaginary*cleanFactor*in[fftIndex][1]);
					*(float complex *) PyArray_GetPtr(dataF, fLoc) *= cexp(2*imaginary*PI* *(double *) PyArray_GetPtr(fq, qLoc) * *(frac + nChan*i + k));
				}
				
				*(unsigned char *) PyArray_GetPtr(validF, vLoc) = (unsigned char) cleanFactor;
			}
			
			PyDimMem_FREE(dLoc);
			PyDimMem_FREE(fLoc);
			PyDimMem_FREE(qLoc);
			PyDimMem_FREE(vLoc);

			fftwf_free(in);
		}
	}
	fftwf_destroy_plan(p);
	fftwf_free(inP);
	free(frac);

	Py_XDECREF(data);
	Py_XDECREF(fq);
	Py_XDECREF(times);
	Py_XDECREF(windowData);

	signalsF = Py_BuildValue("(OO)", PyArray_Return(dataF), PyArray_Return(validF));
	Py_XDECREF(dataF);
	Py_XDECREF(validF);
	
	return signalsF;
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
	PyArrayObject *data1, *data2, *valid1, *valid2, *vis;
	long nStand, nChan, nFFT, nBL;	

	if(!PyArg_ParseTuple(args, "OOOO", &signals1, &signals2, &sigValid1, &sigValid2)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it usable
	data1 = (PyArrayObject *) PyArray_ContiguousFromObject(signals1, NPY_COMPLEX64, 3, 3);
	data2 = (PyArrayObject *) PyArray_ContiguousFromObject(signals2, NPY_COMPLEX64, 3, 3);
	valid1 = (PyArrayObject *) PyArray_ContiguousFromObject(sigValid1, NPY_UINT8, 2, 2);
	valid2 = (PyArrayObject *) PyArray_ContiguousFromObject(sigValid2, NPY_UINT8, 2, 2);

	// Get channel count and number of FFTs stored
	nStand = (long) data1->dimensions[0];
	nChan = (long) data1->dimensions[1];
	nFFT = (long) data1->dimensions[2];
	nBL = (nStand+1)*nStand/2;
	
	// Create the output visibility array and fill with zeros
	npy_intp dims[2];
	dims[0] = (npy_intp) nBL;
	dims[1] = (npy_intp) nChan;
	vis = (PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_COMPLEX64);
	if(vis == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data1);
		Py_XDECREF(data2);
		Py_XDECREF(valid1);
		Py_XDECREF(valid2);
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
	float complex tempVis;
	float complex *a, *b, *v;
	a = (float complex *) data1->data;
	b = (float complex *) data2->data;
	v = (float complex *) vis->data;
	
	// Time-domain blanking control
	long nActVis;
	unsigned char *u1, *u2;
	u1 = (unsigned char *) valid1->data;
	u2 = (unsigned char *) valid2->data;
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(c, f, nActVis, tempVis)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
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
	Py_XDECREF(data1);
	Py_XDECREF(data2);
	Py_XDECREF(valid1);
	Py_XDECREF(valid2);

	output = Py_BuildValue("O", PyArray_Return(vis));
	Py_XDECREF(vis);

	return output;
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


/*
  Module Setup - Function Definitions and Documentation
*/

static PyMethodDef CorrelatorMethods[] = {
	{"FEngineR2", (PyCFunction) FEngineR2, METH_VARARGS|METH_KEYWORDS, FEngineR2_doc}, 
	{"FEngineR3", (PyCFunction) FEngineR3, METH_VARARGS|METH_KEYWORDS, FEngineR3_doc}, 
	{"FEngineC2", (PyCFunction) FEngineC2, METH_VARARGS|METH_KEYWORDS, FEngineC2_doc}, 
	{"FEngineC3", (PyCFunction) FEngineC3, METH_VARARGS|METH_KEYWORDS, FEngineC3_doc}, 
	{"XEngine2",  (PyCFunction) XEngine2,  METH_VARARGS,               XEngine2_doc }, 
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
	PyModule_AddObject(m, "__version__", PyString_FromString("0.4"));
	PyModule_AddObject(m, "__revision__", PyString_FromString("$Rev$"));
	
	// LSL FFTW Wisdom
	pModule = PyImport_ImportModule("lsl.common.paths");
	pDataPath = PyObject_GetAttrString(pModule, "data");
	sprintf(filename, "%s/fftwf_wisdom.txt", PyString_AsString(pDataPath));
	read_wisdom(filename, m);
}


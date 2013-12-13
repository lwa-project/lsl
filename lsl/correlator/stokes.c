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


static PyObject *FPSDR2(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signalsX, *signalsY, *signalsF;
	PyArrayObject *dataX, *dataY, *dataF;
	int nChan = 64;
	int Overlap = 1;
	int Clip = 0;

	long i, j, k, nStand, nSamps, nFFT;
	
	static char *kwlist[] = {"signalsX", "signalsY", "LFFT", "Overlap", "ClipLevel", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OO|iii", kwlist, &signalsX, &signalsY, &nChan, &Overlap, &Clip)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it usable
	dataX = (PyArrayObject *) PyArray_ContiguousFromObject(signalsX, NPY_INT16, 2, 2);
	dataY = (PyArrayObject *) PyArray_ContiguousFromObject(signalsY, NPY_INT16, 2, 2);
	
	// Get the properties of the data
	nStand = (long) dataX->dimensions[0];
	nSamps = (long) dataX->dimensions[1];
	
	// Make sure the dimensions of X and Y agree
	if( dataY->dimensions[0] != nStand ) {
		PyErr_Format(PyExc_TypeError, "X and Y signals have different stand counts");
		Py_XDECREF(dataX);
		Py_XDECREF(dataY);
		return NULL;
	}
	if( dataY->dimensions[1] != nSamps ) {
		PyErr_Format(PyExc_TypeError, "X and Y signals have different sample counts");
		Py_XDECREF(dataX);
		Py_XDECREF(dataY);
		return NULL;
	}
	
	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / 2 / nChan * Overlap - Overlap + 1;
	npy_intp dims[3];
	dims[0] = (npy_intp) 4;
	dims[1] = (npy_intp) nStand;
	dims[2] = (npy_intp) nChan;
	dataF = (PyArrayObject*) PyArray_SimpleNew(3, dims, NPY_DOUBLE);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(dataX);
		Py_XDECREF(dataY);
		return NULL;
	}
	PyArray_FILLWBYTE(dataF, 0);
	
	// Create the FFTW plan                          
	fftwf_complex *inPX, *inPY, *inX, *inY;                          
	inPX = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * 2*nChan);
	inPY = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * 2*nChan);
	fftwf_plan pX;
	fftwf_plan pY;
	pX = fftwf_plan_dft_1d(2*nChan, inPX, inPX, FFTW_FORWARD, FFTW_ESTIMATE);
	pY = fftwf_plan_dft_1d(2*nChan, inPY, inPY, FFTW_FORWARD, FFTW_ESTIMATE);
	
	// Integer delay, FFT, and fractional delay
	long secStart, fftIndex;
	short int *aX, *aY;
	double *b;
	aX = (short int *) dataX->data;
	aY = (short int *) dataY->data;
	b = (double *) dataF->data;
	
	// Time-domain blanking control
	double cleanFactor;
	long nActFFT[nStand];
	
	#ifdef _OPENMP
		#ifdef _MKL
			fftw3_mkl.number_of_user_threads = omp_get_num_threads();
		#endif
		#pragma omp parallel default(shared) private(inX, inY, secStart, i, j, k, fftIndex, cleanFactor)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			nActFFT[i] = 0;
			inX = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * 2*nChan);
			inY = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * 2*nChan);
			
			for(j=0; j<nFFT; j++) {
				cleanFactor = 1.0;
				secStart = (long) (nSamps * i) + 2*nChan*j/Overlap;
				
				for(k=0; k<2*nChan; k++) {
					inX[k][0] = (double) *(aX + secStart + k);
					inX[k][1] = 0.0;
					
					if( Clip && (inX[k][0] >= Clip || inX[k][0] <= -Clip) ) {
						cleanFactor = 0.0;
					}
					
					inY[k][0] = (double) *(aY + secStart + k);
					inY[k][1] = 0.0;
					
					if( Clip && (inY[k][0] >= Clip || inY[k][0] <= -Clip) ) {
						cleanFactor = 0.0;
					}
				}
				
				fftwf_execute_dft(pX, inX, inX);
				fftwf_execute_dft(pY, inY, inY);
				
				for(k=0; k<nChan; k++) {
					fftIndex = k;
					
					// I
					*(b + 0*nChan*nStand + nChan*i + k) += cleanFactor*inX[fftIndex][0]*inX[fftIndex][0];
					*(b + 0*nChan*nStand + nChan*i + k) += cleanFactor*inX[fftIndex][1]*inX[fftIndex][1];
					*(b + 0*nChan*nStand + nChan*i + k) += cleanFactor*inY[fftIndex][0]*inY[fftIndex][0];
					*(b + 0*nChan*nStand + nChan*i + k) += cleanFactor*inY[fftIndex][1]*inY[fftIndex][1];
					
					// Q
					*(b + 1*nChan*nStand + nChan*i + k) += cleanFactor*inX[fftIndex][0]*inX[fftIndex][0];
					*(b + 1*nChan*nStand + nChan*i + k) += cleanFactor*inX[fftIndex][1]*inX[fftIndex][1];
					*(b + 1*nChan*nStand + nChan*i + k) -= cleanFactor*inY[fftIndex][0]*inY[fftIndex][0];
					*(b + 1*nChan*nStand + nChan*i + k) -= cleanFactor*inY[fftIndex][1]*inY[fftIndex][1];
					
					// U
					*(b + 2*nChan*nStand + nChan*i + k) += 2*cleanFactor*inX[fftIndex][0]*inY[fftIndex][0];
					*(b + 2*nChan*nStand + nChan*i + k) += 2*cleanFactor*inX[fftIndex][1]*inY[fftIndex][1];
					
					// V
					*(b + 3*nChan*nStand + nChan*i + k) += 2*cleanFactor*inX[fftIndex][1]*inY[fftIndex][0];
					*(b + 3*nChan*nStand + nChan*i + k) -= 2*cleanFactor*inX[fftIndex][0]*inY[fftIndex][1];
				}
				
				nActFFT[i] += (long) cleanFactor;
			}
			fftwf_free(inX);
			fftwf_free(inY);
		}
	}
	fftwf_destroy_plan(pX);
	fftwf_destroy_plan(pY);
	fftwf_free(inPX);
	fftwf_free(inPY);
	
	// cblas_dscal(nChan*nStand, 1.0/(2*nChan*nFFT), b, 1);
	for(i=0; i<4; i++) {
		for(j=0; j<nStand; j++) {
			cblas_dscal(nChan, 1.0/(2*nChan*nActFFT[j]), (b + i*nChan*nStand + j*nChan), 1);
		}
	}

	Py_XDECREF(dataX);
	Py_XDECREF(dataY);

	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(FPSDR2_doc, \
"Perform a series of Fourier transforms on real-valued data to get the PSD\n\
for the four Stokes parameters: I, Q, U, and V.\n\
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
 * psd: 3-D numpy.double (Stokes parameter (I,Q,U,V) by stands by channels) of PSD data\n\
");


static PyObject *FPSDR3(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signalsX, *signalsY, *signalsF, *window;
	PyArrayObject *dataX, *dataY, *dataF, *windowData;
	int nChan = 64;
	int Overlap = 1;
	int Clip = 0;
	
	long i, j, k, nStand, nSamps, nFFT;
	
	static char *kwlist[] = {"signalsX", "signalsY", "LFFT", "Overlap", "ClipLevel", "window", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OO|iiiO:set_callback", kwlist, &signalsX, &signalsY, &nChan, &Overlap, &Clip, &window)) {
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
	dataX = (PyArrayObject *) PyArray_ContiguousFromObject(signalsX, NPY_INT16, 2, 2);
	dataY = (PyArrayObject *) PyArray_ContiguousFromObject(signalsY, NPY_INT16, 2, 2);
	
	// Calculate the windowing function
	window = Py_BuildValue("(i)", 2*nChan);
	window = PyObject_CallObject(windowFunc, window);
	windowData = (PyArrayObject *) PyArray_ContiguousFromObject(window, NPY_DOUBLE, 1, 1);
	Py_DECREF(window);
	
	// Get the properties of the data
	nStand = (long) dataX->dimensions[0];
	nSamps = (long) dataX->dimensions[1];
	
	// Make sure the dimensions of X and Y agree
	if( dataY->dimensions[0] != nStand ) {
		PyErr_Format(PyExc_TypeError, "X and Y signals have different stand counts");
		Py_XDECREF(dataX);
		Py_XDECREF(dataY);
		return NULL;
	}
	if( dataY->dimensions[1] != nSamps ) {
		PyErr_Format(PyExc_TypeError, "X and Y signals have different sample counts");
		Py_XDECREF(dataX);
		Py_XDECREF(dataY);
		return NULL;
	}

	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / 2 / nChan * Overlap - Overlap + 1;
	npy_intp dims[3];
	dims[0] = (npy_intp) 4;
	dims[1] = (npy_intp) nStand;
	dims[2] = (npy_intp) nChan;
	dataF = (PyArrayObject*) PyArray_SimpleNew(3, dims, NPY_DOUBLE);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(dataX);
		Py_XDECREF(dataY);
		Py_XDECREF(windowData);
		return NULL;
	}
	PyArray_FILLWBYTE(dataF, 0);
	
	// Create the FFTW plan                          
	fftwf_complex *inPX, *inPY, *inX, *inY;                          
	inPX = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * 2*nChan);
	inPY = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * 2*nChan);
	fftwf_plan pX;
	fftwf_plan pY;
	pX = fftwf_plan_dft_1d(2*nChan, inPX, inPX, FFTW_FORWARD, FFTW_ESTIMATE);
	pY = fftwf_plan_dft_1d(2*nChan, inPY, inPY, FFTW_FORWARD, FFTW_ESTIMATE);
	
	// Integer delay, FFT, and fractional delay
	long secStart, fftIndex;
	short int *aX, *aY;
	double *b, *c;
	aX = (short int *) dataX->data;
	aY = (short int *) dataY->data;
	b = (double *) dataF->data;
	c = (double *) windowData->data;
	
	// Time-domain blanking control
	double cleanFactor;
	long nActFFT[nStand];
	
	#ifdef _OPENMP
		#ifdef _MKL
			fftw3_mkl.number_of_user_threads = omp_get_num_threads();
		#endif
		#pragma omp parallel default(shared) private(inX, inY, secStart, i, j, k, fftIndex, cleanFactor)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			nActFFT[i] = 0;
			inX = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * 2*nChan);
			inY = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * 2*nChan);
			
			for(j=0; j<nFFT; j++) {
				cleanFactor = 1.0;
				secStart = nSamps * i + 2*nChan*j/Overlap;
				
				for(k=0; k<2*nChan; k++) {
					inX[k][0] = (double) *(aX + secStart + k) * *(c + k);
					inX[k][1] = 0.0;
					
					if( Clip && (inX[k][0] >= Clip || inX[k][0] <= -Clip) ) {
						cleanFactor = 0.0;
					}
					
					inY[k][0] = (double) *(aY + secStart + k) * *(c + k);
					inY[k][1] = 0.0;
					
					if( Clip && (inY[k][0] >= Clip || inY[k][0] <= -Clip) ) {
						cleanFactor = 0.0;
					}
				}
				
				fftwf_execute_dft(pX, inX, inX);
				fftwf_execute_dft(pY, inY, inY);
				
				for(k=0; k<nChan; k++) {
					fftIndex = k;
					
					// I
					*(b + 0*nChan*nStand + nChan*i + k) += cleanFactor*inX[fftIndex][0]*inX[fftIndex][0];
					*(b + 0*nChan*nStand + nChan*i + k) += cleanFactor*inX[fftIndex][1]*inX[fftIndex][1];
					*(b + 0*nChan*nStand + nChan*i + k) += cleanFactor*inY[fftIndex][0]*inY[fftIndex][0];
					*(b + 0*nChan*nStand + nChan*i + k) += cleanFactor*inY[fftIndex][1]*inY[fftIndex][1];
					
					// Q
					*(b + 1*nChan*nStand + nChan*i + k) += cleanFactor*inX[fftIndex][0]*inX[fftIndex][0];
					*(b + 1*nChan*nStand + nChan*i + k) += cleanFactor*inX[fftIndex][1]*inX[fftIndex][1];
					*(b + 1*nChan*nStand + nChan*i + k) -= cleanFactor*inY[fftIndex][0]*inY[fftIndex][0];
					*(b + 1*nChan*nStand + nChan*i + k) -= cleanFactor*inY[fftIndex][1]*inY[fftIndex][1];
					
					// U
					*(b + 2*nChan*nStand + nChan*i + k) += 2*cleanFactor*inX[fftIndex][0]*inY[fftIndex][0];
					*(b + 2*nChan*nStand + nChan*i + k) += 2*cleanFactor*inX[fftIndex][1]*inY[fftIndex][1];
					
					// V
					*(b + 3*nChan*nStand + nChan*i + k) += 2*cleanFactor*inX[fftIndex][1]*inY[fftIndex][0];
					*(b + 3*nChan*nStand + nChan*i + k) -= 2*cleanFactor*inX[fftIndex][0]*inY[fftIndex][1];
				}
				
				nActFFT[i] += (long) cleanFactor;
			}
			fftwf_free(inX);
			fftwf_free(inY);
		}
	}
	fftwf_destroy_plan(pX);
	fftwf_destroy_plan(pY);
	fftwf_free(inPX);
	fftwf_free(inPY);
	
	// cblas_dscal(nChan*nStand, 1.0/(2*nChan*nFFT), b, 1);
	for(i=0; i<4; i++) {
		for(j=0; j<nStand; j++) {
			cblas_dscal(nChan, 1.0/(2*nChan*nActFFT[j]), (b + i*nChan*nStand + j*nChan), 1);
		}
	}

	Py_XDECREF(dataX);
	Py_XDECREF(dataY);
	Py_XDECREF(windowData);

	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(FPSDR3_doc, \
"Perform a series of Fourier transforms with windows on real-valued data to\n\
get the PSD for the four Stokes parameters: I, Q, U, and V.\n\
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
 * psd: 3-D numpy.double (Stokes parameter (I,Q,U,V) by stands by channels) of PSD data\n\
");


static PyObject *FPSDC2(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signalsX, *signalsY, *signalsF;
	PyArrayObject *dataX, *dataY, *dataF;
	int nChan = 64;
	int Overlap = 1;
	int Clip = 0;

	long i, j, k, nStand, nSamps, nFFT;

	static char *kwlist[] = {"signalsX", "signalsY", "LFFT", "Overlap", "ClipLevel", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OO|iii", kwlist, &signalsX, &signalsY, &nChan, &Overlap, &Clip)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it usable
	dataX = (PyArrayObject *) PyArray_ContiguousFromObject(signalsX, NPY_COMPLEX64, 2, 2);
	dataY = (PyArrayObject *) PyArray_ContiguousFromObject(signalsY, NPY_COMPLEX64, 2, 2);

	// Get the properties of the data
	nStand = (long) dataX->dimensions[0];
	nSamps = (long) dataX->dimensions[1];
	
	// Make sure the dimensions of X and Y agree
	if( dataY->dimensions[0] != nStand ) {
		PyErr_Format(PyExc_TypeError, "X and Y signals have different stand counts");
		Py_XDECREF(dataX);
		Py_XDECREF(dataY);
		return NULL;
	}
	if( dataY->dimensions[1] != nSamps ) {
		PyErr_Format(PyExc_TypeError, "X and Y signals have different sample counts");
		Py_XDECREF(dataX);
		Py_XDECREF(dataY);
		return NULL;
	}
	
	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / nChan * Overlap - Overlap + 1;
	npy_intp dims[3];
	dims[0] = (npy_intp) 4;
	dims[1] = (npy_intp) nStand;
	dims[2] = (npy_intp) nChan;
	dataF = (PyArrayObject*) PyArray_SimpleNew(3, dims, NPY_DOUBLE);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(dataX);
		Py_XDECREF(dataY);
		return NULL;
	}
	PyArray_FILLWBYTE(dataF, 0);

	// Create the FFTW plan
	fftwf_complex *inPX, *inPY, *inX, *inY;
	inPX = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nChan);
	inPY = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nChan);
	fftwf_plan pX, pY;
	pX = fftwf_plan_dft_1d(nChan, inPX, inPX, FFTW_FORWARD, FFTW_ESTIMATE);
	pY = fftwf_plan_dft_1d(nChan, inPY, inPY, FFTW_FORWARD, FFTW_ESTIMATE);

	// Integer delay, FFT, and fractional delay
	long secStart;
	float complex *aX, *aY;
	double *b;
	aX = (float complex *) dataX->data;
	aY = (float complex *) dataY->data;
	b = (double *) dataF->data;
	
	// Time-domain blanking control
	double cleanFactor;
	long nActFFT[nStand];
	
	#ifdef _OPENMP
		#ifdef _MKL
			fftw3_mkl.number_of_user_threads = omp_get_num_threads();
		#endif
		#pragma omp parallel default(shared) private(inX, inY, secStart, i, j, k, cleanFactor)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			nActFFT[i] = 0;
			inX = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nChan);
			inY = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nChan);
			
			for(j=0; j<nFFT; j++) {
				cleanFactor = 1.0;
				secStart = nSamps * i + nChan*j/Overlap;
				
				for(k=0; k<nChan; k++) {
					inX[k][0] = creal(*(aX + secStart + k));
					inX[k][1] = cimag(*(aX + secStart + k));
					
					if( Clip && cabs(*(aX + secStart + k)) >= Clip ) {
						cleanFactor = 0.0;
					}
					
					inY[k][0] = creal(*(aY + secStart + k));
					inY[k][1] = cimag(*(aY + secStart + k));
					
					if( Clip && cabs(*(aY + secStart + k)) >= Clip ) {
						cleanFactor = 0.0;
					}
				}
				
				fftwf_execute_dft(pX, inX, inX);
				fftwf_execute_dft(pY, inY, inY);
				
				for(k=0; k<nChan; k++) {
					// I
					*(b + 0*nChan*nStand + nChan*i + k) += cleanFactor*inX[k][0]*inX[k][0];
					*(b + 0*nChan*nStand + nChan*i + k) += cleanFactor*inX[k][1]*inX[k][1];
					*(b + 0*nChan*nStand + nChan*i + k) += cleanFactor*inY[k][0]*inY[k][0];
					*(b + 0*nChan*nStand + nChan*i + k) += cleanFactor*inY[k][1]*inY[k][1];
					
					// Q
					*(b + 1*nChan*nStand + nChan*i + k) += cleanFactor*inX[k][0]*inX[k][0];
					*(b + 1*nChan*nStand + nChan*i + k) += cleanFactor*inX[k][1]*inX[k][1];
					*(b + 1*nChan*nStand + nChan*i + k) -= cleanFactor*inY[k][0]*inY[k][0];
					*(b + 1*nChan*nStand + nChan*i + k) -= cleanFactor*inY[k][1]*inY[k][1];
					
					// U
					*(b + 2*nChan*nStand + nChan*i + k) += 2*cleanFactor*inX[k][0]*inY[k][0];
					*(b + 2*nChan*nStand + nChan*i + k) += 2*cleanFactor*inX[k][1]*inY[k][1];
					
					// V
					*(b + 3*nChan*nStand + nChan*i + k) += 2*cleanFactor*inX[k][1]*inY[k][0];
					*(b + 3*nChan*nStand + nChan*i + k) -= 2*cleanFactor*inX[k][0]*inY[k][1];
				}
				
				nActFFT[i] += (long) cleanFactor;
			}
			fftwf_free(inX);
			fftwf_free(inY);
		}
	}
	fftwf_destroy_plan(pX);
	fftwf_destroy_plan(pY);
	fftwf_free(inPX);
	fftwf_free(inPY);

	// Shift and scale FFTs
	double *temp, *temp2;
	temp2 = (double *) malloc(sizeof(double)*nChan/2);
	for(i=0; i<4; i++) {
		for(j=0; j<nStand; j++) {
			temp = b + nChan*nStand*i + nChan*j;
			memcpy(temp2, temp, sizeof(double)*nChan/2);
			memmove(temp, temp+nChan/2, sizeof(double)*nChan/2);
			memcpy(temp+nChan/2, temp2, sizeof(double)*nChan/2);
			
			cblas_dscal(nChan, 1.0/(nChan*nActFFT[j]), (b + i*nChan*nStand + j*nChan), 1);
		}
	}
	free(temp2);

	Py_XDECREF(dataX);
	Py_XDECREF(dataY);
	
	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(FPSDC2_doc, \
"Perform a series of Fourier transforms on complex-valued data to get the\n\
PSD for the four Stokes parameters: I, Q, U, and V.\n\
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
 * psd: 3-D numpy.double (Stokes parameter (I,Q,U,V) by stands by channels) of PSD data\n\
");


static PyObject *FPSDC3(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signalsX, *signalsY, *signalsF, *window;
	PyArrayObject *dataX, *dataY, *dataF, *windowData;
	int nChan = 64;
	int Overlap = 1;
	int Clip = 0;

	long i, j, k, nStand, nSamps, nFFT;

	static char *kwlist[] = {"signalsX", "signalsY", "LFFT", "Overlap", "ClipLevel", "window", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OO|iiiO:set_callback", kwlist, &signalsX, &signalsY, &nChan, &Overlap, &Clip, &window)) {
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
	dataX = (PyArrayObject *) PyArray_ContiguousFromObject(signalsX, NPY_COMPLEX64, 2, 2);
	dataY = (PyArrayObject *) PyArray_ContiguousFromObject(signalsY, NPY_COMPLEX64, 2, 2);
	
	// Calculate the windowing function
	window = Py_BuildValue("(i)", nChan);
	window = PyObject_CallObject(windowFunc, window);
	windowData = (PyArrayObject *) PyArray_ContiguousFromObject(window, NPY_DOUBLE, 1, 1);
	Py_DECREF(window);

	// Get the properties of the data
	nStand = (long) dataX->dimensions[0];
	nSamps = (long) dataX->dimensions[1];

	// Make sure the dimensions of X and Y agree
	if( dataY->dimensions[0] != nStand ) {
		PyErr_Format(PyExc_TypeError, "X and Y signals have different stand counts");
		Py_XDECREF(dataX);
		Py_XDECREF(dataY);
		return NULL;
	}
	if( dataY->dimensions[1] != nSamps ) {
		PyErr_Format(PyExc_TypeError, "X and Y signals have different sample counts");
		Py_XDECREF(dataX);
		Py_XDECREF(dataY);
		return NULL;
	}
	
	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / nChan * Overlap - Overlap + 1;
	npy_intp dims[3];
	dims[0] = (npy_intp) 4;
	dims[1] = (npy_intp) nStand;
	dims[2] = (npy_intp) nChan;
	dataF = (PyArrayObject*) PyArray_SimpleNew(3, dims, NPY_DOUBLE);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(dataX);
		Py_XDECREF(dataY);
		Py_XDECREF(windowData);
		return NULL;
	}
	PyArray_FILLWBYTE(dataF, 0);

	// Create the FFTW plan
	fftwf_complex *inPX, *inPY, *inX, *inY;
	inPX = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nChan);
	inPY = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nChan);
	fftwf_plan pX, pY;
	pX = fftwf_plan_dft_1d(nChan, inPX, inPX, FFTW_FORWARD, FFTW_ESTIMATE);
	pY = fftwf_plan_dft_1d(nChan, inPY, inPY, FFTW_FORWARD, FFTW_ESTIMATE);

	// Integer delay, FFT, and fractional delay
	long secStart;
	float complex *aX, *aY;
	double *b, *c;
	aX = (float complex *) dataX->data;
	aY = (float complex *) dataY->data;
	b = (double *) dataF->data;
	c = (double *) windowData->data;
	
	// Time-domain blanking control
	double cleanFactor;
	long nActFFT[nStand];
	
	#ifdef _OPENMP
		#ifdef _MKL
			fftw3_mkl.number_of_user_threads = omp_get_num_threads();
		#endif
		#pragma omp parallel default(shared) private(inX, inY, secStart, i, j, k, cleanFactor)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			nActFFT[i] = 0;
			inX = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nChan);
			inY = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nChan);
			
			for(j=0; j<nFFT; j++) {
				cleanFactor = 1.0;
				secStart = nSamps * i + nChan*j/Overlap;
				
				for(k=0; k<nChan; k++) {
					inX[k][0] = creal(*(aX + secStart + k)) * *(c + k);
					inX[k][1] = cimag(*(aX + secStart + k)) * *(c + k);
					
					if( Clip && cabs(*(aX + secStart + k)) >= Clip ) {
						cleanFactor = 0.0;
					}
					
					inY[k][0] = creal(*(aY + secStart + k)) * *(c + k);
					inY[k][1] = cimag(*(aY + secStart + k)) * *(c + k);
					
					if( Clip && cabs(*(aY + secStart + k)) >= Clip ) {
						cleanFactor = 0.0;
					}
				}
				
				fftwf_execute_dft(pX, inX, inX);
				fftwf_execute_dft(pY, inY, inY);
				
				for(k=0; k<nChan; k++) {
					// I
					*(b + 0*nChan*nStand + nChan*i + k) += cleanFactor*inX[k][0]*inX[k][0];
					*(b + 0*nChan*nStand + nChan*i + k) += cleanFactor*inX[k][1]*inX[k][1];
					*(b + 0*nChan*nStand + nChan*i + k) += cleanFactor*inY[k][0]*inY[k][0];
					*(b + 0*nChan*nStand + nChan*i + k) += cleanFactor*inY[k][1]*inY[k][1];
					
					// Q
					*(b + 1*nChan*nStand + nChan*i + k) += cleanFactor*inX[k][0]*inX[k][0];
					*(b + 1*nChan*nStand + nChan*i + k) += cleanFactor*inX[k][1]*inX[k][1];
					*(b + 1*nChan*nStand + nChan*i + k) -= cleanFactor*inY[k][0]*inY[k][0];
					*(b + 1*nChan*nStand + nChan*i + k) -= cleanFactor*inY[k][1]*inY[k][1];
					
					// U
					*(b + 2*nChan*nStand + nChan*i + k) += 2*cleanFactor*inX[k][0]*inY[k][0];
					*(b + 2*nChan*nStand + nChan*i + k) += 2*cleanFactor*inX[k][1]*inY[k][1];
					
					// V
					*(b + 3*nChan*nStand + nChan*i + k) += 2*cleanFactor*inX[k][1]*inY[k][0];
					*(b + 3*nChan*nStand + nChan*i + k) -= 2*cleanFactor*inX[k][0]*inY[k][1];
				}
				
				nActFFT[i] += (long) cleanFactor;
			}
			fftwf_free(inX);
			fftwf_free(inY);
		}
	}
	fftwf_destroy_plan(pX);
	fftwf_destroy_plan(pY);
	fftwf_free(inPX);
	fftwf_free(inPY);

	// Shift and scale FFTs
	double *temp, *temp2;
	temp2 = (double *) malloc(sizeof(double)*nChan/2);
	for(i=0; i<4; i++) {
		for(j=0; j<nStand; j++) {
			temp = b + nChan*nStand*i + nChan*j;
			memcpy(temp2, temp, sizeof(double)*nChan/2);
			memmove(temp, temp+nChan/2, sizeof(double)*nChan/2);
			memcpy(temp+nChan/2, temp2, sizeof(double)*nChan/2);
			
			cblas_dscal(nChan, 1.0/(nChan*nActFFT[j]), (b + i*nChan*nStand + j*nChan), 1);
		}
	}
	free(temp2);
	
	Py_XDECREF(dataX);
	Py_XDECREF(dataY);
	Py_XDECREF(windowData);

	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(FPSDC3_doc, \
"Perform a series of Fourier transforms with windows on complex-valued data\n\
to get the PSD for the four Stokes parameters: I, Q, U, and V.\n\
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
 * psd: 3-D numpy.double (Stokes parameter (I,Q,U,V) by stands by channels) of PSD data\n\
");


/*
  Cross-Multiplication And Accumulation Function ("X Engines")
    1. XEngine2 - XMAC two collections of signals
*/

static PyObject *XEngine2(PyObject *self, PyObject *args) {
	PyObject *signalsX, *signalsY, *sigValidX, *sigValidY, *output;
	PyArrayObject *dataX, *dataY, *validX, *validY, *vis;
	long nStand, nChan, nFFT, nBL;	

	if(!PyArg_ParseTuple(args, "OOOO", &signalsX, &signalsY, &sigValidX, &sigValidY)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it usable
	dataX = (PyArrayObject *) PyArray_ContiguousFromObject(signalsX, NPY_COMPLEX64, 3, 3);
	dataY = (PyArrayObject *) PyArray_ContiguousFromObject(signalsY, NPY_COMPLEX64, 3, 3);
	validX = (PyArrayObject *) PyArray_ContiguousFromObject(sigValidX, NPY_UINT8, 2, 2);
	validY = (PyArrayObject *) PyArray_ContiguousFromObject(sigValidY, NPY_UINT8, 2, 2);

	// Get channel count and number of FFTs stored
	nStand = (long) dataX->dimensions[0];
	nChan = (long) dataX->dimensions[1];
	nFFT = (long) dataX->dimensions[2];
	nBL = (nStand+1)*nStand/2;
	
	// Create the output visibility array and fill with zeros
	npy_intp dims[3];
	dims[0] = (npy_intp) 4;
	dims[1] = (npy_intp) nBL;
	dims[2] = (npy_intp) nChan;
	vis = (PyArrayObject*) PyArray_SimpleNew(3, dims, NPY_COMPLEX64);
	if(vis == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(dataX);
		Py_XDECREF(dataY);
		Py_XDECREF(validX);
		Py_XDECREF(validY);
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
	float complex tempVis1, tempVis2;
	float complex *a, *b, *v;
	a = (float complex *) dataX->data;
	b = (float complex *) dataY->data;
	v = (float complex *) vis->data;
	
	// Time-domain blanking control
	long nActVis;
	unsigned char *u1, *u2;
	u1 = (unsigned char *) validX->data;
	u2 = (unsigned char *) validY->data;
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(c, f, nActVis, tempVis1, tempVis2)
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
				// I
				cblas_cdotc_sub(nFFT, (a + mapper[bl][1]*nChan*nFFT + c*nFFT), 1, (a + mapper[bl][0]*nChan*nFFT + c*nFFT), 1, &tempVis1);
				cblas_cdotc_sub(nFFT, (b + mapper[bl][1]*nChan*nFFT + c*nFFT), 1, (b + mapper[bl][0]*nChan*nFFT + c*nFFT), 1, &tempVis2);
				*(v + 0*nBL*nChan + bl*nChan + c) = (tempVis1 + tempVis2) / nActVis;
				
				// Q
				*(v + 1*nBL*nChan + bl*nChan + c) = (tempVis1 - tempVis2) / nActVis;
				
				// U
				cblas_cdotc_sub(nFFT, (b + mapper[bl][1]*nChan*nFFT + c*nFFT), 1, (a + mapper[bl][0]*nChan*nFFT + c*nFFT), 1, &tempVis1);
				cblas_cdotc_sub(nFFT, (a + mapper[bl][0]*nChan*nFFT + c*nFFT), 1, (b + mapper[bl][1]*nChan*nFFT + c*nFFT), 1, &tempVis2);
				*(v + 2*nBL*nChan + bl*nChan + c) = (tempVis1 + tempVis2) / nActVis;
				
				// V
				*(v + 3*nBL*nChan + bl*nChan + c) = (tempVis1 - tempVis2) / nActVis / _Complex_I;
			}
		}
	}
	Py_XDECREF(dataX);
	Py_XDECREF(dataY);
	Py_XDECREF(validX);
	Py_XDECREF(validY);

	output = Py_BuildValue("O", PyArray_Return(vis));
	Py_XDECREF(vis);

	return output;
}

PyDoc_STRVAR(XEngine2_doc, \
"Perform all XMACs for a data stream out of the F engine using OpenMP that\n\
creates the four Stokes parameters: I, Q, U, and V.\n\
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
  * visibility: 3-D numpy.cdouble (Stokes parameter (I,Q,U,V by baseline by\n\
  channel) array of cross-correlated and averaged visibility data.\n\
");


/*
  Module Setup - Function Definitions and Documentation
*/

static PyMethodDef StokesMethods[] = {
	{"FPSDR2",   (PyCFunction) FPSDR2,    METH_VARARGS|METH_KEYWORDS, FPSDR2_doc   }, 
	{"FPSDR3",   (PyCFunction) FPSDR3,    METH_VARARGS|METH_KEYWORDS, FPSDR3_doc   }, 
	{"FPSDC2",   (PyCFunction) FPSDC2,    METH_VARARGS|METH_KEYWORDS, FPSDC2_doc   }, 
	{"FPSDC3",   (PyCFunction) FPSDC3,    METH_VARARGS|METH_KEYWORDS, FPSDC3_doc   }, 
	{"XEngine2", (PyCFunction) XEngine2,  METH_VARARGS,               XEngine2_doc }, 
	{NULL,       NULL,      0,                          NULL         }
};

PyDoc_STRVAR(stokes_doc, \
"Extension to take X and Y timeseries data and create the four Stokes\n\
parameters.\n\
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
\n\
Also included is an X-Engine for use with the lsl.correlator._core module to\n\
perform cross-correlations for the stokes parameters.\n\
\n\
See the inidividual functions for more details.");


/*
  Module Setup - Initialization
*/

PyMODINIT_FUNC init_stokes(void) {
	char filename[256];
	PyObject *m, *pModule, *pDataPath;

	// Module definitions and functions
	m = Py_InitModule3("_stokes", StokesMethods, stokes_doc);
	import_array();
	
	// Version and revision information
	PyModule_AddObject(m, "__version__", PyString_FromString("0.1"));
	PyModule_AddObject(m, "__revision__", PyString_FromString("$Rev$"));
	
	// LSL FFTW Wisdom
	pModule = PyImport_ImportModule("lsl.common.paths");
	pDataPath = PyObject_GetAttrString(pModule, "data");
	sprintf(filename, "%s/fftwf_wisdom.txt", PyString_AsString(pDataPath));
	read_wisdom(filename, m);
}

#include "Python.h"
#include <math.h>
#include <stdio.h>
#include <complex.h>
#include <cblas.h>
#include <fftw3.h>
#include <stdlib.h>

#ifdef _OPENMP
	#include <omp.h>
	
	// OpenMP scheduling method
	#ifndef OMP_SCHEDULER
	#define OMP_SCHEDULER dynamic
	#endif
#endif

#include "numpy/arrayobject.h"
#include "numpy/npy_math.h"


#if PY_MAJOR_VERSION >= 3
	#define PyCapsule_Type PyCObject_Type
	#define PyString_FromString PyUnicode_FromString
char* PyString_AsString(PyObject *ob) {
	PyObject *enc;
	char *cstr;
	enc = PyUnicode_AsEncodedString(ob, "utf-8", "Error");
	if( enc == NULL ) {
		PyErr_Format(PyExc_ValueError, "Cannot encode string");
		return NULL;
	}
	cstr = PyBytes_AsString(enc);
	Py_XDECREF(enc);
	return cstr;
}
#endif


/*
 Load in FFTW wisdom.  Based on the read_wisdom function in PRESTO.
*/

void read_wisdom(char *filename, PyObject *m) {
	int status = 0;
	FILE *wisdomfile;
	
	wisdomfile = fopen(filename, "r");
	if( wisdomfile != NULL ) {
		status = fftwf_import_wisdom_from_file(wisdomfile);
		fclose(wisdomfile);
	}
	PyModule_AddObject(m, "useWisdom", PyBool_FromLong(status));
}


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
		return sin(x*NPY_PI)/(x*NPY_PI);
	}
}


/*
  Complex magnitude squared functions
*/

double cabs2(double complex z) {
	return creal(z)*creal(z) + cimag(z)*cimag(z);
}

float cabs2f(float complex z) {
	return crealf(z)*crealf(z) + cimagf(z)*cimagf(z);
}


/*
  FFT Functions
    1. FPSDR - FFT a real-valued collection of signals
    2. FPSDC - FFT a complex-valued collection of signals
*/

static PyObject *FPSDR(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signalsX, *signalsY, *signalsF, *window=Py_None;
	PyArrayObject *dataX=NULL, *dataY=NULL, *dataF=NULL, *windowData=NULL;
	int nChan = 64;
	int Overlap = 1;
	int Clip = 0;
	
	long i, j, k, nStand, nSamps, nFFT;
	
	static char *kwlist[] = {"signalsX", "signalsY", "LFFT", "Overlap", "ClipLevel", "window", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OO|iiiO:set_callback", kwlist, &signalsX, &signalsY, &nChan, &Overlap, &Clip, &window)) {
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
	dataX = (PyArrayObject *) PyArray_ContiguousFromObject(signalsX, NPY_INT16, 2, 2);
	dataY = (PyArrayObject *) PyArray_ContiguousFromObject(signalsY, NPY_INT16, 2, 2);
	if( dataX == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input array signalsX to 2-D int16");
		goto fail;
	}
	if( dataY == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input array signalsY to 2-D int16");
		goto fail;
	}
	
	// Get the properties of the data
	nStand = (long) PyArray_DIM(dataX, 0);
	nSamps = (long) PyArray_DIM(dataX, 1);
	
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
		window = Py_BuildValue("(i)", 2*nChan);
		window = PyObject_CallObject(windowFunc, window);
		windowData = (PyArrayObject *) PyArray_ContiguousFromObject(window, NPY_DOUBLE, 1, 1);
		Py_DECREF(window);
	}
	
	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / ((2*nChan)/Overlap) - (2*nChan)/((2*nChan)/Overlap) + 1;
	npy_intp dims[3];
	dims[0] = (npy_intp) 4;
	dims[1] = (npy_intp) nStand;
	dims[2] = (npy_intp) nChan;
	dataF = (PyArrayObject*) PyArray_ZEROS(3, dims, NPY_DOUBLE, 0);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		goto fail;
	}
	
	Py_BEGIN_ALLOW_THREADS
	
	// Create the FFTW plan                          
	float *inP, *inX, *inY;                          
	float complex *outP, *outX, *outY;
	inP = (float *) fftwf_malloc(sizeof(float) * 2*nChan);
	outP = (float complex *) fftwf_malloc(sizeof(float complex) * (nChan+1));
	fftwf_plan p;
	p = fftwf_plan_dft_r2c_1d(2*nChan, inP, outP, FFTW_ESTIMATE);
	
	// Data indexing and access
	long secStart;
	short int *aX, *aY;
	double *b, *c;
	aX = (short int *) PyArray_DATA(dataX);
	aY = (short int *) PyArray_DATA(dataY);
	b = (double *) PyArray_DATA(dataF);
	if( windowData != NULL ) {
		c = (double *) PyArray_DATA(windowData);
	}
	
	// Time-domain blanking control
	double cleanFactor;
	long nActFFT;
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(inX, inY, outX, outY, i, j, k, secStart, cleanFactor, nActFFT)
	#endif
	{
		inX = (float *) fftwf_malloc(sizeof(float) * 2*nChan);
		inY = (float *) fftwf_malloc(sizeof(float) * 2*nChan);
		outX = (float complex *) fftwf_malloc(sizeof(float complex) * (nChan+1));
		outY = (float complex *) fftwf_malloc(sizeof(float complex) * (nChan+1));
		
		#ifdef _OPENMP
			#pragma omp for schedule(OMP_SCHEDULER)
		#endif
		for(i=0; i<nStand; i++) {
			nActFFT = 0;
			
			for(j=0; j<nFFT; j++) {
				cleanFactor = 1.0;
				secStart = nSamps * i + 2*nChan*j/Overlap;
				
				for(k=0; k<2*nChan; k++) {
					inX[k] = (float) *(aX + secStart + k);
					inY[k] = (float) *(aY + secStart + k);
					
					if( Clip && ( fabsf(inX[k]) >= Clip || fabsf(inY[k]) >= Clip ) ) {
						cleanFactor = 0.0;
					}
					
					if( windowData != NULL ) {
						inX[k] *= *(c + k);
						inY[k] *= *(c + k);
					}
				}
				
				fftwf_execute_dft_r2c(p, inX, outX);
				fftwf_execute_dft_r2c(p, inY, outY);
				
				for(k=0; k<nChan; k++) {
					// I
					*(b + 0*nChan*nStand + nChan*i + k) += cleanFactor*cabs2f(outX[k]);
					*(b + 0*nChan*nStand + nChan*i + k) += cleanFactor*cabs2f(outY[k]);
					
					// Q
					*(b + 1*nChan*nStand + nChan*i + k) += cleanFactor*cabs2f(outX[k]);
					*(b + 1*nChan*nStand + nChan*i + k) -= cleanFactor*cabs2f(outY[k]);
					
					// U
					*(b + 2*nChan*nStand + nChan*i + k) += 2*cleanFactor*crealf(outX[k])*crealf(outY[k]);
					*(b + 2*nChan*nStand + nChan*i + k) += 2*cleanFactor*cimagf(outX[k])*cimagf(outY[k]);
					
					// V
					*(b + 3*nChan*nStand + nChan*i + k) +=2*cleanFactor*cimagf(outX[k])*crealf(outY[k]);
					*(b + 3*nChan*nStand + nChan*i + k) -=2*cleanFactor*crealf(outX[k])*cimagf(outY[k]);
				}
				
				nActFFT += (long) cleanFactor;
			}
			
			// Scale FFTs
			for(j=0; j<4; j++) {
				cblas_dscal(nChan, 1.0/(2*nChan*nActFFT), (b + j*nChan*nStand + nChan*i), 1);
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

PyDoc_STRVAR(FPSDR_doc, \
"Perform a series of Fourier transforms with windows on real-valued data to\n\
get the PSD for the four Stokes parameters: I, Q, U, and V.\n\
\n\
Input arguments are:\n\
 * signals: 2-D numpy.int16 (stands by samples) array of data to FFT\n\
\n\
Input keywords are:\n\
 * LFFT: number of FFT channels to make (default=64)\n\
 * Overlap: number of overlapped FFTs to use (default=1)\n\
 * window: Callable Python function for generating the window or None for no window\n\
 * ClipLevel: count value of 'bad' data.  FFT windows with instantaneous powers\n\
   greater than or equal to this value greater are zeroed.  Setting the ClipLevel\n\
   to zero disables time-domain blanking\n\
\n\
Outputs:\n\
 * psd: 3-D numpy.double (Stokes parameter (I,Q,U,V) by stands by channels) of PSD data\n\
");


static PyObject *FPSDC(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signalsX, *signalsY, *signalsF, *window=Py_None;
	PyArrayObject *dataX=NULL, *dataY=NULL, *dataF=NULL, *windowData=NULL;
	int nChan = 64;
	int Overlap = 1;
	int Clip = 0;
	
	long i, j, k, nStand, nSamps, nFFT;
	
	static char *kwlist[] = {"signalsX", "signalsY", "LFFT", "Overlap", "ClipLevel", "window", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OO|iiiO:set_callback", kwlist, &signalsX, &signalsY, &nChan, &Overlap, &Clip, &window)) {
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
	dataX = (PyArrayObject *) PyArray_ContiguousFromObject(signalsX, NPY_COMPLEX64, 2, 2);
	dataY = (PyArrayObject *) PyArray_ContiguousFromObject(signalsY, NPY_COMPLEX64, 2, 2);
	if( dataX == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input array signalsX to 2-D complex64");
		goto fail;
	}
	if( dataY == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input array signalsY to 2-D complex64");
		goto fail;
	}
	
	// Get the properties of the data
	nStand = (long) PyArray_DIM(dataX, 0);
	nSamps = (long) PyArray_DIM(dataX, 1);
	
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
		window = Py_BuildValue("(i)", nChan);
		window = PyObject_CallObject(windowFunc, window);
		windowData = (PyArrayObject *) PyArray_ContiguousFromObject(window, NPY_DOUBLE, 1, 1);
		Py_DECREF(window);
	}
	
	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / (nChan/Overlap) - nChan/(nChan/Overlap) + 1;
	npy_intp dims[3];
	dims[0] = (npy_intp) 4;
	dims[1] = (npy_intp) nStand;
	dims[2] = (npy_intp) nChan;
	dataF = (PyArrayObject*) PyArray_ZEROS(3, dims, NPY_DOUBLE, 0);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		goto fail;
	}
	
	Py_BEGIN_ALLOW_THREADS
	
	// Create the FFTW plan
	float complex *inP, *inX, *inY;
	inP = (float complex*) fftwf_malloc(sizeof(float complex) * nChan);
	fftwf_plan p;
	p = fftwf_plan_dft_1d(nChan, inP, inP, FFTW_FORWARD, FFTW_ESTIMATE);
	
	// Data indexing and access
	long secStart;
	float complex *aX, *aY;
	double *b, *c, *temp2;
	aX = (float complex *) PyArray_DATA(dataX);
	aY = (float complex *) PyArray_DATA(dataY);
	b = (double *) PyArray_DATA(dataF);
	if( windowData != NULL ) {
		c = (double *) PyArray_DATA(windowData);
	}
	
	// Time-domain blanking control
	double cleanFactor;
	long nActFFT;
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(inX, inY, i, j, k, secStart, cleanFactor, nActFFT, temp2)
	#endif
	{
		inX = (float complex*) fftwf_malloc(sizeof(float complex) * nChan);
		inY = (float complex*) fftwf_malloc(sizeof(float complex) * nChan);
		temp2 = (double *) malloc(sizeof(double)*(nChan/2+nChan%2));
		
		#ifdef _OPENMP
			#pragma omp for schedule(OMP_SCHEDULER)
		#endif
		for(i=0; i<nStand; i++) {
			nActFFT = 0;
			
			for(j=0; j<nFFT; j++) {
				cleanFactor = 1.0;
				secStart = nSamps * i + nChan*j/Overlap;
				
				for(k=0; k<nChan; k++) {
					inX[k] = *(aX + secStart + k);
					inY[k] = *(aY + secStart + k);
					
					if( Clip && ( cabsf(inX[k]) >= Clip || cabsf(inY[k]) >= Clip ) ) {
						cleanFactor = 0.0;
					}
					
					if( windowData != NULL ) {
						inX[k] *= *(c + k);
						inY[k] *= *(c + k);
					}
				}
				
				fftwf_execute_dft(p, inX, inX);
				fftwf_execute_dft(p, inY, inY);
				
				for(k=0; k<nChan; k++) {
					// I
					*(b + 0*nChan*nStand + nChan*i + k) += cleanFactor*cabs2f(inX[k]);
					*(b + 0*nChan*nStand + nChan*i + k) += cleanFactor*cabs2f(inY[k]);
					
					// Q
					*(b + 1*nChan*nStand + nChan*i + k) += cleanFactor*cabs2f(inX[k]);
					*(b + 1*nChan*nStand + nChan*i + k) -= cleanFactor*cabs2f(inY[k]);
					
					// U
					*(b + 2*nChan*nStand + nChan*i + k) += 2*cleanFactor*crealf(inX[k])*crealf(inY[k]);
					*(b + 2*nChan*nStand + nChan*i + k) += 2*cleanFactor*cimagf(inX[k])*cimagf(inY[k]);
					
					// V
					*(b + 3*nChan*nStand + nChan*i + k) +=2*cleanFactor*cimagf(inX[k])*crealf(inY[k]);
					*(b + 3*nChan*nStand + nChan*i + k) -=2*cleanFactor*crealf(inX[k])*cimagf(inY[k]);
				}
				
				nActFFT += (long) cleanFactor;
			}
			
			for(j=0; j<4; j++) {
				// Shift FFTs
				memcpy(temp2, (b + j*nChan*nStand + nChan*i), sizeof(double)*(nChan/2+nChan%2));
				memmove((b + j*nChan*nStand + nChan*i), (b + j*nChan*nStand + nChan*i)+nChan/2+nChan%2, sizeof(double)*nChan/2);
				memcpy((b + j*nChan*nStand + nChan*i)+nChan/2, temp2, sizeof(double)*(nChan/2+nChan%2));
				
				// Scale FFTs
				cblas_dscal(nChan, 1.0/(nActFFT*nChan), (b + j*nChan*nStand + nChan*i), 1);
			}
		}
		
		fftwf_free(inX);
		fftwf_free(inY);
		free(temp2);
	}
	fftwf_destroy_plan(p);
	fftwf_free(inP);
	
	Py_END_ALLOW_THREADS
	
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

PyDoc_STRVAR(FPSDC_doc, \
"Perform a series of Fourier transforms with windows on complex-valued data\n\
to get the PSD for the four Stokes parameters: I, Q, U, and V.\n\
\n\
Input arguments are:\n\
 * signals: 2-D numpy.complex64 (stands by samples) array of data to FFT\n\
\n\
Input keywords are:\n\
 * LFFT: number of FFT channels to make (default=64)\n\
 * Overlap: number of overlapped FFTs to use (default=1)\n\
 * window: Callable Python function for generating the window or None for no window\n\
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
	if( dataX == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input signalsX array to 3-D complex64");
		goto fail;
	}
	if( dataY == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input signalsY array to 3-D complex64");
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
	float complex tempVis1, tempVis2;
	float complex *a, *b, *v;
	a = (float complex *) PyArray_DATA(dataX);
	b = (float complex *) PyArray_DATA(dataY);
	v = (float complex *) PyArray_DATA(vis);
	
	// Time-domain blanking control
	long nActVis;
	unsigned char *u1, *u2;
	u1 = (unsigned char *) PyArray_DATA(validX);
	u2 = (unsigned char *) PyArray_DATA(validY);
	
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
	{"FPSDR",    (PyCFunction) FPSDR,    METH_VARARGS|METH_KEYWORDS, FPSDR_doc   }, 
	{"FPSDC",    (PyCFunction) FPSDC,    METH_VARARGS|METH_KEYWORDS, FPSDC_doc   }, 
	{"XEngine2", (PyCFunction) XEngine2, METH_VARARGS,               XEngine2_doc}, 
	{NULL,       NULL,                   0,                          NULL        }
};

PyDoc_STRVAR(stokes_doc, \
"Extension to take X and Y timeseries data and create the four Stokes\n\
parameters.\n\
\n\
The functions defined in this module are:\n\
  * FPSDR -  FFT and integrate function for computing a series of overlapped\n\
    Fourier transforms for a real-valued (TBW) signal from a collection of\n\
    stands all at once.\n\
  * FPSDC - FFT and integrate function for computing a series of overlapped\n\
    Fourier transforms for a complex-valued (TBN and DRX) signal from a \n\
    collection of stands/beams all at once.\n\
\n\
Also included is an X-Engine for use with the lsl.correlator._core module to\n\
perform cross-correlations for the stokes parameters.\n\
\n\
See the inidividual functions for more details.");


/*
  Module Setup - Initialization
*/

#if PY_MAJOR_VERSION >= 3
	#define MOD_ERROR_VAL NULL
	#define MOD_SUCCESS_VAL(val) val
	#define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
	#define MOD_DEF(ob, name, methods, doc) \
	   static struct PyModuleDef moduledef = { \
	      PyModuleDef_HEAD_INIT, name, doc, -1, methods, }; \
	   ob = PyModule_Create(&moduledef);
#else
	#define MOD_ERROR_VAL
	#define MOD_SUCCESS_VAL(val)
	#define MOD_INIT(name) PyMODINIT_FUNC init##name(void)
	#define MOD_DEF(ob, name, methods, doc) \
	   ob = Py_InitModule3(name, methods, doc);
#endif

MOD_INIT(_stokes) {
	char filename[256];
	PyObject *m, *pModule, *pDataPath=NULL;
	
	Py_Initialize();
	
	// Module definitions and functions
	MOD_DEF(m, "_stokes", StokesMethods, stokes_doc);
	if( m == NULL ) {
		return MOD_ERROR_VAL;
	}
	import_array();
	
	// Version and revision information
	PyModule_AddObject(m, "__version__", PyString_FromString("0.3"));
	PyModule_AddObject(m, "__revision__", PyString_FromString("$Rev$"));
	
	// LSL FFTW Wisdom
	pModule = PyImport_ImportModule("lsl.common.paths");
	if( pModule != NULL ) {
		pDataPath = PyObject_GetAttrString(pModule, "data");
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

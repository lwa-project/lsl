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
	#define PyString_AsString PyBytes_AsString
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
	PyObject *signals, *signalsF, *window=Py_None;
	PyArrayObject *data=NULL, *dataF=NULL, *windowData=NULL;
	int nChan = 64;
	int Overlap = 1;
	int Clip = 0;
	
	long i, j, k, nStand, nSamps, nFFT;
	
	static char *kwlist[] = {"signals", "LFFT", "Overlap", "ClipLevel", "window", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|iiiO:set_callback", kwlist, &signals, &nChan, &Overlap, &Clip, &window)) {
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
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_INT16, 2, 2);
	if( data == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input array signals to 2-D int16");
		goto fail;
	}
	
	// Get the properties of the data
	nStand = (long) PyArray_DIM(data, 0);
	nSamps = (long) PyArray_DIM(data, 1);
	
	// Calculate the windowing function
	if( windowFunc != Py_None ) {
		window = Py_BuildValue("(i)", 2*nChan);
		window = PyObject_CallObject(windowFunc, window);
		windowData = (PyArrayObject *) PyArray_ContiguousFromObject(window, NPY_DOUBLE, 1, 1);
		Py_DECREF(window);
	}
	
	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / ((2*nChan)/Overlap) - (2*nChan)/((2*nChan)/Overlap) + 1;
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChan;
	dataF = (PyArrayObject*) PyArray_ZEROS(2, dims, NPY_DOUBLE, 0);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		goto fail;
	}
	
	Py_BEGIN_ALLOW_THREADS
	
	// Create the FFTW plan                          
	float *inP, *in;                          
	float complex *outP, *out;
	inP = (float *) fftwf_malloc(sizeof(float) * 2*nChan);
	outP = (float complex *) fftwf_malloc(sizeof(float complex) * (nChan+1));
	fftwf_plan p;
	p = fftwf_plan_dft_r2c_1d(2*nChan, inP, outP, FFTW_ESTIMATE);
	
	// Data indexing and access
	long secStart;
	short int *a;
	double *b, *c;
	a = (short int *) PyArray_DATA(data);
	b = (double *) PyArray_DATA(dataF);
	if( windowData != NULL ) {
		c = (double *) PyArray_DATA(windowData);
	}
	
	// Time-domain blanking control
	double cleanFactor;
	long nActFFT;
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(in, out, i, j, k, secStart, cleanFactor, nActFFT)
	#endif
	{
		in = (float *) fftwf_malloc(sizeof(float) * 2*nChan);
		out = (float complex *) fftwf_malloc(sizeof(float complex) * (nChan+1));
		
		#ifdef _OPENMP
			#pragma omp for schedule(OMP_SCHEDULER)
		#endif
		for(i=0; i<nStand; i++) {
			nActFFT = 0;
			
			for(j=0; j<nFFT; j++) {
				cleanFactor = 1.0;
				secStart = nSamps * i + 2*nChan*j/Overlap;
				
				for(k=0; k<2*nChan; k++) {
					in[k] = (float) *(a + secStart + k);
					
					if( Clip && fabsf(in[k]) >= Clip ) {
						cleanFactor = 0.0;
					}
					
					if( windowData != NULL ) {
						in[k] *= *(c + k);
					}
				}
				
				fftwf_execute_dft_r2c(p, in, out);
				
				for(k=0; k<nChan; k++) {
					*(b + nChan*i + k) += cleanFactor*cabs2f(out[k]);
				}
				
				nActFFT += (long) cleanFactor;
			}
			
			// Scale FFTs
			cblas_dscal(nChan, 1.0/(2*nChan*nActFFT), (b + i*nChan), 1);
		}
		
		fftwf_free(in);
		fftwf_free(out);
	}
	fftwf_destroy_plan(p);
	fftwf_free(inP);
	fftwf_free(outP);
	
	Py_END_ALLOW_THREADS
	
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

PyDoc_STRVAR(FPSDR_doc, \
"Perform a series of Fourier transforms with windows on real-valued data to\n\
get the PSD.\n\
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
 * psd: 2-D numpy.double (stands by channels) of PSD data\n\
");


static PyObject *FPSDC(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *signalsF, *window=Py_None;
	PyArrayObject *data=NULL, *dataF=NULL, *windowData=NULL;
	int nChan = 64;
	int Overlap = 1;
	int Clip = 0;

	long i, j, k, nStand, nSamps, nFFT;

	static char *kwlist[] = {"signals", "LFFT", "Overlap", "ClipLevel", "window", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|iiiO:set_callback", kwlist, &signals, &nChan, &Overlap, &Clip, &window)) {
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
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_COMPLEX64, 2, 2);
	if( data == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input array signals to 2-D complex64");
		goto fail;
	}
	
	// Get the properties of the data
	nStand = (long) PyArray_DIM(data, 0);
	nSamps = (long) PyArray_DIM(data, 1);
	
	// Calculate the windowing function
	if( windowFunc != Py_None ) {
		window = Py_BuildValue("(i)", nChan);
		window = PyObject_CallObject(windowFunc, window);
		windowData = (PyArrayObject *) PyArray_ContiguousFromObject(window, NPY_DOUBLE, 1, 1);
		Py_DECREF(window);
	}
	
	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / (nChan/Overlap) - nChan/(nChan/Overlap) + 1;
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChan;
	dataF = (PyArrayObject*) PyArray_ZEROS(2, dims, NPY_DOUBLE, 0);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
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
	float complex *a;
	double *b, *c, *temp2;
	a = (float complex *) PyArray_DATA(data);
	b = (double *) PyArray_DATA(dataF);
	if( windowData != NULL ) {
		c = (double *) PyArray_DATA(windowData);
	}
	
	// Time-domain blanking control
	double cleanFactor;
	long nActFFT;
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(in, i, j, k, secStart, cleanFactor, nActFFT, temp2)
	#endif
	{
		in = (float complex*) fftwf_malloc(sizeof(float complex) * nChan);
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
					in[k] = *(a + secStart + k);
					
					if( Clip && cabsf(in[k]) >= Clip ) {
						cleanFactor = 0.0;
					}
					
					if( windowData != NULL ) {
						in[k] *= *(c + k);
					}
				}
				
				fftwf_execute_dft(p, in, in);
				
				for(k=0; k<nChan; k++) {
					*(b + nChan*i + k) += cleanFactor*cabs2f(in[k]);
				}
				
				nActFFT += (long) cleanFactor;
			}
			
			// Shift FFTs
			memcpy(temp2, b + nChan*i, sizeof(double)*(nChan/2+nChan%2));
			memmove(b + nChan*i, b + nChan*i+nChan/2+nChan%2, sizeof(double)*nChan/2);
			memcpy(b + nChan*i+nChan/2, temp2, sizeof(double)*(nChan/2+nChan%2));
			
			// Scale FFTs
			cblas_dscal(nChan, 1.0/(nActFFT*nChan), (b + i*nChan), 1);
		}
		
		fftwf_free(in);
		free(temp2);
	}
	fftwf_destroy_plan(p);
	fftwf_free(inP);
	
	Py_END_ALLOW_THREADS
	
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

PyDoc_STRVAR(FPSDC_doc, \
"Perform a series of Fourier transforms with windows on complex-valued data\n\
to get the PSD.\n\
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
 * psd: 2-D numpy.double (stands by channels) of PSD data\n\
");


/*
  Module Setup - Function Definitions and Documentation
*/

static PyMethodDef SpecMethods[] = {
	{"FPSDR", (PyCFunction) FPSDR, METH_VARARGS|METH_KEYWORDS, FPSDR_doc}, 
	{"FPSDC", (PyCFunction) FPSDC, METH_VARARGS|METH_KEYWORDS, FPSDC_doc}, 
	{NULL,    NULL,                0,                          NULL     }
};

PyDoc_STRVAR(spec_doc, \
"Extension to take timeseries data and convert it to the frequency domain.\n\
\n\
The functions defined in this module are:\n\
  * FPSDR -  FFT and integrate function for computing a series of overlapped\n\
    Fourier transforms for a real-valued (TBW) signal from a collection of\n\
    stands all at once.\n\
  * FPSDC - FFT and integrate function for computing a series of overlapped\n\
    Fourier transforms for a complex-valued (TBN and DRX) signal from a \n\
    collection of stands/beams all at once.\n\
\n\
See the inidividual functions for more details.\n\
\n\
.. versionchanged:: 1.0.1\n\
\tRemoved the polyphase filterbank versions of the four core functions.\n\
");



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

MOD_INIT(_spec) {
	char filename[256];
	PyObject *m, *pModule, *pDataPath;

	// Module definitions and functions
	MOD_DEF(m, "_spec", SpecMethods, spec_doc);
	if( m == NULL ) {
		return MOD_ERROR_VAL;
	}
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
	
	#if PY_MAJOR_VERSION >= 3
		return m;
	#endif
}

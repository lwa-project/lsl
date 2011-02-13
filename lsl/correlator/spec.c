#include "Python.h"
#include <math.h>
#include <stdio.h>
#include <fftw3.h>
#include <stdlib.h>
#include <complex.h>

#ifdef _OPENMP
        #include <omp.h>
#endif

#include "numpy/arrayobject.h"

#define PI 3.1415926535898
#define imaginary _Complex_I


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


static PyObject *FPSDR2(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *signalsF;
	PyArrayObject *data, *dataF;
	int nChan = 64;
	int Overlap = 1;

	long i, j, k, nStand, nSamps, nFFT;
	
	static char *kwlist[] = {"signals", "LFFT", "Overlap", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|ii", kwlist, &signals, &nChan, &Overlap)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it useable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, PyArray_DOUBLE, 2, 2);
	
	// Get the properties of the data
	nStand = (long) data->dimensions[0];
	nSamps = (long) data->dimensions[1];
	
	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / 2 / nChan * Overlap - Overlap + 1;
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) (nChan - 1);
	dataF = (PyArrayObject*) PyArray_SimpleNew(2, dims, PyArray_DOUBLE);
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
	double *a, *b;
	a = (double *) data->data;
	b = (double *) dataF->data;
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(in, secStart, j, k, fftIndex)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * 2*nChan);
			
			for(j=0; j<nFFT; j++) {
				secStart = (long) (2*nChan*((float) j)/Overlap);
				
				for(k=0; k<2*nChan; k++) {
					in[k][0] = *(a + nSamps*i + secStart + k);
					in[k][1] = 0.0;
				}
				
				fftw_execute_dft(p, in, in);
				
				for(k=0; k<(nChan-1); k++) {
					fftIndex = k + 1;
					*(b + (nChan-1)*i + k) += in[fftIndex][0]*in[fftIndex][0] / 2 / nChan;
					*(b + (nChan-1)*i + k) += in[fftIndex][1]*in[fftIndex][1] / 2 / nChan;
				}
			}

			fftw_free(in);
		}
	}
	fftw_destroy_plan(p);
	fftw_free(inP);

	Py_XDECREF(data);

	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(FPSDR2_doc, "Perform a series of Fourier transforms on real-valued data to get the PSD");


static PyObject *FPSDR3(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *signalsF, *window;
	PyArrayObject *data, *dataF, *windowData;
	int nChan = 64;
	int Overlap = 1;

	long i, j, k, nStand, nSamps, nFFT;
	
	static char *kwlist[] = {"signals", "LFFT", "Overlap", "window", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|iiO", kwlist, &signals, &nChan, &Overlap,&window)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it useable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, PyArray_DOUBLE, 2, 2);
	windowData = (PyArrayObject *) PyArray_ContiguousFromObject(window, PyArray_DOUBLE, 1, 1);
	
	// Check data dimensions
	if( windowData->dimensions[0] != 2*nChan) {
		PyErr_Format(PyExc_TypeError, "window has a different channel count than twice nChan");
		Py_XDECREF(data);
		Py_XDECREF(windowData);
		return NULL;
	}
	
	// Get the properties of the data
	nStand = (long) data->dimensions[0];
	nSamps = (long) data->dimensions[1];

	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / 2 / nChan * Overlap - Overlap + 1;
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) (nChan - 1);
	dataF = (PyArrayObject*) PyArray_SimpleNew(2, dims, PyArray_DOUBLE);
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
	double *a, *b, *c;
	a = (double *) data->data;
	b = (double *) dataF->data;
	c = (double *) windowData->data;
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(in, secStart, j, k, fftIndex)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * 2*nChan);
			
			for(j=0; j<nFFT; j++) {
				secStart = (long) (2*nChan*((float) j)/Overlap);
				
				for(k=0; k<2*nChan; k++) {
					in[k][0] = *(a + nSamps*i + secStart + k) * *(c + k);
					in[k][1] = 0.0;
				}
				
				fftw_execute_dft(p, in, in);
				
				for(k=0; k<(nChan-1); k++) {
					fftIndex = k + 1;
					*(b + (nChan-1)*i + k) += in[fftIndex][0]*in[fftIndex][0] / 2 / nChan;
					*(b + (nChan-1)*i + k) += in[fftIndex][1]*in[fftIndex][1] / 2 / nChan;
				}
			}

			fftw_free(in);
		}
	}
	fftw_destroy_plan(p);
	fftw_free(inP);

	Py_XDECREF(data);
	Py_XDECREF(windowData);

	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(FPSDR3_doc, "Perform a series of Fourier transforms with windows on real-valued data to get the PSD.");


static PyObject *FPSDC2(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *signalsF;
	PyArrayObject *data, *dataF;
	int nChan = 64;
	int Overlap = 1;

	long i, j, k, nStand, nSamps, nFFT;

	static char *kwlist[] = {"signals", "LFFT", "Overlap", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|ii", kwlist, &signals, &nChan, &Overlap)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it useable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, PyArray_CDOUBLE, 2, 2);

	// Get the properties of the data
	nStand = (long) data->dimensions[0];
	nSamps = (long) data->dimensions[1];

	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / nChan * Overlap - Overlap + 1;
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) (nChan - 1);
	dataF = (PyArrayObject*) PyArray_SimpleNew(2, dims, PyArray_DOUBLE);
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
	double complex *a;
	double *b;
	a = (double complex *) data->data;
	b = (double *) dataF->data;
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(in, secStart, j, k, fftIndex)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nChan);
			
			for(j=0; j<nFFT; j++) {
				secStart = (long) (nChan*((float) j)/Overlap);
				
				for(k=0; k<nChan; k++) {
					in[k][0] = creal( *(a + nSamps*i + secStart + k) );
					in[k][0] = cimag( *(a + nSamps*i + secStart + k) );
				}
				
				fftw_execute_dft(p, in, in);
				
				for(k=0; k<(nChan-1); k++) {
					fftIndex = ((k+1) + nChan/2) % nChan;
					*(b + (nChan-1)*i + k) += in[fftIndex][0]*in[fftIndex][0] / nChan + in[fftIndex][1]*in[fftIndex][1] / nChan;
				}
			}

			fftw_free(in);
		}
	}
	fftw_destroy_plan(p);
	fftw_free(inP);

	Py_XDECREF(data);

	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(FPSDC2_doc, "Perform a series of Fourier transforms on complex-valued data to get the PSD");


static PyObject *FPSDC3(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *signalsF, *window;
	PyArrayObject *data, *dataF, *windowData;
	int nChan = 64;
	int Overlap = 1;
	double SampleRate = 1.0e5;

	long i, j, k, nStand, nSamps, nFFT;

	static char *kwlist[] = {"signals", "LFFT", "Overlap", "window", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|iiO", kwlist, &signals, &nChan, &Overlap, &window)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it useable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, PyArray_CDOUBLE, 2, 2);
	windowData = (PyArrayObject *) PyArray_ContiguousFromObject(window, PyArray_DOUBLE, 1, 1);
	
	// Check data dimensions
	if( windowData->dimensions[0] != nChan) {
		PyErr_Format(PyExc_TypeError, "window has a different channel count than nChan");
		Py_XDECREF(data);
		Py_XDECREF(windowData);
		return NULL;
	}

	// Get the properties of the data
	nStand = (long) data->dimensions[0];
	nSamps = (long) data->dimensions[1];

	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / nChan * Overlap - Overlap + 1;
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) (nChan - 1);
	dataF = (PyArrayObject*) PyArray_SimpleNew(2, dims, PyArray_CDOUBLE);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
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
	double complex *a;
	double *b, *c;
	a = (double complex *) data->data;
	b = (double *) dataF->data;
	c = (double *) windowData->data;
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(in, secStart, j, k, fftIndex)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nChan);
			
			for(j=0; j<nFFT; j++) {
				secStart = (long) (nChan*((float) j)/Overlap);
				
				for(k=0; k<nChan; k++) {
					in[k][0] = creal( *(a + nSamps*i + secStart + k) ) * *(c + k);
					in[k][0] = cimag( *(a + nSamps*i + secStart + k) ) * *(c + k);
				}
				
				fftw_execute_dft(p, in, in);
				
				for(k=0; k<(nChan-1); k++) {
					fftIndex = ((k+1) + nChan/2) % nChan;
					*(b + (nChan-1)*i + k) += in[fftIndex][0]*in[fftIndex][0] / nChan + in[fftIndex][1]*in[fftIndex][1] / nChan;
				}
			}

			fftw_free(in);
		}
	}
	fftw_destroy_plan(p);
	fftw_free(inP);

	Py_XDECREF(data);
	Py_XDECREF(windowData);

	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(FPSDC3_doc, "Perform a series of Fourier transforms with windows on complex-valued data to get the PSD");


static PyObject *PPSDR2(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *signalsF;
	PyArrayObject *data, *dataF;
	int nChan = 64;
	int nTaps = 4;
	int Overlap = 1;

	long i, j, k, nStand, nSamps, nFFT;
	
	static char *kwlist[] = {"signals", "LFFT", "Overlap", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|ii", kwlist, &signals, &nChan, &Overlap)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it useable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, PyArray_DOUBLE, 2, 2);
	
	// Get the properties of the data
	nStand = (long) data->dimensions[0];
	nSamps = (long) data->dimensions[1];

	// Compute the filterbank window for the correct numer of taps
	double fbWindow[2*nChan*nTaps];
	double complex tempFB[nChan-1];
	for(i=0; i<2*nChan*nTaps; i++) {
		fbWindow[i] = sinc((double) (i - nTaps*nChan + 0.5)/2/nChan) / 2 / nChan;
	}
	for(i=0; i<(nChan-1); i++) {
			tempFB[i] = 0.0;
	}
	
	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / 2 / nChan * Overlap - Overlap + 1;
	nFFT = (nFFT / nTaps) * nTaps;
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) (nChan - 1);
	dataF = (PyArrayObject*) PyArray_SimpleNew(2, dims, PyArray_DOUBLE);
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
	double *a, *b;
	a = (double *) data->data;
	b = (double *) dataF->data;
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(in, secStart, tempFB, j, k, fftIndex)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * 2*nChan);
			
			for(j=0; j<nFFT; j++) {
				if(j % nTaps == 0) {
					for(k=0; k<(nChan-1); k++) {
						*(b + (nChan-1)*i + k) += cabs(tempFB[k]) / 2 / nChan;
						tempFB[k] = 0;
					}
				}

				secStart = (long) (2*nChan*((float) j)/Overlap);
				
				for(k=0; k<2*nChan; k++) {
					in[k][0] = *(a + nSamps*i + secStart + k) * fbWindow[2*nChan*(j % nTaps) + k];
					in[k][1] = 0.0;
				}
				
				fftw_execute_dft(p, in, in);
				
				for(k=0; k<(nChan-1); k++) {
					fftIndex = k + 1;
					tempFB[k] = in[fftIndex][0] + imaginary * in[fftIndex][1];
				}
			}

			fftw_free(in);
		}
	}
	fftw_destroy_plan(p);
	fftw_free(inP);

	Py_XDECREF(data);

	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(PPSDR2_doc, "Perform a series of filter bank transforms on real-valued data to get the PSD");


static PyObject *PPSDR3(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *signalsF, *window;
	PyArrayObject *data, *dataF, *windowData;
	int nChan = 64;
	int nTaps = 4;
	int Overlap = 1;

	long i, j, k, nStand, nSamps, nFFT;
	
	static char *kwlist[] = {"signals", "LFFT", "Overlap", "window", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|iiO", kwlist, &signals, &nChan, &Overlap, &window)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it useable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, PyArray_DOUBLE, 2, 2);
	windowData = (PyArrayObject *) PyArray_ContiguousFromObject(window, PyArray_DOUBLE, 1, 1);
	
	// Check data dimensions
	if( windowData->dimensions[0] != 2*nTaps*nChan) {
		PyErr_Format(PyExc_TypeError, "window has a different channel count than 2*nTaps*nChan");
		Py_XDECREF(data);
		Py_XDECREF(windowData);
		return NULL;
	}
	
	// Get the properties of the data
	nStand = (long) data->dimensions[0];
	nSamps = (long) data->dimensions[1];

	// Compute the filterbank window for the correct numer of taps
	npy_intp *qLoc;
	double fbWindow[2*nChan*nTaps];
	double complex tempFB[nChan-1];
	qLoc = PyDimMem_NEW(1);
	for(i=0; i<2*nChan*nTaps; i++) {
		qLoc[0] = (npy_intp) i;
		fbWindow[i] = sinc((double) (i - nTaps*nChan + 0.5)/2/nChan) / 2 / nChan;
		fbWindow[i] *= *(double *) PyArray_GetPtr(windowData, qLoc);
	}
	for(i=0; i<(nChan-1); i++) {
			tempFB[i] = 0.0;
	}
	PyDimMem_FREE(qLoc);

	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / 2 / nChan * Overlap - Overlap + 1;
	nFFT = (nFFT / nTaps) * nTaps;
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) (nChan - 1);
	dataF = (PyArrayObject*) PyArray_SimpleNew(2, dims, PyArray_DOUBLE);
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
	double *a, *b, *c;
	a = (double *) data->data;
	b = (double *) dataF->data;
	c = (double *) windowData->data;
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(in, secStart, tempFB, j, k, fftIndex)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * 2*nChan);
			
			for(j=0; j<nFFT; j++) {
				if(j % nTaps == 0) {
					for(k=0; k<(nChan-1); k++) {
						*(b + (nChan-1)*i + k) += cabs(tempFB[k]) / 2 / nChan;
						tempFB[k] = 0;
					}
				}

				secStart = (long) (2*nChan*((float) j)/Overlap);
				
				for(k=0; k<2*nChan; k++) {
					in[k][0] = *(a + nSamps*i + secStart + k) * fbWindow[2*nChan*(j % nTaps) + k];
					in[k][1] = 0.0;
				}
				
				fftw_execute_dft(p, in, in);
				
				for(k=0; k<(nChan-1); k++) {
					fftIndex = k + 1;
					tempFB[k] = in[fftIndex][0] + imaginary * in[fftIndex][1];
				}
			}
			fftw_free(in);
		}
	}
	fftw_destroy_plan(p);
	fftw_free(inP);

	Py_XDECREF(data);
	Py_XDECREF(windowData);

	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(PPSDR3_doc, "Perform a series of filter bank transforms with windows on real-valued data to get the PSD.");


static PyObject *PPSDC2(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *signalsF;
	PyArrayObject *data, *dataF;
	int nChan = 64;
	int nTaps = 4;
	int Overlap = 1;

	long i, j, k, nStand, nSamps, nFFT;

	static char *kwlist[] = {"signals", "LFFT", "Overlap", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|ii", kwlist, &signals, &nChan, &Overlap)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it useable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, PyArray_CDOUBLE, 2, 2);

	// Get the properties of the data
	nStand = (long) data->dimensions[0];
	nSamps = (long) data->dimensions[1];

	// Compute the filterbank window for the correct numer of taps
	double fbWindow[nChan*nTaps];
	double complex tempFB[nChan-1];
	for(i=0; i<nChan*nTaps; i++) {
		fbWindow[i] = sinc((double) (i - nTaps*nChan/2 + 0.5)/nChan) / nChan;
	}
	for(i=0; i<(nChan-1); i++) {
		tempFB[i] = 0.0;
	}

	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / nChan * Overlap - Overlap + 1;
	nFFT = (nFFT / nTaps) * nTaps;
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) (nChan - 1);
	dataF = (PyArrayObject*) PyArray_SimpleNew(2, dims, PyArray_DOUBLE);
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
	double complex *a;
	double *b;
	a = (double complex *) data->data;
	b = (double *) dataF->data;
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(in, secStart, tempFB, j, k, fftIndex)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nChan);
			
			for(j=0; j<nFFT; j++) {
				if(j % nTaps == 0) {
					for(k=0; k<(nChan-1); k++) {
						*(b + (nChan-1)*i + k) += cabs(tempFB[k]) / nChan;
						tempFB[k] = 0.0;
					}
				}

				secStart = (long) (nChan*((float) j)/Overlap);
				
				for(k=0; k<nChan; k++) {
					in[k][0] = creal( *(a + nSamps*i + secStart + k) ) * fbWindow[nChan*(j % nTaps) + k];
					in[k][1] = cimag( *(a + nSamps*i + secStart + k) ) * fbWindow[nChan*(j % nTaps) + k];
				}
				
				fftw_execute_dft(p, in, in);
				
				for(k=0; k<(nChan-1); k++) {
					fftIndex = ((k+1) + nChan/2) % nChan;
					tempFB[k] = in[fftIndex][0] + imaginary * in[fftIndex][1];
				}
			}

			fftw_free(in);
		}
	}
	fftw_destroy_plan(p);
	fftw_free(inP);

	Py_XDECREF(data);

	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(PPSDC2_doc, "Perform a series of filter bank transforms on complex-valued data to get the PSD");


static PyObject *PPSDC3(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *signalsF, *window;
	PyArrayObject *data, *dataF, *windowData;
	int nChan = 64;
	int nTaps = 4;
	int Overlap = 1;

	long i, j, k, nStand, nSamps, nFFT;

	static char *kwlist[] = {"signals","LFFT", "Overlap", "window", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|iiO", kwlist, &signals, &nChan, &Overlap, &window)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it useable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, PyArray_CDOUBLE, 2, 2);
	windowData = (PyArrayObject *) PyArray_ContiguousFromObject(window, PyArray_DOUBLE, 1, 1);
	
	// Check data dimensions
	if( windowData->dimensions[0] != nTaps*nChan) {
		PyErr_Format(PyExc_TypeError, "window has a different channel count than nTaps*nChan");
		Py_XDECREF(data);
		Py_XDECREF(windowData);
		return NULL;
	}

	// Get the properties of the data
	nStand = (long) data->dimensions[0];
	nSamps = (long) data->dimensions[1];

	// Compute the filterbank window for the correct numer of taps
	npy_intp *qLoc;
	double fbWindow[nChan*nTaps];
	double complex tempFB[nChan-1];
	qLoc = PyDimMem_NEW(1);
	for(i=0; i<nChan*nTaps; i++) {
		qLoc[0] = (npy_intp) i;
		fbWindow[i] = sinc((double) (i - nTaps*nChan/2 + 0.5)/nChan) / nChan;
		fbWindow[i] *= *(double *) PyArray_GetPtr(windowData, qLoc);
	}
	for(i=0; i<(nChan-1); i++) {
			tempFB[i] = 0.0;
	}
	PyDimMem_FREE(qLoc);

	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / nChan * Overlap - Overlap + 1;
	nFFT = (nFFT / nTaps) * nTaps;
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) (nChan - 1);
	dataF = (PyArrayObject*) PyArray_SimpleNew(2, dims, PyArray_CDOUBLE);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
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
	double complex *a;
	double *b;
	a = (double complex *) data->data;
	b = (double *) dataF->data;
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(in, secStart, j, k, fftIndex)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nStand; i++) {
			in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nChan);
			
			for(j=0; j<nFFT; j++) {
				if(j % nTaps == 0) {
					for(k=0; k<(nChan-1); k++) {
						*(b + (nChan-1)*i + k) += cabs(tempFB[k]) / nChan;
						tempFB[k] = 0;
					}
				}

				secStart = (long) (nChan*((float) j)/Overlap);
				
				for(k=0; k<nChan; k++) {
					in[k][0] = creal( *(a + nSamps*i + secStart + k) ) * fbWindow[nChan*(j % nTaps) + k];
					in[k][1] = cimag( *(a + nSamps*i + secStart + k) ) * fbWindow[nChan*(j % nTaps) + k];
				}
				
				fftw_execute_dft(p, in, in);
				
				for(k=0; k<(nChan-1); k++) {
					fftIndex = ((k+1) + nChan/2) % nChan;
					tempFB[k] = in[fftIndex][0] + imaginary * in[fftIndex][1];
				}
			}

			fftw_free(in);
		}
	}
	fftw_destroy_plan(p);
	fftw_free(inP);

	Py_XDECREF(data);
	Py_XDECREF(windowData);

	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalsF;
}

PyDoc_STRVAR(PPSDC3_doc, "Perform a series of filter bank transforms with windows on complex-valued data to get the PSD");


/*
  Module Setup - Function Definitions and Documentation
*/

static PyMethodDef SpecMethods[] = {
	{"FPSDR2",  FPSDR2,  METH_KEYWORDS, FPSDR2_doc}, 
	{"FPSDR3",  FPSDR3,  METH_KEYWORDS, FPSDR3_doc}, 
	{"FPSDC2",  FPSDC2,  METH_KEYWORDS, FPSDC2_doc}, 
	{"FPSDC3",  FPSDC3,  METH_KEYWORDS, FPSDC3_doc}, 
	{"PPSDR2",  PPSDR2,  METH_KEYWORDS, PPSDR2_doc}, 
	{"PPSDR3",  PPSDR3,  METH_KEYWORDS, PPSDR3_doc}, 
	{"PPSDC2",  PPSDC2,  METH_KEYWORDS, PPSDC2_doc}, 
	{"PPSDC3",  PPSDC3,  METH_KEYWORDS, PPSDC3_doc}, 
	{NULL, NULL, 0, NULL}
};

PyDoc_STRVAR(spec_doc, \
"Test extension to replace lsl.correlator.fx.calcSpectra");


/*
  Module Setup - Initialization
*/

PyMODINIT_FUNC init_spec(void) {
	(void) Py_InitModule3("_spec", SpecMethods, spec_doc);
	import_array();
}

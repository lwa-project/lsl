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

static PyObject *FEngineR(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signal, *freq, *delay, *signalF;
	PyArrayObject *data, *fq, *times, *dataF;
	int nChan = 64;
	int Overlap = 1;
	double SampleRate = 196.0e6;

	long i, j, nFFT, start;
	
	static char *kwlist[] = {"signal", "freq", "delay", "LFFT", "Overlap", "SampleRate", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OOO|iid", kwlist, &signal, &freq, &delay, &nChan, &Overlap, &SampleRate)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it useable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signal, PyArray_DOUBLE, 1, 1);
	fq = (PyArrayObject *) PyArray_ContiguousFromObject(freq, PyArray_DOUBLE, 1, 1);
	times = (PyArrayObject *) PyArray_ContiguousFromObject(delay, PyArray_DOUBLE, 1, 1);
	
	// Check data dimensions
	if(nChan != (fq->dimensions[0]+1)) {
		PyErr_Format(PyExc_RuntimeError, "freq has a different channel count than nChan");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		return NULL;
	}
	
	if(fq->dimensions[0] != times->dimensions[0]) {
		PyErr_Format(PyExc_TypeError, "freq and delays have different channel counts");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		return NULL;
	}
	
	// Compute the interger sample offset and the fractional sample delay
	npy_intp *tLoc;
	tLoc = PyDimMem_NEW(1);
	tLoc[0] = (npy_intp) fq->dimensions[0] / 2;
	start = (long) round(*(double *) PyArray_GetPtr(times, tLoc) * SampleRate);
	double frac[nChan];
	for(i=0; i<nChan; i++) {
		tLoc[0] = (npy_intp) i;
		frac[i] = *(double *) PyArray_GetPtr(times, tLoc) - start/SampleRate;
	}
	PyDimMem_FREE(tLoc);

	// Find out how large the output array needs to be and initialize it
	nFFT = (data->dimensions[0] - start) / 2 / nChan * Overlap - Overlap + 1;
	npy_intp dims[2];
	dims[0] = nChan - 1;
	dims[1] = nFFT;
	dataF = (PyArrayObject*) PyArray_SimpleNew(2, dims, PyArray_CDOUBLE);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		return NULL;
	}

	// Create the FFTW plan
	fftw_complex *in, *out;
	fftw_plan p;
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 2*nChan);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 2*nChan);
	p = fftw_plan_dft_1d(2*nChan, in, out, FFTW_FORWARD, FFTW_MEASURE);
	
	// Part 1:  Integer delay, FFT, and fractional delay
	long secStart, fftIndex;
	npy_intp *dLoc, *fLoc, *qLoc;
	dLoc = PyDimMem_NEW(1);
	fLoc = PyDimMem_NEW(2);
	qLoc = PyDimMem_NEW(1);
	for(i=0; i<nFFT; i++) {
		secStart = start + ((long) (2*nChan*((float) i)/Overlap));
		for(j=0; j<2*nChan; j++) {
			dLoc[0] = secStart + j;
			in[j][0] = *(double *) PyArray_GetPtr(data, dLoc);
			in[j][1] = 0.0;
			out[j][0] = 0.0;
			out[j][1] = 0.0;
		}
		
		fftw_execute(p);
		
		fLoc[1] = (npy_intp) i;
		for(j=0; j<(nChan-1); j++) {
			fLoc[0] = (npy_intp) j;
			qLoc[0] = (npy_intp) j;
			fftIndex = j + 1;
			*(double complex *) PyArray_GetPtr(dataF, fLoc) = out[fftIndex][0] + imaginary*out[fftIndex][1];
			*(double complex *) PyArray_GetPtr(dataF, fLoc) *= cexp(-2*imaginary*PI* *(double *) PyArray_GetPtr(fq, qLoc) * frac[j]);
		}
	}
	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
	
	PyDimMem_FREE(dLoc);
	PyDimMem_FREE(fLoc);
	PyDimMem_FREE(qLoc);

	Py_XDECREF(data);
	Py_XDECREF(fq);
	Py_XDECREF(times);

	signalF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalF;
}

PyDoc_STRVAR(FEngineR_doc, "Perform a series of overlaped Fourier transforms on real-valued data.");


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

	// Find out how large the output array needs to be and initialize it
	nFFT = (nSamps - startMax) / 2 / nChan * Overlap - Overlap + 1;
	npy_intp dims[3];
	dims[0] = nStand;
	dims[1] = nChan - 1;
	dims[2] = nFFT;
	dataF = (PyArrayObject*) PyArray_SimpleNew(3, dims, PyArray_CDOUBLE);
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
		#pragma omp parallel default(shared) private(dLoc, fLoc, qLoc, in, j, k, fftIndex)
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
					dLoc[1] = secStart + k;
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

PyDoc_STRVAR(FEngineR2_doc, "Perform a series of overlaped Fourier transforms on real-valued data using OpenMP.");


static PyObject *FEngineC(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signal, *freq, *delay, *signalF;
	PyArrayObject *data, *fq, *times, *dataF;
	int nChan = 64;
	int Overlap = 1;
	double SampleRate = 1.0e5;

	long i, j, nFFT, start;

	static char *kwlist[] = {"signal", "freq", "delay", "LFFT", "Overlap", "SampleRate", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OOO|iid", kwlist, &signal, &freq, &delay, &nChan, &Overlap, &SampleRate)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it useable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signal, PyArray_CDOUBLE, 1, 1);
	fq = (PyArrayObject *) PyArray_ContiguousFromObject(freq, PyArray_DOUBLE, 1, 1);
	times = (PyArrayObject *) PyArray_ContiguousFromObject(delay, PyArray_DOUBLE, 1, 1);
	
	// Check data dimensions
	if(nChan != (fq->dimensions[0]+1)) {
		PyErr_Format(PyExc_RuntimeError, "freq has a different channel count than nChan");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		return NULL;
	}
	
	if(fq->dimensions[0] != times->dimensions[0]) {
		PyErr_Format(PyExc_TypeError, "freq and delays have different channel counts");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		return NULL;
	}

	// Compute the interger sample offset and the fractional sample delay
	npy_intp *tLoc;
	tLoc = PyDimMem_NEW(1);
	tLoc[0] = (npy_intp) fq->dimensions[0] / 2;
	start = (long) round(*(double *) PyArray_GetPtr(times, tLoc) * SampleRate);
	double frac[nChan];
	for(i=0; i<nChan; i++) {
		tLoc[0] = (npy_intp) i;
		frac[i] = *(double *) PyArray_GetPtr(times, tLoc) - start/SampleRate;
	}
	PyDimMem_FREE(tLoc);

	// Find out how large the output array needs to be and initialize it
	nFFT = data->dimensions[0] / nChan * Overlap - Overlap + 1;
	npy_intp dims[2];
	dims[0] = nChan - 1;
	dims[1] = nFFT;
	dataF = (PyArrayObject*) PyArray_SimpleNew(2, dims, PyArray_CDOUBLE);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		Py_XDECREF(fq);
		Py_XDECREF(times);
		return NULL;
	}

	// Create the FFTW plan
	fftw_complex *in, *out;
	fftw_plan p;
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nChan);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nChan);
	p = fftw_plan_dft_1d(nChan, in, out, FFTW_FORWARD, FFTW_MEASURE);

	// Part 1:  Integer delay, FFT, and fractional delay
	long secStart, fftIndex;
	npy_intp *dLoc, *fLoc, *qLoc;
	dLoc = PyDimMem_NEW(1);
	fLoc = PyDimMem_NEW(2);
	qLoc = PyDimMem_NEW(1);
	for(i=0; i<nFFT; i++) {
		secStart = start + ((long) (nChan*((float) i)/Overlap));
		for(j=0; j<nChan; j++) {
			dLoc[0] = secStart + j;
			in[j][0] = creal(*(double complex *) PyArray_GetPtr(data, dLoc));
			in[j][1] = cimag(*(double complex *) PyArray_GetPtr(data, dLoc));
			out[j][0] = 0.0;
			out[j][1] = 0.0;
		}
		
		fftw_execute(p);
		
		fLoc[1] = (npy_intp) i;
		for(j=0; j<(nChan-1); j++) {
			fLoc[0] = (npy_intp) j;
			qLoc[0] = (npy_intp) j;
			fftIndex = ((j+1) + nChan/2) % nChan;
			if( fftIndex < nChan/2 ) {
					fftIndex += 1;
			}
			*(double complex *) PyArray_GetPtr(dataF, fLoc) = (out[fftIndex][0] + imaginary*out[fftIndex][1]);
			*(double complex *) PyArray_GetPtr(dataF, fLoc) *= cexp(-2*imaginary*PI* *(double *) PyArray_GetPtr(fq, qLoc) * frac[j]);
		}
	}
	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
	
	PyDimMem_FREE(dLoc);
	PyDimMem_FREE(fLoc);
	PyDimMem_FREE(qLoc);

	Py_XDECREF(data);
	Py_XDECREF(fq);
	Py_XDECREF(times);

	signalF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);

	return signalF;
}

PyDoc_STRVAR(FEngineC_doc, "Perform a series of overlaped Fourier transforms on complex-valued data.");


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
	double frac[nStand][nChan];
	tLoc = PyDimMem_NEW(2);
	for(i=0; i<nStand; i++) {
		tLoc[0] = (npy_intp) i;
		tLoc[1] = (npy_intp) (nChan / 2);
		start[i] = (long) round(*(double *) PyArray_GetPtr(times, tLoc) * SampleRate);

		for(j=0; j<nChan; j++) {
			tLoc[1] = (npy_intp) j;
			frac[i][j] = *(double *) PyArray_GetPtr(times, tLoc) - (double) start[i]/SampleRate;
		}
	}
	PyDimMem_FREE(tLoc);

	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / nChan * Overlap - Overlap + 1;
	npy_intp dims[3];
	dims[0] = nStand;
	dims[1] = nChan - 1;
	dims[2] = nFFT;
	dataF = (PyArrayObject*) PyArray_SimpleNew(3, dims, PyArray_CDOUBLE);
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
		#pragma omp parallel default(shared) private(dLoc, fLoc, qLoc, in, j, k, fftIndex)
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
					dLoc[1] = secStart + k;
					in[k][0] = creal(*(double complex *) PyArray_GetPtr(data, dLoc));
					in[k][1] = cimag(*(double complex *) PyArray_GetPtr(data, dLoc));
				}	
			
				fftw_execute_dft(p, in, in);
			
				fLoc[2] = (npy_intp) j;
				for(k=0; k<(nChan-1); k++) {
					fLoc[1] = (npy_intp) k;
					qLoc[0] = (npy_intp) k;
					fftIndex = ((k+1) + nChan/2) % nChan;
					if( fftIndex < nChan/2 ) {
						fftIndex += 1;
					}
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

PyDoc_STRVAR(FEngineC2_doc, "Perform a series of overlaped Fourier transforms on complex-valued data using OpenMP.");


static PyObject *XEngine(PyObject *self, PyObject *args) {
	PyObject *signal1, *signal2, *output;
	PyArrayObject *data1, *data2, *vis;
	long nChan, nFFT1, nFFT2, nFFT;

	if(!PyArg_ParseTuple(args, "OO", &signal1, &signal2)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it useable
	data1 = (PyArrayObject *) PyArray_ContiguousFromObject(signal1, PyArray_CDOUBLE, 2, 2);
	data2 = (PyArrayObject *) PyArray_ContiguousFromObject(signal2, PyArray_CDOUBLE, 2, 2);

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
	int dims[1];
	dims[0] = nChan;
	vis = (PyArrayObject*) PyArray_SimpleNew(1, dims, PyArray_CDOUBLE);
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

PyDoc_STRVAR(XEngine_doc, "Perform XMAC on two data streams out of the F engine.\n\
NOTE: This function is *slower* than a pure numpy version of the same function.\n");


static PyObject *XEngine2(PyObject *self, PyObject *args) {
	PyObject *signals, *output;
	PyArrayObject *data, *vis;
	long nStand, nChan, nFFT, nBL;	

	if(!PyArg_ParseTuple(args, "O", &signals)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it useable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, PyArray_CDOUBLE, 3, 3);

	// Get channel count and number of FFTs stored
	nStand = (long) data->dimensions[0];
	nChan = (long) data->dimensions[1];
	nFFT = (long) data->dimensions[2];
	nBL = (nStand+1)*nStand/2;
	
	// Create the output visibility array and fill with zeros
	int dims[2];
	dims[0] = nStand*(nStand+1)/2;
	dims[1] = nChan;
	vis = (PyArrayObject*) PyArray_SimpleNew(2, dims, PyArray_CDOUBLE);
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
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(vLoc, dLoc1, dLoc2, c, f)
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
			
			for(f=0; f<nFFT; f++) {
				dLoc1[2] = (npy_intp) f;
				dLoc2[2] = (npy_intp) f;
			
				for(c=0; c<nChan; c++) {
					vLoc[1] = (npy_intp) c;
					dLoc1[1] = (npy_intp) c;
					dLoc2[1] = (npy_intp) c;
				
					*(double complex *) PyArray_GetPtr(vis, vLoc) += *(double complex *) PyArray_GetPtr(data, dLoc1) * conj(*(double complex *) PyArray_GetPtr(data, dLoc2)) / nFFT;
				}
			}
			
			PyDimMem_FREE(vLoc);
			PyDimMem_FREE(dLoc1);
			PyDimMem_FREE(dLoc2);
		}
	}
	
	Py_XDECREF(data);

	output = Py_BuildValue("O", PyArray_Return(vis));
	Py_XDECREF(vis);

	return output;
}

PyDoc_STRVAR(XEngine2_doc, "Perform all XMACs for a data stream out of the F engine using OpenMP.\n\
Inputs:\n\
  * signals - stand x channel x fft numpy array of frequency-domain signals\n\
Outputs:\n\
  * vis - cross-power spectra for every baseline (include autocorrelations) for\n\
    the input data\n");


static PyMethodDef CorrelatorMethods[] = {
	{"FEngineR", FEngineR, METH_VARARGS | METH_KEYWORDS, FEngineR_doc}, 
	{"FEngineR2", FEngineR2, METH_VARARGS | METH_KEYWORDS, FEngineR2_doc}, 
	{"FEngineC", FEngineC, METH_VARARGS | METH_KEYWORDS, FEngineC_doc}, 
	{"FEngineC2", FEngineC2, METH_VARARGS | METH_KEYWORDS, FEngineC2_doc}, 
	{"XEngine", XEngine, METH_VARARGS, XEngine_doc}, 
	{"XEngine2", XEngine2, METH_VARARGS, XEngine2_doc}, 
	{NULL, NULL, 0, NULL}
};

PyDoc_STRVAR(correlator_doc, "C-based F and X engines for the LWA software FX correlator.  These function\n\
are meant to provide an alternative to the lsl.correlator.fx.correlate function and \n\
provide a much-needed speed boost to cross-correlation.\n\
\n\
The function defined in this module are:\n\
  * FEngineR - F-engine for computing a series of overlapped Fourier transforms with\n\
    delay corrections for a real-valued (TBW) signal from a single stand\n\
  * FEngineR2 - Similar to FEngineR, but works with a collection of stands all at\n\
    once\n\
  * FEngineC - F-engine for computing a series of overlapped Fourier transforms with\n\
    delay corrections for a complex-valued (TBN) signal from a single stand\n\
  * FEngineC2 - Similar to FEngineC, but works with a collection of stands all at\n\
    once\n\
  * XEngine - Cross-multipy and accumulate cross-power spectra using the inputs for\n\
    two stands\n\
  * XEngine2 - Similar to XEngine, but works with a collection of stands all at\n\
    once\n\
See the inidividual functions for more details.");


PyMODINIT_FUNC init_core(void) {
	(void) Py_InitModule3("_core", CorrelatorMethods, correlator_doc);
	import_array();
}

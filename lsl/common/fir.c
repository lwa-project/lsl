#include "Python.h"
#include <math.h>
#include <stdio.h>
#include <cblas.h>
#include <fftw3.h>
#include <stdlib.h>
#include <complex.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

#include "numpy/arrayobject.h"


/*
 applyFIR - Given a pointer to a 16-bit integer data stream, the number of data samples to 
            process, a pointer a set of 16-bit FIR coefficients, the number of taps, and a
            pointer to a 32-bit float output stream, apply the FIR coefficents.
*/

void applyFIR(short int *data, long nSamps, short int *coeff, long nTaps, float *output) {
	long i, j;
	
	memset(output, 0, sizeof(int)*nSamps);
	
	for(i=0; i<nTaps; i++) {
		for(j=0; j<=i; j++) {
			*(output + i) += (float) *(coeff + j) * *(data + i - j);
		}
		*(output + i) /= 32767.0;
	}
	
	for(i=nTaps; i<nSamps; i++) {
		for(j=0; j<nTaps; j++) {
			*(output + i) += (float) *(coeff + j) * *(data + i - j);
		}
		*(output + i) /= 32767.0;
	}
}


/*
 applyFIRDelayed - Like applyFIR but also applies a sample delay before filtering.
*/

void applyFIRDelayed(short int *data, long nSamps, short int *coeff, long nTaps, long sampleDelay, float *output) {
	long i, j;
	
	memset(output, 0, sizeof(int)*nSamps);
	
	for(i=0; i<sampleDelay; i++) {
		*(output + i) = 0;
	}
	
	for(i=sampleDelay; i<(nTaps+sampleDelay); i++) {
		for(j=0; j<=(i-sampleDelay); j++) {
			*(output + i) += (float) *(coeff + j) * *(data + i - sampleDelay - j);
		}
		*(output + i) /= 32767.0;
	}
	
	for(i=(nTaps+sampleDelay); i<nSamps; i++) {
		for(j=0; j<nTaps; j++) {
			*(output + i) += (float) *(coeff + j) * *(data + i - sampleDelay - j);
		}
		*(output + i) /= 32767.0;
	}
}


/*
 integerFIR - Function for taking a numpy.int16 data set and applying a FIR
              filter.  This function is similar to scipy.signal.lfilter but
              is more memory efficient since it doesn't require floats.
*/

static PyObject *integerFIR(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *filter, *output;
	PyArrayObject *data, *coeff, *dataF;
	
	long nSamps, nTaps;
	
	if(!PyArg_ParseTuple(args, "OO", &signals, &filter)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it usable
	data  = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_INT16, 1, 1);
	coeff = (PyArrayObject *) PyArray_ContiguousFromObject(filter, NPY_INT16, 1, 1);
	
	// Get sample and tap counts
	nSamps = (long) data->dimensions[0];
	nTaps  = (long) coeff->dimensions[0];
	
	// Create the output data holders
	npy_intp dims[1];
	dims[0] = (npy_intp) nSamps;
	dataF = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_INT32);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		Py_XDECREF(coeff);
		Py_XDECREF(dataF);
		return NULL;
	}
	PyArray_FILLWBYTE(dataF, 0);
	
	// Go
	short int *a, *b;
	float *c;
	a = (short int *) data->data;
	b = (short int *) coeff->data;
	c = (float *) dataF->data;
	applyFIR(a, nSamps, b, nTaps, c);
	
	
	Py_XDECREF(data);
	Py_XDECREF(coeff);
	
	output = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);
	
	return output;
}

PyDoc_STRVAR(integerFIR_doc, \
"Given a 1-D numpy.int16 array of data values and a numpy.int16 array of FIR\n\
coefficients, apply the coefficients to the data.\n\
\n\
Inputs arguments are:\n\
 * data: 1-D numpy.int16 array of data\n\
 * coeffs: 1-D numpy.int16 array of FIR coefficients\n\
\n\
Outputs:\n\
 * result: 1-D numpy.float32 array of the filtered data\n\
");


/*
 integerFIRDelayed - Function similar to integerFIR but adds in a FIFO delay segment.
*/

static PyObject *integerFIRDelayed(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *filter, *output;
	PyArrayObject *data, *coeff, *dataF;
	
	long nSamps, nTaps, sampleDelay;
	
	if(!PyArg_ParseTuple(args, "OOl", &signals, &filter, &sampleDelay)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it usable
	data  = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_INT16, 1, 1);
	coeff = (PyArrayObject *) PyArray_ContiguousFromObject(filter, NPY_INT16, 1, 1);
	
	// Get sample and tap counts
	nSamps = (long) data->dimensions[0];
	nTaps  = (long) coeff->dimensions[0];
	
	// Create the output data holders
	npy_intp dims[1];
	dims[0] = (npy_intp) nSamps;
	dataF = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		Py_XDECREF(coeff);
		Py_XDECREF(dataF);
		return NULL;
	}
	PyArray_FILLWBYTE(dataF, 0);
	
	// Go
	short int *a, *b;
	float *c;
	a = (short int *) data->data;
	b = (short int *) coeff->data;
	c = (float *) dataF->data;
	applyFIRDelayed(a, nSamps, b, nTaps, sampleDelay, c);
	
	
	Py_XDECREF(data);
	Py_XDECREF(coeff);
	
	output = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);
	
	return output;
}

PyDoc_STRVAR(integerFIRDelayed_doc, \
"Given a 1-D numpy.int16 array of data values, a numpy.int16 array of FIR\n\
coefficients, and a FIFO sample delay, delay the signal and apply the\n\
coefficients to the data.\n\
\n\
Inputs arguments are:\n\
 * data: 1-D numpy.int16 array of data\n\
 * coeffs: 1-D numpy.int16 array of FIR coefficients\n\
 * sampleDelay: interger number of samples to delay the signal (must be >=0)\n\
\n\
Outputs:\n\
 * result: 1-D numpy.float32 array of the delayed and filtered data\n\
");


static PyObject *integerBeamformer(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *filters, *courses, *fines, *gains, *output;
	PyArrayObject *data, *filter, *course, *fine, *gain, *dataFX, *dataFY;
	
	long i, j, k, nStand, nSamps, nFilts, nTaps;
	
	if(!PyArg_ParseTuple(args, "OOOOO", &signals, &filters, &courses, &fines, &gains)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Bring the data into C and make it usable
	data   = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_INT16, 2, 2);
	filter = (PyArrayObject *) PyArray_ContiguousFromObject(filters, NPY_INT16, 3, 3);
	course = (PyArrayObject *) PyArray_ContiguousFromObject(courses, NPY_INT16, 1, 1);
	fine   = (PyArrayObject *) PyArray_ContiguousFromObject(fines,   NPY_INT16, 1, 1);
	gain   = (PyArrayObject *) PyArray_ContiguousFromObject(gains,   NPY_INT16, 2, 2);
	
	// Check data dimensions
	if( data->dimensions[0] != filter->dimensions[0] ) {
		PyErr_Format(PyExc_TypeError, "signals and FIR filters have different input counts");
		Py_XDECREF(data);
		Py_XDECREF(filter);
		Py_XDECREF(course);
		Py_XDECREF(fine);
		Py_XDECREF(gain);
		return NULL;
	}
	if( data->dimensions[0] != course->dimensions[0] ) {
		PyErr_Format(PyExc_TypeError, "signals and course delays have different input counts");
		Py_XDECREF(data);
		Py_XDECREF(filter);
		Py_XDECREF(course);
		Py_XDECREF(fine);
		Py_XDECREF(gain);
		return NULL;
	}
	if( data->dimensions[0] != fine->dimensions[0] ) {
		PyErr_Format(PyExc_TypeError, "signals and find delays have different input counts");
		Py_XDECREF(data);
		Py_XDECREF(filter);
		Py_XDECREF(course);
		Py_XDECREF(fine);
		Py_XDECREF(gain);
		return NULL;
	}
	if( data->dimensions[0]/2 != gain->dimensions[0] ) {
		PyErr_Format(PyExc_TypeError, "signals and gains have different input counts");
		Py_XDECREF(data);
		Py_XDECREF(filter);
		Py_XDECREF(course);
		Py_XDECREF(fine);
		Py_XDECREF(gain);
		return NULL;
	}
	if( gain->dimensions[1] != 4 ) {
		PyErr_Format(PyExc_TypeError, "seconds dimension of gains must be four");
		Py_XDECREF(data);
		Py_XDECREF(filter);
		Py_XDECREF(course);
		Py_XDECREF(fine);
		Py_XDECREF(gain);
		return NULL;
	}
	
	// Get sample and tap counts
	nStand = (long) data->dimensions[0];
	nSamps = (long) data->dimensions[1];
	nFilts = (long) filter->dimensions[1];
	nTaps  = (long) filter->dimensions[2];
	
	long maxCourse = 0;
	short int *c;
	c = (short int *) course->data;
	for(i=0; i<nStand/2; i++) {
		k = 2*i;
		if( *(c + k) > maxCourse ) {
			maxCourse = (long) *(c + k);
		}
		k = 2*i + 1;
		if( *(c + k) > maxCourse ) {
			maxCourse = (long) *(c + k);
		}
	}
	maxCourse += nTaps/2 - 1;
	
	// Create the output data holders
	npy_intp dims[1];
	dims[0] = (npy_intp) (nSamps - maxCourse);
	dataFX = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	if(dataFX == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array for X polarization");
		Py_XDECREF(data);
		Py_XDECREF(filter);
		Py_XDECREF(course);
		Py_XDECREF(fine);
		Py_XDECREF(gain);
		Py_XDECREF(dataFX);
		return NULL;
	}
	PyArray_FILLWBYTE(dataFX, 0);
	
	dataFY = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	if(dataFY == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array for Y polarization");
		Py_XDECREF(data);
		Py_XDECREF(filter);
		Py_XDECREF(course);
		Py_XDECREF(fine);
		Py_XDECREF(gain);
		Py_XDECREF(dataFX);
		Py_XDECREF(dataFY);
		return NULL;
	}
	PyArray_FILLWBYTE(dataFY, 0);
	
	// Go
	short int *d, *f, *w, *g;
	float *x, *y;
	d = (short int *) data->data;
	f = (short int *) filter->data;
	c = (short int *) course->data;
	w = (short int *) fine->data;
	g = (short int *) gain->data;
	x = (int *) dataFX->data;
	y = (int *) dataFY->data;
	
	float *t1, *t2, *tX, *tY;
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(j, k, t1, t2, tX, tY)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(dynamic)
		#endif
		for(i=0; i<nStand/2; i++) {
			/* 
			  Skip stands with all gains of zero.
			*/
			if( *(g + 4*i + 0) == 0 && *(g + 4*i + 1) == 0 && *(g + 4*i + 2) == 0 && *(g + 4*i + 3) == 0 ) {
				continue;
			}
			
			/*
			  Allocate temporary arrays
			    * t1 -> output of integer delay and FIR filter for X pol.
			    * t2 -> output of integer delay and FIR filter for Y pol.
			    * tX -> temporary beam for X pol.
			    * tY -> temporary beam for Y pol.
			*/
			t1 = (int *) malloc(sizeof(float)*nSamps);
			t2 = (int *) malloc(sizeof(float)*nSamps);
			tX = (int *) malloc(sizeof(float)*nSamps);
			tY = (int *) malloc(sizeof(float)*nSamps);
			
			/*
			  Delay + FIR
			*/
			k = 2*i;
			applyFIRDelayed((d+k*nSamps), nSamps, (f + k*nFilts*nTaps + *(w+k)*nTaps), nTaps, *(c+k), t1);
			k = 2*i + 1;
			applyFIRDelayed((d+k*nSamps), nSamps, (f + k*nFilts*nTaps + *(w+k)*nTaps), nTaps, *(c+k), t2);
			
			for(j=0; j<nSamps; j++) {
				*(tX + j) = *(t1 + j) * *(g + 4*i + 0) / 32767.0;
				*(tY + j) = *(t1 + j) * *(g + 4*i + 1) / 32767.0;
				
				*(tX + j) += *(t2 + j) * *(g + 4*i + 2) / 32767.0;
				*(tY + j) += *(t2 + j) * *(g + 4*i + 3) / 32767.0;
			}
			
			/*
			  Add the partial sums to the full beam
			*/
			#ifdef _OPENMP
				#pragma omp critical
			#endif
			{
				for(j=0; j<(nSamps-maxCourse); j++) {
					*(x + j) += *(tX + j + maxCourse);
					*(y + j) += *(tY + j + maxCourse);
				}
			}
			
			/*
			  Cleanup
			*/
			free(t1);
			free(t2);
			free(tX);
			free(tY);
		}
	}
	
	Py_XDECREF(data);
	Py_XDECREF(filter);
	Py_XDECREF(course);
	Py_XDECREF(fine);
	Py_XDECREF(gain);
	
	output = Py_BuildValue("(OO)", PyArray_Return(dataFX), PyArray_Return(dataFY));
	Py_XDECREF(dataFX);
	Py_XDECREF(dataFY);
	
	return output;
}

PyDoc_STRVAR(integerBeamformer_doc, \
"Given a 2-D numpy.int16 array (stands by samples) of data values, 3-D array of FIR\n\
filter coefficients (stands by filters by taps), a 1-D numpy.int16 array of course\n\
(FIFO) delays, a 1-D numpy.int16 array of fine delay FIR filters, and a 2-D array of\n\
gains (stands by [XX, XY, YX, YY]), apply the delays and sum the signals.\n\
\n\
Inputs arguments are:\n\
 * data: 2-D numpy.int16 array of data (stands by samples)\n\
 * coeffs: 3-D numpy.int16 array of FIR coefficients (stands by filters by taps)\n\
 * course: 1-D numpy.int16 array of FIFO delays in samples\n\
 * fine: 1-D numpy.int16 array of which FIR filter to apply for fine delay\n\
 * gain: 2-D numpy.int16 arry of gains (stands by [XX, XY, YX, YY]), where XX is the X\n\
         contribution to the output X pol., XY is the X contribution to the output Y\n\
         pol., YX is the Y contribtion to the output X pol., and YY is the Y contribu-\n\
         tion to the output Y pol.\n\
\n\
Outputs:\n\
 * results: two element tuple (output X, outpuY) of the beamformer sum.  Each element is\n\
            a 1-D numpy.float32 array.\n\
\n\
.. note::\n\
\tThe structure of data is assumed to be that the polarizations are ordered, e.g., the X\n\
\tpolarization of stand 1 is immediately followed by the Y polarization.\n\
");


static PyMethodDef FIRMethods[] = {
	{"integer16",         integerFIR,        METH_VARARGS, integerFIR_doc       }, 
	{"integer16Delayed",  integerFIRDelayed, METH_VARARGS, integerFIRDelayed_doc}, 
	{"integerBeamformer", integerBeamformer, METH_VARARGS, integerBeamformer_doc}, 
	{NULL,                NULL,              0,            NULL                 }
};

PyDoc_STRVAR(FIRMethods_doc, \
"This module contains a collection of function to speed up FIR filtering of TBW\n\
data (represented as numpy.int16 arrays) and the SoftwareDP class.  The funtions\n\
provided in this module are:\n\
 * integer16: Apply a FIR filter to numpy.int16 data,\n\
 * integer16Delayed: Apply a FIFO delay and a FIR filter to numpy.int16 data, and\n\
 * integerBeamformer: Software implementation of the DP beamformer.\n\
");


/*
  Module Setup - Initialization
*/

PyMODINIT_FUNC init_fir(void) {
	PyObject *m;

	// Module definitions and functions
	m = Py_InitModule3("_fir", FIRMethods, FIRMethods_doc);
	import_array();
	
	// Version and revision information
	PyModule_AddObject(m, "__version__", PyString_FromString("0.1"));
	PyModule_AddObject(m, "__revision__", PyString_FromString("$Rev$"));
	
}

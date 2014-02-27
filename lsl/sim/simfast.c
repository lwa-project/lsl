#include "Python.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

#include "numpy/arrayobject.h"

#define PI 3.1415926535898
#define imaginary _Complex_I


static PyObject *FastVis(PyObject *self, PyObject *args) {
	PyObject *antarray, *bls, *output, *temp, *temp2, *temp3;
	PyArrayObject *freq, *ha, *dec, *flux, *uvwF, *visF, *tempA;
	long int i, j, k;
	long int nAnt, nSrc, nFreq, nBL, chanMin, chanMax;
	double lat, pcAz, pcEl, pcHA, pcDec;
	
	chanMin = 0;
	chanMax = -1;
	pcAz = 0.0;
	pcEl = 90.0;
	if(!PyArg_ParseTuple(args, "OOlldd", &antarray, &bls, &chanMin, &chanMax, &pcAz, &pcEl)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	// Bring the data into C and make it usable
	/* Site latitude */
	temp = PyObject_GetAttrString(antarray, "lat");
	lat = PyFloat_AsDouble(temp);
	Py_DECREF(temp);
	
	/* Frequencies (GHz) */
	temp = PyObject_CallMethod(antarray, "get_afreqs", NULL);
	freq = (PyArrayObject *) PyArray_ContiguousFromObject(temp, NPY_DOUBLE, 1, 2);
	Py_DECREF(temp);
	temp = PyArray_Squeeze(freq);
	Py_XDECREF(freq);
	freq = (PyArrayObject *) PyArray_ContiguousFromObject(temp, NPY_DOUBLE, 1, 1);
	Py_DECREF(temp);
	
	/* Source cache and properties */
	/** Cache **/
	temp = PyObject_GetAttrString(antarray, "_cache");
	/** Hour angle **/
	temp2 = PyDict_GetItemString(temp, "s_ha");
	if( temp2 == NULL ) {
		PyErr_Format(PyExc_TypeError, "Cannot find HA array 's_ha' in the simulation cache");
		Py_XDECREF(freq);
		return NULL;
	}
	ha   = (PyArrayObject *) PyArray_ContiguousFromObject(temp2, NPY_DOUBLE, 1, 1);
	/** Declination **/
	temp2 = PyDict_GetItemString(temp, "s_dec");
	if( temp2 == NULL ) {
		PyErr_Format(PyExc_TypeError, "Cannot find dec. array 's_dec' in the simulation cache");
		Py_XDECREF(freq);
		Py_XDECREF(ha);
		return NULL;
	}
	dec  = (PyArrayObject *) PyArray_ContiguousFromObject(temp2, NPY_DOUBLE, 1, 1);
	/** Flux density as a function of frequency **/
	temp2 = PyDict_GetItemString(temp, "jys");
	if( temp2 == NULL ) {
		PyErr_Format(PyExc_TypeError, "Cannot find flux density array 'jys' in the simulation cache");
		Py_XDECREF(freq);
		Py_XDECREF(ha);
		Py_XDECREF(dec);
		return NULL;
	}
	flux = (PyArrayObject *) PyArray_ContiguousFromObject(temp2, NPY_COMPLEX128, 2, 2);
	Py_DECREF(temp);
	
	/* Pointing center */
	pcAz *= PI / 180.0;
	pcEl *= PI / 180.0;
	if( pcEl == PI/2 ) {
		pcEl -= 1e-8;
	}
	/** Conversion to hour angle and declination **/
	pcHA = atan2(sin(pcAz-PI), (cos(pcAz-PI)*sin(lat) + tan(pcEl)*cos(lat)));
	pcDec = asin(sin(lat)*sin(pcEl) - cos(lat)*cos(pcEl)*cos(pcAz-PI));
	
	// Check data dimensions
	if(ha->dimensions[0] != dec->dimensions[0]) {
		PyErr_Format(PyExc_TypeError, "Source hour angle and declination arrays do not contain the same number of elements");
		Py_XDECREF(freq);
		Py_XDECREF(ha);
		Py_XDECREF(dec);
		Py_XDECREF(flux);
		return NULL;
	}
	
	if(flux->dimensions[0] != ha->dimensions[0]) {
		PyErr_Format(PyExc_TypeError, "Source flux dimensions do not agree with number of sources");
		Py_XDECREF(freq);
		Py_XDECREF(ha);
		Py_XDECREF(dec);
		Py_XDECREF(flux);
		return NULL;
	}
	
	if(flux->dimensions[1] != freq->dimensions[0]) {
		PyErr_Format(PyExc_TypeError, "Source flux dimensions do not agree with number of channels");
		Py_XDECREF(freq);
		Py_XDECREF(ha);
		Py_XDECREF(dec);
		Py_XDECREF(flux);
		return NULL;
	}
	
	// Dimensions
	temp = PyObject_CallMethod(antarray, "__len__", NULL);
	nAnt = (long int) PyInt_AsLong(temp);
	temp2 = PyObject_CallMethod(bls, "__len__", NULL);
	nBL = (long int) PyInt_AsLong(temp2);
	nFreq = (long int) freq->dimensions[0];
	nSrc = (long int) ha->dimensions[0];
	Py_DECREF(temp);
	Py_DECREF(temp2);
	if( chanMax < chanMin ) {
		chanMax = nFreq;
	} else {
		chanMax += 1;
	}
	printf("Found nAnt: %li\n      nBL: %li\n      nFreq: %li\n      nSrc: %li\n", nAnt, nBL, nFreq, nSrc);
	printf("Channel Range: %li to %li\n", chanMin, chanMax);
	
	// Find out how large the output array needs to be and initialize it
	npy_intp dims1[3];
	dims1[0] = (npy_intp) nBL;
	dims1[1] = (npy_intp) 3;
	dims1[2] = (npy_intp) nFreq;
	uvwF = (PyArrayObject*) PyArray_SimpleNew(3, dims1, NPY_DOUBLE);
	if(uvwF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create uvw output array");
		Py_XDECREF(freq);
		Py_XDECREF(ha);
		Py_XDECREF(dec);
		Py_XDECREF(flux);
		return NULL;
	}
	PyArray_FILLWBYTE(uvwF, 0);
	
	npy_intp dims2[2];
	dims2[0] = (npy_intp) nBL;
	dims2[1] = (npy_intp) nFreq;
	visF = (PyArrayObject*) PyArray_SimpleNew(2, dims2, NPY_COMPLEX64);
	if(visF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create visibility output array");
		Py_XDECREF(freq);
		Py_XDECREF(ha);
		Py_XDECREF(dec);
		Py_XDECREF(flux);
		Py_XDECREF(uvwF);
		return NULL;
	}
	PyArray_FILLWBYTE(visF, 0);
	
	// Load in the antenna equatorial positions
	double *pos, *t;
	pos = (double *) malloc(nAnt*3*sizeof(double));
	for(i=0; i<nAnt; i++) {
		temp = PyInt_FromLong(i);
		temp2 = PyObject_GetItem(antarray, temp);
		temp3 = PyObject_GetAttrString(temp2, "pos");
		Py_DECREF(temp);
		Py_DECREF(temp2);
		
		tempA = (PyArrayObject *) PyArray_ContiguousFromObject(temp3, NPY_DOUBLE, 1, 1);
		t = (double *) tempA->data;
		
		for(j=0; j<3; j++) {
			*(pos + 3*i + j) = *(t + j);
		}
		
		Py_XDECREF(tempA);
		Py_DECREF(temp3);
	}
	
	// Load in baseline pairs
	int *bll;
	bll = (int *) malloc(nBL*2*sizeof(int));
	for(i=0; i<nBL; i++) {
		temp = PyInt_FromLong(i);
		temp2 = PyObject_GetItem(bls, temp);
		Py_DECREF(temp);
		
		for(j=0; j<2; j++) {
			temp = PyInt_FromLong(j);
			temp3 = PyObject_GetItem(temp2, temp);
			*(bll + 2*i + j) = (int) PyInt_AsLong(temp3);
			Py_DECREF(temp);
			Py_DECREF(temp3);
		}
		
		Py_DECREF(temp2);
	}
	
	// Equatorial to topocentric baseline conversion basis for the phase center
	double pcsinHA, pccosHA, pcsinDec, pccosDec;
	pcsinHA = sin(pcHA);
	pccosHA = cos(pcHA);
	pcsinDec = sin(pcDec);
	pccosDec = cos(pcDec);
	
	// Setup variables for the loop
	int a1, a2;
	double blx, bly, blz, x, y, z, u, v, w;
	double tempHA, tempDec;
	float complex *tempVis;
	double *a, *b, *c, *e;
	double complex *d;
	float complex *f;
	a = (double *) freq->data;
	b = (double *) ha->data;
	c = (double *) dec->data;
	d = (double complex *) flux->data;
	e = (double *) uvwF->data;
	f = (float complex *) visF->data;
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(a1, a2, blx, bly, blz, tempHA, tempDec, tempVis, x, y, z, u, v, w, i, j, k)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(static)
		#endif
		for(i=0; i<nBL; i++) {
			// Antenna indicies for the baseline
			a1 = *(bll + 2*i + 0);
			a2 = *(bll + 2*i + 1);
			
			// Baseline in equatorial coordinates
			blx = *(pos + 3*a1 + 0) - *(pos + 3*a2 + 0);
			bly = *(pos + 3*a1 + 1) - *(pos + 3*a2 + 1);
			blz = *(pos + 3*a1 + 2) - *(pos + 3*a2 + 2);
			
			// Baseline visibility
			tempVis = (float complex *) malloc(nFreq*sizeof(float complex));
			memset(tempVis, 0, nFreq*sizeof(float complex));
			
			for(j=0; j<nSrc; j++) {
				if( cabs(*(d + nFreq*j + nFreq/2)) <= 0 ) {
					continue;
				}
				
				// Source pointing
				tempHA = *(b + j);
				tempDec = *(c + j);
				
				// Baseline to topocentric coordinates - z only
				z =  cos(tempDec)*cos(tempHA)*blx - cos(tempDec)*sin(tempHA)*bly + sin(tempDec)*blz;
				
				for(k=chanMin; k<chanMax; k++) {
					// Compute w
					w = *(a + k) * z;
					
					// Compute the contribution of this source to the baseline visibility (with the conjugation)
					*(tempVis + k) += *(d + nFreq*j + k) * cexp(2*PI*imaginary*w);
				}
			}
			
			// Zenith pointing
			x =  pcsinHA*blx +          pccosHA*bly;
			y = -pcsinDec*pccosHA*blx + pcsinDec*pcsinHA*bly + pccosDec*blz;
			z =  pccosDec*pccosHA*blx - pccosDec*pcsinHA*bly + pcsinDec*blz;
			
			for(k=chanMin; k<chanMax; k++) {
				// Compute u, v, and w for a zenith pointing (hardcoded for LWA1)
				u = *(a + k) * x;
				v = *(a + k) * y;
				w = *(a + k) * z;
				
				// Save
				*(e + i*3*nFreq + 0*nFreq + k) = u;
				*(e + i*3*nFreq + 1*nFreq + k) = v;
				*(e + i*3*nFreq + 2*nFreq + k) = w;
				
				// Phase to zenith
				*(f + i*nFreq + k) = *(tempVis + k) * cexp(-2*PI*imaginary*w);
			}
			
			free(tempVis);
		}
	}
	
	// Cleanup
	free(pos);
	free(bll);
	Py_XDECREF(freq);
	Py_XDECREF(ha);
	Py_XDECREF(dec);
	Py_XDECREF(flux);
	
	output = Py_BuildValue("(OO)", PyArray_Return(uvwF), PyArray_Return(visF));
	Py_XDECREF(uvwF);
	Py_XDECREF(visF);
	
	return output;
}

PyDoc_STRVAR(FastVis_doc, \
"\n\
");


/*
  Module Setup - Function Definitions and Documentation
*/

static PyMethodDef SimMethods[] = {
	{"FastVis", (PyCFunction) FastVis, METH_VARARGS, FastVis_doc}, 
	{NULL,      NULL,                  0,            NULL       }
};

PyDoc_STRVAR(sim_doc, \
"\n\
");


/*
  Module Setup - Initialization
*/

PyMODINIT_FUNC init_simfast(void) {
	PyObject *m;

	// Module definitions and functions
	m = Py_InitModule3("_simfast", SimMethods, sim_doc);
	import_array();
	
	// Version and revision information
	PyModule_AddObject(m, "__version__", PyString_FromString("0.1"));
	PyModule_AddObject(m, "__revision__", PyString_FromString("$Rev$"));
}


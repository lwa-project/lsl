#include "Python.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

#include "numpy/arrayobject.h"

#include "protos.h"

#define PI 3.1415926535898
#define imaginary _Complex_I


static PyObject *FastVis(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *antarray, *bls, *output, *temp, *temp2, *temp3;
	PyArrayObject *freq, *ha, *dec, *flux, *shape, *uvwF, *visF, *tempA;
	long int i, j, k;
	long int resolve, nAnt, nSrc, nFreq, nBL, chanMin, chanMax;
	double lat, pcAz, pcEl, pcHA, pcDec;
	
	chanMin = 0;
	chanMax = -1;
	pcAz = 0.0;
	pcEl = 90.0;
	resolve = 0;
	static char *kwlist[] = {"aa", "bls", "chanMin", "chanMax", "pcAz", "pcEl", "resolve_src", NULL};
	if( !PyArg_ParseTupleAndKeywords(args, kwds, "OOlldd|i", kwlist, &antarray, &bls, &chanMin, &chanMax, &pcAz, &pcEl, &resolve) ) {
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
	if( PyArray_NDIM(freq) == 2 ) {
		// Flatten 2-D arrays since the first dimension is one
		Py_DECREF(temp);
		temp = PyArray_Flatten(freq, NPY_ANYORDER);
		Py_XDECREF(freq);
		freq = (PyArrayObject *) PyArray_ContiguousFromAny(temp, NPY_DOUBLE, 1, 1);
	}
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
	/** Source shape **/
	temp2 = PyDict_GetItemString(temp, "s_shp");
	if( temp2 == NULL ) {
		PyErr_Format(PyExc_TypeError, "Cannot find source shape array 's_shp' in the simulation cache");
		Py_XDECREF(freq);
		Py_XDECREF(ha);
		Py_XDECREF(dec);
		Py_XDECREF(flux);
		return NULL;
	}
	shape = (PyArrayObject *) PyArray_ContiguousFromObject(temp2, NPY_DOUBLE, 2, 2);
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
		Py_XDECREF(shape);
		return NULL;
	}
	
	if(flux->dimensions[0] != ha->dimensions[0]) {
		PyErr_Format(PyExc_TypeError, "Source flux dimensions do not agree with number of sources");
		Py_XDECREF(freq);
		Py_XDECREF(ha);
		Py_XDECREF(dec);
		Py_XDECREF(flux);
		Py_XDECREF(shape);
		return NULL;
	}
	
	if(shape->dimensions[0] != 3) {
		PyErr_Format(PyExc_TypeError, "Source shape dimensions do not agree with number of required parameters");
		Py_XDECREF(freq);
		Py_XDECREF(ha);
		Py_XDECREF(dec);
		Py_XDECREF(flux);
		Py_XDECREF(shape);
		return NULL;
	}
	if(shape->dimensions[1] != ha->dimensions[0]) {
		PyErr_Format(PyExc_TypeError, "Source shape dimensions do not agree with number of sources");
		Py_XDECREF(freq);
		Py_XDECREF(ha);
		Py_XDECREF(dec);
		Py_XDECREF(flux);
		Py_XDECREF(shape);
		return NULL;
	}
	
	if(flux->dimensions[1] != freq->dimensions[0]) {
		PyErr_Format(PyExc_TypeError, "Source flux dimensions do not agree with number of channels");
		Py_XDECREF(freq);
		Py_XDECREF(ha);
		Py_XDECREF(dec);
		Py_XDECREF(flux);
		Py_XDECREF(shape);
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
	/*
	printf("Found nAnt: %li\n      nBL: %li\n      nFreq: %li\n      nSrc: %li\n", nAnt, nBL, nFreq, nSrc);
	printf("Channel Range: %li to %li\n", chanMin, chanMax);
	*/
	
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
		Py_XDECREF(shape);
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
		Py_XDECREF(shape);
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
	double tempHA, tempDec, tempA0, tempA1, tempTheta, tempX;
	float complex *tempVis;
	double *a, *b, *c, *e, *g;
	double complex *d;
	float complex *f;
	a = (double *) freq->data;
	b = (double *) ha->data;
	c = (double *) dec->data;
	d = (double complex *) flux->data;
	e = (double *) uvwF->data;
	f = (float complex *) visF->data;
	g = (double *) shape->data;
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(a1, a2, blx, bly, blz, tempHA, tempDec, tempA0, tempA1, tempTheta, tempX, tempVis, x, y, z, u, v, w, i, j, k)
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
				// Source pointing
				tempHA = *(b + j);
				tempDec = *(c + j);
				
				// Shape
				tempA0 = *(g + 0*nSrc + j);
				tempA1 = *(g + 1*nSrc + j);
				tempTheta = *(g + 2*nSrc + j);
				
				// Baseline to topocentric coordinates
				x =  sin(tempHA)*blx +              cos(tempHA)*bly;
				y = -sin(tempDec)*cos(tempHA)*blx + sin(tempDec)*sin(tempHA)*bly + cos(tempDec)*blz;
				z =  cos(tempDec)*cos(tempHA)*blx - cos(tempDec)*sin(tempHA)*bly + sin(tempDec)*blz;
				
				for(k=chanMin; k<chanMax; k++) {
					// Compute w
					u = *(a + k) * x;
					v = *(a + k) * y;
					w = *(a + k) * z;
					
					// Correction for the source shape
					if( resolve && tempA0 != 0.0 && tempA1 != 0.0 ) {
						tempX  = tempA0*(u*cos(tempTheta) - v*sin(tempTheta)) * tempA0*(u*cos(tempTheta) - v*sin(tempTheta));
						tempX += tempA1*(u*sin(tempTheta) + v*cos(tempTheta)) * tempA1*(u*sin(tempTheta) + v*cos(tempTheta));
						tempX = 2.0*PI * sqrt(tempX);
						
						if( tempX != 0.0 ) {
							tempX = 2.0 * j1(tempX)/tempX;
						} else {
							tempX = 1.0;
						}
					} else {
						tempX = 1.0;
					}
					
					// Compute the contribution of this source to the baseline visibility (with the conjugation)
					*(tempVis + k) += tempX * *(d + nFreq*j + k) * cexp(2*PI*imaginary*w);
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
	Py_XDECREF(shape);
	
	output = Py_BuildValue("(OO)", PyArray_Return(uvwF), PyArray_Return(visF));
	Py_XDECREF(uvwF);
	Py_XDECREF(visF);
	
	return output;
}

PyDoc_STRVAR(FastVis_doc, \
"Fast visibility simulation package based on the AIPY amp.AntennaArray.sum()\n\
method.  This function differs from sim() in the sense that it does not\n\
support the ionospheric refraction terms.  However, it is implemtned using\n\
OpenMP and should be significantly faster for larger simulations.\n\
\n\
Inputs arguements are:\n\
  * aa: AntennaArray instances generated by lsl.sim.vis.buildAntennaArray()\n\
  * bls: A list of baseline pairs to compute visibilities for\n\
  * chanMin: The first frequency channel to calculate the visibilities for\n\
  * chanMax: The last frequency channel to calculate the visbilities for\n\
  * pcAz: The azimuth of the phase center in degrees\n\
  * pcEl: The elevation of the phase center in degrees\n\
\n\
Input keywords are:\n\
  * resolve_src: Boolean of whether or not source sizes should be used in\n\
    the simulation.  If this is set to False, all sources are treated as\n\
    unresolved.\n\
\n\
Outputs:\n\
  * uvw: A 3-D numpy.float64 array of uvw coordinates (baselines by (u,v,w)\n\
    by channels)\n\
  * vis: A 2-D numpy.complex64 array of complex visibilities (baselines by\n\
    channels)\n\
");


/*
  Module Setup - Function Definitions and Documentation
*/

static PyMethodDef SimMethods[] = {
	{"FastVis", (PyCFunction) FastVis, METH_VARARGS|METH_KEYWORDS, FastVis_doc}, 
	{NULL,      NULL,                  0,                          NULL       }
};

PyDoc_STRVAR(sim_doc, \
"C-based visibility simulation engine.  These functions are meant to provide\n\
an alternative to the AIPY simulation methods and a much-needed speed boost\n\
to simulation-heavy tasks.\n\
\n\
The functions defined in this modulae are:\n\
  * FastVis - Compute uvw coordinates and visibilities for the provided array\n\
    and baseline list.\n\
\n\
See the inidividual functions for more details.\n\
\n\
.. versionadded:: 1.0.1");


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
	PyModule_AddObject(m, "__revision__", PyString_FromString("$Rev: 1639 $"));
}


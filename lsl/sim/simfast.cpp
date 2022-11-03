#include "Python.h"
#include <cmath>
#include <complex>

#ifdef _OPENMP
    #include <omp.h>
    
    // OpenMP scheduling method
    #ifndef OMP_SCHEDULER
    #define OMP_SCHEDULER dynamic
    #endif
#endif

#include "numpy/arrayobject.h"
#include "numpy/npy_math.h"

#include "../common/py3_compat.h"

#include "protos.h"


/*
  Complex types
*/

typedef std::complex<float> Complex32;
typedef std::complex<double> Complex64;


/*
  Macro for 2*pi
*/

#define TPI (2*NPY_PI*Complex64(0,1))


template<typename FluxType, typename OutType>
void compute_visibility(long nBL,
                        long nChan,
                        long nSrc,
                        long chanMin,
                        long chanMax,
                        long resolve,
                        double phs_ha,
                        double phs_dec,
                        int const* bl_pair,
                        double const* ant_pos,
                        double const* freq,
                        double const* src_ha,
                        double const* src_dec,
                        FluxType const* src_flux,
                        double const* src_shape,
                        double* uvw,
                        OutType* vis) {
    // Setup
    long int i, j, k;
    int a1, a2;
    double blx, bly, blz;
    Complex64 *tempVis;
    double tempHA, tempDec;
    double tempA0, tempA1, tempTheta, tempX;
    double x, y, z;
    double u, v, w;
    
    Py_BEGIN_ALLOW_THREADS
    
    // Equatorial to topocentric baseline conversion basis for the phase center
    double pcsinHA, pccosHA, pcsinDec, pccosDec;
    pcsinHA = sin(phs_ha);
    pccosHA = cos(phs_ha);
    pcsinDec = sin(phs_dec);
    pccosDec = cos(phs_dec);
    
    #ifdef _OPENMP
        #pragma omp parallel default(shared) private(a1, a2, blx, bly, blz, tempHA, tempDec, tempA0, tempA1, tempTheta, tempX, tempVis, x, y, z, u, v, w, i, j, k)
    #endif
    {
        tempVis = (Complex64 *) malloc(nChan*sizeof(Complex64));
        
        #ifdef _OPENMP
            #pragma omp for schedule(OMP_SCHEDULER)
        #endif
        for(i=0; i<nBL; i++) {
            // Antenna indicies for the baseline
            a1 = *(bl_pair + 2*i + 0);
            a2 = *(bl_pair + 2*i + 1);
            
            // Baseline in equatorial coordinates
            blx = *(ant_pos + 3*a1 + 0) - *(ant_pos + 3*a2 + 0);
            bly = *(ant_pos + 3*a1 + 1) - *(ant_pos + 3*a2 + 1);
            blz = *(ant_pos + 3*a1 + 2) - *(ant_pos + 3*a2 + 2);
            
            // Baseline visibility
            memset(tempVis, 0, nChan*sizeof(Complex64));
            
            for(j=0; j<nSrc; j++) {
                // Source pointing
                tempHA = *(src_ha + j);
                tempDec = *(src_dec + j);
                
                // Shape
                tempA0 = *(src_shape + 0*nSrc + j);
                tempA1 = *(src_shape + 1*nSrc + j);
                tempTheta = *(src_shape + 2*nSrc + j);
                
                // Baseline to topocentric coordinates
                x =  sin(tempHA)*blx +              cos(tempHA)*bly;
                y = -sin(tempDec)*cos(tempHA)*blx + sin(tempDec)*sin(tempHA)*bly + cos(tempDec)*blz;
                z =  cos(tempDec)*cos(tempHA)*blx - cos(tempDec)*sin(tempHA)*bly + sin(tempDec)*blz;
                
                for(k=chanMin; k<chanMax; k++) {
                    // Compute w
                    u = *(freq + k) * x;
                    v = *(freq + k) * y;
                    w = *(freq + k) * z;
                    
                    // Correction for the source shape
                    if( resolve && tempA0 != 0.0 && tempA1 != 0.0 ) {
                        tempX  = tempA0*(u*cos(tempTheta) - v*sin(tempTheta)) * tempA0*(u*cos(tempTheta) - v*sin(tempTheta));
                        tempX += tempA1*(u*sin(tempTheta) + v*cos(tempTheta)) * tempA1*(u*sin(tempTheta) + v*cos(tempTheta));
                        tempX = 2.0*NPY_PI * sqrt(tempX);
                        
                        if( tempX != 0.0 ) {
                            tempX = 2.0 * j1(tempX)/tempX;
                        } else {
                            tempX = 1.0;
                        }
                    } else {
                        tempX = 1.0;
                    }
                    
                    // Compute the contribution of this source to the baseline visibility (with the conjugation)
                    *(tempVis + k) += tempX * ((Complex64) *(src_flux + nChan*j + k)) * exp(TPI*w);
                }
            }
            
            // Zenith pointing
            x =  pcsinHA*blx +          pccosHA*bly;
            y = -pcsinDec*pccosHA*blx + pcsinDec*pcsinHA*bly + pccosDec*blz;
            z =  pccosDec*pccosHA*blx - pccosDec*pcsinHA*bly + pcsinDec*blz;
            
            for(k=chanMin; k<chanMax; k++) {
                // Compute u, v, and w for a zenith pointing (hardcoded for LWA1)
                u = *(freq + k) * x;
                v = *(freq + k) * y;
                w = *(freq + k) * z;
                
                // Save
                *(uvw + i*3*nChan + 0*nChan + k) = u;
                *(uvw + i*3*nChan + 1*nChan + k) = v;
                *(uvw + i*3*nChan + 2*nChan + k) = w;
                
                // Phase to zenith
                *(vis + i*nChan + k) = ((OutType) (*(tempVis + k) * exp(-TPI*w)));
            }
        }
        
        free(tempVis);
    }
    
    Py_END_ALLOW_THREADS
}


static PyObject *FastVis(PyObject *self, PyObject *args, PyObject *kwds) {
    PyObject *antarray, *bls, *output, *temp, *temp2, *temp3;
    PyArrayObject *freq=NULL, *ha=NULL, *dec=NULL, *flux=NULL, *shape=NULL, *uvwF=NULL, *visF=NULL, *tempA=NULL;
    long int i, j;
    long int resolve, nAnt, nSrc, nFreq, nBL, chanMin, chanMax;
    double lat, pcAz, pcEl, pcHA, pcDec;
    
    chanMin = 0;
    chanMax = -1;
    pcAz = 0.0;
    pcEl = 90.0;
    resolve = 0;
    char const* kwlist[] = {"aa", "bls", "chan_min", "chan_max", "pc_az", "pc_el", "resolve_src", NULL};
    if( !PyArg_ParseTupleAndKeywords(args, kwds, "OOlldd|i", const_cast<char **>(kwlist), &antarray, &bls, &chanMin, &chanMax, &pcAz, &pcEl, &resolve) ) {
        PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
        goto fail;
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
        goto fail;
    }
    ha   = (PyArrayObject *) PyArray_ContiguousFromObject(temp2, NPY_DOUBLE, 1, 1);
    /** Declination **/
    temp2 = PyDict_GetItemString(temp, "s_dec");
    if( temp2 == NULL ) {
        PyErr_Format(PyExc_TypeError, "Cannot find dec. array 's_dec' in the simulation cache");
        goto fail;
    }
    dec  = (PyArrayObject *) PyArray_ContiguousFromObject(temp2, NPY_DOUBLE, 1, 1);
    /** Flux density as a function of frequency **/
    temp2 = PyDict_GetItemString(temp, "jys");
    if( temp2 == NULL ) {
        PyErr_Format(PyExc_TypeError, "Cannot find flux density array 'jys' in the simulation cache");
        goto fail;
    }
    flux = (PyArrayObject *) PyArray_ContiguousFromObject(temp2,
                                                          PyArray_TYPE((PyArrayObject *) temp2),
                                                          2, 2);
    /** Source shape **/
    temp2 = PyDict_GetItemString(temp, "s_shp");
    if( temp2 == NULL ) {
        PyErr_Format(PyExc_TypeError, "Cannot find source shape array 's_shp' in the simulation cache");
        goto fail;
    }
    shape = (PyArrayObject *) PyArray_ContiguousFromObject(temp2, NPY_DOUBLE, 2, 2);
    Py_DECREF(temp);
    
    /* Pointing center */
    pcAz *= NPY_PI / 180.0;
    pcEl *= NPY_PI / 180.0;
    if( pcEl == NPY_PI/2 ) {
        pcEl -= 1e-8;
    }
    /** Conversion to hour angle and declination **/
    pcHA = atan2(sin(pcAz-NPY_PI), (cos(pcAz-NPY_PI)*sin(lat) + tan(pcEl)*cos(lat)));
    pcDec = asin(sin(lat)*sin(pcEl) - cos(lat)*cos(pcEl)*cos(pcAz-NPY_PI));
    
    // Check data dimensions
    if(PyArray_DIM(ha, 0) != PyArray_DIM(dec, 0)) {
        PyErr_Format(PyExc_RuntimeError, "Source hour angle and declination arrays do not contain the same number of elements");
        goto fail;
    }
    
    if(PyArray_DIM(flux, 0) != PyArray_DIM(ha, 0)) {
        PyErr_Format(PyExc_RuntimeError, "Source flux dimensions do not agree with number of sources");
        goto fail;
    }
    
    if(PyArray_DIM(shape, 0) != 3) {
        PyErr_Format(PyExc_RuntimeError, "Source shape dimensions do not agree with number of required parameters");
        goto fail;
    }
    if(PyArray_DIM(shape, 1) != PyArray_DIM(ha, 0)) {
        PyErr_Format(PyExc_RuntimeError, "Source shape dimensions do not agree with number of sources");
        goto fail;
    }
    
    if(PyArray_DIM(flux, 1) != PyArray_DIM(freq, 0)) {
        PyErr_Format(PyExc_RuntimeError, "Source flux dimensions do not agree with number of channels");
        goto fail;
    }
    
    // Dimensions
    temp = PyObject_CallMethod(antarray, "__len__", NULL);
    nAnt = (long int) PyInt_AsLong(temp);
    temp2 = PyObject_CallMethod(bls, "__len__", NULL);
    nBL = (long int) PyInt_AsLong(temp2);
    nFreq = (long int) PyArray_DIM(freq, 0);
    nSrc = (long int) PyArray_DIM(ha, 0);
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
    uvwF = (PyArrayObject*) PyArray_ZEROS(3, dims1, NPY_DOUBLE, 0);
    if(uvwF == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create uvw output array");
        goto fail;
    }
    
    npy_intp dims2[2];
    dims2[0] = (npy_intp) nBL;
    dims2[1] = (npy_intp) nFreq;
    visF = (PyArrayObject*) PyArray_ZEROS(2, dims2, NPY_COMPLEX64, 0);
    if(visF == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create visibility output array");
        goto fail;
    }
    
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
        t = (double *) PyArray_DATA(tempA);
        
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
    
#define LAUNCH_COMPUTE_VIS(IterType) \
        compute_visibility<IterType>(nBL, nFreq, nSrc, chanMin, chanMax, resolve, pcHA, pcDec, \
                                     (int *) bll, \
                                     (double *) pos, \
                                     (double *) PyArray_DATA(freq), \
                                     (double *) PyArray_DATA(ha), \
                                     (double *) PyArray_DATA(dec), \
                                     (IterType *) PyArray_DATA(flux), \
                                     (double *) PyArray_DATA(shape), \
                                     (double *) PyArray_DATA(uvwF), \
                                     (Complex32 *) PyArray_DATA(visF));
    
    switch( PyArray_TYPE(flux) ){
        case( NPY_FLOAT32    ): LAUNCH_COMPUTE_VIS(float);     break;
        case( NPY_FLOAT64    ): LAUNCH_COMPUTE_VIS(double);    break;
        case( NPY_COMPLEX64  ): LAUNCH_COMPUTE_VIS(Complex32); break;
        case( NPY_COMPLEX128 ): LAUNCH_COMPUTE_VIS(Complex64); break;
        default: PyErr_Format(PyExc_RuntimeError, "Unsupport input data type"); goto fail;
    }
        
#undef LAUNCH_COMPUTE_VIS
        
    // Cleanup
    free(pos);
    free(bll);
    
    output = Py_BuildValue("(OO)", PyArray_Return(uvwF), PyArray_Return(visF));
    
    Py_XDECREF(freq);
    Py_XDECREF(ha);
    Py_XDECREF(dec);
    Py_XDECREF(flux);
    Py_XDECREF(shape);
    Py_XDECREF(uvwF);
    Py_XDECREF(visF);
    
    return output;
    
fail:
    Py_XDECREF(freq);
    Py_XDECREF(ha);
    Py_XDECREF(dec);
    Py_XDECREF(flux);
    Py_XDECREF(shape);
    Py_XDECREF(uvwF);
    Py_XDECREF(visF);
    
    return NULL;
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
 * chan_min: The first frequency channel to calculate the visibilities for\n\
 * chan_max: The last frequency channel to calculate the visbilities for\n\
 * pc_az: The azimuth of the phase center in degrees\n\
 * pc_el: The elevation of the phase center in degrees\n\
\n\
Input keywords are:\n\
* resolve_src: Boolean of whether or not source sizes should be used in\n\
               the simulation.  If this is set to False, all sources are\n\
               treated as unresolved.\n\
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

MOD_INIT(_simfast) {
    PyObject *m, *all;
    
    Py_Initialize();
    
    // Module definitions and functions
    MOD_DEF(m, "_simfast", SimMethods, sim_doc);
    if( m == NULL ) {
        return MOD_ERROR_VAL;
    }
    import_array();
    
    // Version and revision information
    PyModule_AddObject(m, "__version__", PyString_FromString("0.1"));
    
    // Function listings
    all = PyList_New(0);
    PyList_Append(all, PyString_FromString("FastVis"));
    PyModule_AddObject(m, "__all__", all);
    
    return MOD_SUCCESS_VAL(m);
}

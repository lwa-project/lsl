#include "Python.h"
#include <cmath>
#include <complex>
#include <fftw3.h>

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
#include "../correlator/common.h"


// Maximum number of w-planes to project
#define MAX_W_PLANES 1024

// Gridding convolution size on a side
#define GRID_KERNEL_SIZE 7


double signed_sqrt(double data) {
    if( data >= 0 ){
        return sqrt(data);
    } else {
        return -sqrt(-data);
    }
}


double gridding_kernel_point(long i, 
                             long j, 
                             double ci, 
                             double cj) {
    double v;
    
    v = 1.0 / 2.0 / NPY_PI / 0.5 / 0.5;
    v *= exp(-(i-ci)*(i-ci)/2.0/0.5/0.5 + -(j-cj)*(j-cj)/2.0/0.5/0.5);
    
    return v;
}


void w_projection_kernel(long nPixSide, 
                         double uvRes, 
                         double w, 
                         Complex32* kernel) {
    long i,j,l,m;
    double cL, cM;
    
    memset(kernel, 0, sizeof(Complex32)*nPixSide*nPixSide);
    
    for(i=0; i<nPixSide; i++) {
        m = i;
        if( m > nPixSide/2 ) {
            m = m - nPixSide;
        }
        cM = ((double) m) / nPixSide / uvRes;
        
        for(j=0; j<nPixSide; j++) {
            l = j;
            if( l > nPixSide/2 ) {
                l = nPixSide - l;
            }
            cL = ((double) l) / nPixSide / uvRes;
            
            if( cL*cL + cM*cM < 1.0 ) {
                *(kernel + nPixSide*i + j) = exp(-TPI*w*(sqrt(1 - cL*cL - cM*cM)-1));
            }
        }
    }
}


int compute_planes(int nVis,
                   double wRes,
                   double const* w,
                   long* planeStart,
                   long* planeStop) {
    int nPlanes = 0;
    long i;
    double csw, psw;
    
    *(planeStart + nPlanes) = 0;
    psw = signed_sqrt( *(w+0) );
    for(i=1; i<nVis; i++) {
        csw = signed_sqrt( *(w+i) );
        if( csw >= psw + wRes) {
            *(planeStop + nPlanes) = i;
            nPlanes += 1;
            *(planeStart + nPlanes) = i;
            psw = csw;
        }
    }
    if( *(planeStop + nPlanes - 1) < nVis) {
        *(planeStop + nPlanes) = nVis;
        nPlanes++;
    }
    
    /*
    for(i=0; i<nPlanes; i++) {
        printf("Plane %li: %li to %li (%.1f%%)\n", i, *(planeStart + i), *(planeStop + i), \
            ((float) *(planeStop + i) - *(planeStart + i)) / nVis * 100.0);
    }
    */
    
    return nPlanes;
}


void compute_gridding(long nVis,
                      long nPixSide,
                      double uvRes,
                      double wRes,
                      double const* u,
                      double const* v,
                      double const* w,
                      Complex32 const* vis,
                      Complex32 const* wgt,
                      Complex32* uv,
                      Complex32* bm) {
    // Setup
    long i, j, l, m;
    
    Py_BEGIN_ALLOW_THREADS
    
    // Figure out the w-planes to use
    int nPlanes = 0;
    long *planeStart, *planeStop;
    planeStart = (long *) malloc(MAX_W_PLANES*sizeof(long));
    planeStop = (long *) malloc(MAX_W_PLANES*sizeof(long));
    nPlanes = compute_planes(nVis, wRes, w, planeStart, planeStop);
    
    long secStart, secStop;
    double avgW, ci, cj, temp, temp2;
    long pi, pj;
    Complex32 *suv, *sbm, *kern;
    static float norm = (float) 1.0 / (nPixSide * nPixSide * nPixSide * nPixSide);
    
    // FFT setup
    Complex32* inP;
    inP = (Complex32*) fftwf_malloc(sizeof(Complex32) * nPixSide*nPixSide);
    fftwf_plan pF, pR;
    pF = fftwf_plan_dft_2d(nPixSide, nPixSide, \
                           reinterpret_cast<fftwf_complex*>(inP), \
                           reinterpret_cast<fftwf_complex*>(inP), \
                           FFTW_FORWARD, FFTW_ESTIMATE);
    pR = fftwf_plan_dft_2d(nPixSide, nPixSide, \
                           reinterpret_cast<fftwf_complex*>(inP), \
                           reinterpret_cast<fftwf_complex*>(inP), \
                           FFTW_BACKWARD, FFTW_ESTIMATE);
    
    // Go!
    #ifdef _OPENMP
        #pragma omp parallel default(shared) private(suv, sbm, kern, i, j, l, m, secStart, secStop, avgW, ci, cj, pi, pj, temp, temp2)
    #endif
    {
        // Initialize the sub-grids and the w projection kernel
        suv  = (Complex32*) fftwf_malloc(sizeof(Complex32) * nPixSide*nPixSide);
        sbm  = (Complex32*) fftwf_malloc(sizeof(Complex32) * nPixSide*nPixSide);
        kern = (Complex32*) fftwf_malloc(sizeof(Complex32) * nPixSide*nPixSide);
        
        #ifdef _OPENMP
            #pragma omp for schedule(OMP_SCHEDULER)
        #endif
        for(j=0; j<nPlanes; j++) {
            // Zero out the sub-grids
            memset(suv, 0, sizeof(Complex32)*nPixSide*nPixSide);
            memset(sbm, 0, sizeof(Complex32)*nPixSide*nPixSide);
            
            // Extract the plane index boundaries
            secStart = *(planeStart + j);
            secStop  = *(planeStop  + j);
            
            // Grid and determine the average value for w for this plane
            avgW = 0.0;
            for(i=secStart; i<secStop; i++) {
                avgW += *(w + i);
                
                // Pixel in the uv plane this corresponds to
                ci = *(v + i) / -uvRes;
                if( ci < 0.0 ) {
                    ci += nPixSide;
                }
                cj = *(u + i) / uvRes;
                if( cj < 0.0 ) {
                    cj += nPixSide;
                }
                
                for(m=0; m<GRID_KERNEL_SIZE; m++) {
                    pi = (long) ci + m - GRID_KERNEL_SIZE/2;
                    temp = 1.0 / 2.0 / NPY_PI / 0.5 / 0.5;
                    temp *= exp(-(pi-ci)*(pi-ci)/2.0/0.5/0.5);
                    
                    pi %= nPixSide;
                    if( pi < 0 ) {
                        pi += nPixSide;
                    }
                    
                    for(l=0; l<GRID_KERNEL_SIZE; l++) {
                        //pi = (long) ci + m - GRID_KERNEL_SIZE/2;
                        pj = (long) cj + l - GRID_KERNEL_SIZE/2;
                        
                        temp2 = temp * exp(-(pj-cj)*(pj-cj)/2.0/0.5/0.5);
                        
                        pj %= nPixSide;
                        if( pj < 0 ) {
                            pj += nPixSide;
                        }
                        
                        *(suv + nPixSide*pi + pj) += *(vis + i) * (float) temp2;
                        *(sbm + nPixSide*pi + pj) += *(wgt + i) * (float) temp2;
                    }
                }
            }
            avgW /= (double) (secStop - secStart);
            
            // Build the w projection kernel
            w_projection_kernel(nPixSide, uvRes, avgW, kern);
            
            // Project
            fftwf_execute_dft(pF, \
                              reinterpret_cast<fftwf_complex*>(suv), \
                              reinterpret_cast<fftwf_complex*>(suv));
            fftwf_execute_dft(pF, \
                              reinterpret_cast<fftwf_complex*>(sbm), \
                              reinterpret_cast<fftwf_complex*>(sbm));
            for(i=0; i<nPixSide*nPixSide; i++) {
                *(suv + i) *= *(kern + i) * norm;
                *(sbm + i) *= *(kern + i) * norm;
            }
            fftwf_execute_dft(pR, \
                              reinterpret_cast<fftwf_complex*>(suv), \
                              reinterpret_cast<fftwf_complex*>(suv));
            fftwf_execute_dft(pR, \
                              reinterpret_cast<fftwf_complex*>(sbm), \
                              reinterpret_cast<fftwf_complex*>(sbm));
            
            #ifdef _OPENMP
            #pragma omp critical
            #endif
            {
                for(i=0; i<nPixSide*nPixSide; i++) {
                    *(uv + i) += *(suv+i);
                    *(bm + i) += *(suv+i);
                }
            }
        }
        
        fftwf_free(suv);
        fftwf_free(sbm);
        fftwf_free(kern);
    }
    
    // Cleanup
    fftwf_destroy_plan(pF);
    fftwf_destroy_plan(pR);
    fftwf_free(inP);
    
    free(planeStart);
    free(planeStop);
    
    Py_END_ALLOW_THREADS
}


static PyObject *WProjection(PyObject *self, PyObject *args, PyObject *kwds) {
    PyObject *uVec, *vVec, *wVec, *visVec, *wgtVec, *output;
    PyArrayObject *uu=NULL, *vv=NULL, *ww=NULL, *vd=NULL, *wd=NULL, *uvPlane=NULL, *bmPlane=NULL;
    long uvSize = 80;
    double uvRes = 0.5;
    double wRes = 0.1;
    
    char const* kwlist[] = {"u", "v", "w", "data", "wgt", "uvSize", "uvRes", "wRes", NULL};
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "OOOOO|ldd", const_cast<char **>(kwlist), &uVec, &vVec, &wVec, &visVec, &wgtVec, &uvSize, &uvRes, &wRes)) {
        PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
        goto fail;
    }
    
    // Bring the data into C and make it usable
    uu = (PyArrayObject *) PyArray_ContiguousFromObject(uVec, NPY_FLOAT64, 1, 1);
    vv = (PyArrayObject *) PyArray_ContiguousFromObject(vVec, NPY_FLOAT64, 1, 1);
    ww = (PyArrayObject *) PyArray_ContiguousFromObject(wVec, NPY_FLOAT64, 1, 1);
    vd = (PyArrayObject *) PyArray_ContiguousFromObject(visVec, NPY_COMPLEX64, 1, 1);
    wd = (PyArrayObject *) PyArray_ContiguousFromObject(wgtVec, NPY_COMPLEX64, 1, 1);
    if( uu == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input u array to 1-D float64");
        goto fail;
    }
    if( vv == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input v array to 1-D float64");
        goto fail;
    }
    if( ww == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input w array to 1-D float64");
        goto fail;
    }
    if( vd == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input data array to 1-D complex64");
        goto fail;
    }
    if( wd == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input wgt array to 1-D complex64");
        goto fail;
    }
    
    // Figure out how much data we have to work with
    long nVis;
    nVis = (long) PyArray_DIM(uu, 0);
    
    // Compute the size of the uv plane
    long nPixSide;
    nPixSide = round(uvSize / uvRes);
    
    // Create the output arrays
    npy_intp dims[2];
    dims[0] = (npy_intp) nPixSide;
    dims[1] = (npy_intp) nPixSide;
    /* uv plane */
    uvPlane = (PyArrayObject *) PyArray_SimpleNew(2, dims, NPY_COMPLEX64);
    if( uvPlane == NULL ) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array - uv data");
        goto fail;
    }
    PyArray_FILLWBYTE(uvPlane, 0);
    /* beam */
    bmPlane = (PyArrayObject *) PyArray_SimpleNew(2, dims, NPY_COMPLEX64);
    if( bmPlane == NULL ) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array - beam data");
        goto fail;
    }
    PyArray_FILLWBYTE(bmPlane, 0);
    
    // Get pointers to the data we need
    double *u, *v, *w;
    Complex32 *vis, *wgt, *uv, *bm;
    u = (double *) PyArray_DATA(uu);
    v = (double *) PyArray_DATA(vv);
    w = (double *) PyArray_DATA(ww);
    vis = (Complex32 *) PyArray_DATA(vd);
    wgt = (Complex32 *) PyArray_DATA(wd);
    uv = (Complex32 *) PyArray_DATA(uvPlane);
    bm = (Complex32 *) PyArray_DATA(bmPlane);
    
    // Grid
    compute_gridding(nVis, nPixSide, uvRes, wRes, u, v, w, vis, wgt, uv, bm);
    
    Py_XDECREF(uu);
    Py_XDECREF(vv);
    Py_XDECREF(ww);
    
    output = Py_BuildValue("(OO)", PyArray_Return(uvPlane), PyArray_Return(bmPlane));
    Py_XDECREF(uvPlane);
    Py_XDECREF(bmPlane);
    
    return output;
    
fail:
    Py_XDECREF(uu);
    Py_XDECREF(vv);
    Py_XDECREF(ww);
    Py_XDECREF(vd);
    Py_XDECREF(wd);
    Py_XDECREF(uvPlane);
    Py_XDECREF(bmPlane);
    
    return NULL;
}

PyDoc_STRVAR(WProjection_doc, \
"w-projection gridder for uv data based on the aipy.img.ImgW class and 'Wide-\n\
field Imaging Problems in Radio Astronomy' (Cornwell et al. 2005).\n\
\n\
Input arguments are:\n\
 * u: 1-D numpy.float64 array of u coordinates\n\
 * v: 1-D numpy.float64 array of v coordinates\n\
 * w: 1-D numpy.float64 array of w coordinates\n\
 * data: 1-D numpy.complex64 array of visibility data\n\
 * wgt: 1-D numpy.complex64 array of visibility weight data\n\
\n\
Input keywords are:\n\
 * uvSize: Basis size of the uv plane\n\
 * uvRes: Resolution of the uv plane\n\
 * wRes: Resolution in w for projection\n\
\n\
Outputs are:\n\
 * uvPlane: 2-D numpy.complex64 of the gridded and projected uv plane\n\
 * bmPlane: 2-D numpy.complex64 of the gridded and projected synthesized\n\
            beam\n\
\n\
.. note::\n\
     All of the input arrays are assumed to be sorted by increasing w\n\
\n\
.. note::\n\
     The final size of the uv plane is determined by the ratio of uvSize to\n\
     uvRes\n\
");


/*
  Module Setup - Function Definitions and Documentation
*/

static PyMethodDef GridderMethods[] = {
    {"WProjection", (PyCFunction) WProjection, METH_VARARGS|METH_KEYWORDS, WProjection_doc},
    {NULL,          NULL,                      0,                          NULL           }
};

PyDoc_STRVAR(gridder_doc, \
"Extension that implements uv plane gridding with w-projection.\n\
\n\
The functions defined in this module are:\n\
 * WProjection - w-projection gridder of visibility data\n\
\n\
See the inidividual functions for more details.");


/*
  Module Setup - Initialization
*/

MOD_INIT(_gridder) {
    char filename[256];
    PyObject *m, *pModule, *pDataPath=NULL;
    
    Py_Initialize();
    
    // Module definitions and functions
    MOD_DEF(m, "_gridder", GridderMethods, gridder_doc);
    if( m == NULL ) {
        return MOD_ERROR_VAL;
    }
    import_array();
    
    // Version and revision information
    PyModule_AddObject(m, "__version__", PyString_FromString("0.1"));
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

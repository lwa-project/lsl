#include "Python.h"
#include <cmath>
#include <math.h>
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

// Gridding convolution size on a side and oversampling factor
#define GRID_KERNEL_SIZE 7    // should be an odd number
#define GRID_KERNEL_OVERSAMPLE 64


double signed_sqrt(double data) {
    if( data >= 0 ){
        return sqrt(data);
    } else {
        return -sqrt(-data);
    }
}


// Modified Bessel function of the first find
double iv0(double x) {
    double d = 0.0, ds = 1.0, sum = 1.0;
    do {
        d += 2.0;
        ds *= x*x/(d*d);
        sum += ds;
    } while(ds > sum*1e-8);
    return sum;
}


void kaiser_bessel_1d_kernel_filler(double *kernel1D) {
    int i;
    double x, v, scaleFactor;
    
    scaleFactor = iv0(8.6);
    for(i=0; i<GRID_KERNEL_SIZE*GRID_KERNEL_OVERSAMPLE/2+1; i++) {
        x = ((double) i) / GRID_KERNEL_OVERSAMPLE;
        v = sinc(x) / scaleFactor;
        v *= iv0(8.6 * sqrt(1.0-(2*x/GRID_KERNEL_SIZE*2*x/GRID_KERNEL_SIZE)));
        // Deal with NaNs
        if(v != v) {
            if(i == 0) {
                v = 1.0;
            } else {
                v = 0.0;
            }
        }
        *(kernel1D + i) = v;
    }
}


template<typename OutType>
void w_projection_kernel(long nPixSide, 
                         double uvRes, 
                         double w, 
                         OutType* kernel) {
    long i,j,l,m;
    double cL, cM;
    
    memset(kernel, 0, sizeof(OutType)*nPixSide*nPixSide);
    
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


template<typename InType, typename OutType>
void compute_kernel_correction(long nPixSide,
                               InType *kernel,
                               OutType *corr) {
    long i, j;
    OutType temp, temp2;
    OutType *corr_full;
    corr_full = (OutType*) malloc(nPixSide*GRID_KERNEL_OVERSAMPLE * sizeof(OutType));
    memset(corr_full, 0, sizeof(OutType)*nPixSide*GRID_KERNEL_OVERSAMPLE);
    
    // Copy the kernel over
    for(i=0; i<nPixSide*GRID_KERNEL_OVERSAMPLE; i++) {
        if( i < GRID_KERNEL_SIZE*GRID_KERNEL_OVERSAMPLE/2 ) {
            *(corr_full + i) = *(kernel + i);
        } else if( nPixSide*GRID_KERNEL_OVERSAMPLE - i < GRID_KERNEL_SIZE*GRID_KERNEL_OVERSAMPLE/2) {
            *(corr_full + i) = *(kernel + nPixSide*GRID_KERNEL_OVERSAMPLE - i);
        }
    }
    
    // Inverse transform
    fftwf_plan pB;
    pB = fftwf_plan_r2r_1d(nPixSide*GRID_KERNEL_OVERSAMPLE/2,
                           corr_full, corr_full,
                           FFTW_REDFT01, FFTW_ESTIMATE);
    fftwf_execute(pB);
    fftwf_destroy_plan(pB);
    
    // Select what to keep
    for(i=0; i<nPixSide; i++) {
        if( i < nPixSide/2 ) {
            temp = *(corr_full + i) / GRID_KERNEL_OVERSAMPLE;
        } else {
            temp = *(corr_full + nPixSide - i) / GRID_KERNEL_OVERSAMPLE;
        }
        
        for(j=0; j<nPixSide; j++) {
            if( j < nPixSide/2 ) {
                temp2 = *(corr_full + j) / GRID_KERNEL_OVERSAMPLE;
            } else {
                temp2 = *(corr_full + nPixSide - j) / GRID_KERNEL_OVERSAMPLE;
            }
            
            *(corr + nPixSide*i + j) = temp * temp2;
        }
    }
    
    // Cleanup
    free(corr_full);
}


template<typename InType, typename OutCompType, typename OutRealType>
void compute_gridding(long nVis,
                      long nPixSide,
                      double uvRes,
                      double wRes,
                      double const* u,
                      double const* v,
                      double const* w,
                      InType const* vis,
                      InType const* wgt,
                      OutCompType* uv,
                      OutCompType* bm,
                      OutRealType* corr) {
    // Setup
    long i, j, l, m;
    
    Py_BEGIN_ALLOW_THREADS
    
    // Figure out the w-planes to use
    int nPlanes = 0;
    long *planeStart, *planeStop;
    planeStart = (long *) malloc(MAX_W_PLANES*sizeof(long));
    planeStop = (long *) malloc(MAX_W_PLANES*sizeof(long));
    nPlanes = compute_planes(nVis, wRes, w, planeStart, planeStop);
    
    // Fill in the 1-D gridding kernel
    double *kernel1D;
    kernel1D = (double *) malloc((GRID_KERNEL_SIZE*GRID_KERNEL_OVERSAMPLE/2+1)*sizeof(double));
    kaiser_bessel_1d_kernel_filler(kernel1D);
    
    long secStart, secStop;
    double avgW, ci, cj, temp, temp2;
    long pi, pj, gi, gj;
    Complex32 *suv, *sbm, *kern;
    static float norm = (float) 1.0 / (nPixSide * nPixSide);
    
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
        #pragma omp parallel default(shared) private(suv, sbm, kern, i, j, l, m, secStart, secStop, avgW, ci, cj, pi, pj, gi, gj, temp, temp2)
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
                
                for(m=-GRID_KERNEL_SIZE/2; m<GRID_KERNEL_SIZE/2+1; m++) {
                    pi = (long) round(ci) + m;
                    gi = (long) ((ci - pi)*GRID_KERNEL_OVERSAMPLE);
                    if(gi < 0) {
                        gi = -gi;
                    }
                    
                    temp = *(kernel1D + gi);
                    
                    if( pi < 0 ) {
                        pi += nPixSide;
                    } else if( pi >= nPixSide) {
                        pi -= nPixSide;
                    }
                    
                    for(l=-GRID_KERNEL_SIZE/2; l<GRID_KERNEL_SIZE/2+1; l++) {
                        pj = (long) round(cj) + l;
                        gj = (long) ((cj - pj)*GRID_KERNEL_OVERSAMPLE);
                        if(gj < 0) {
                            gj = -gj;
                        }
                        
                        temp2 = temp * *(kernel1D + gj);
                        
                        if( pj < 0 ) {
                            pj += nPixSide;
                        } else if( pj >= nPixSide) {
                            pj -= nPixSide;
                        }
                        
                        *(suv + nPixSide*pi + pj) += (Complex32) *(vis + i) * (float) temp2;
                        *(sbm + nPixSide*pi + pj) += (Complex32) *(wgt + i) * (float) temp2;
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
                    *(bm + i) += *(sbm+i);
                }
            }
        }
        
        fftwf_free(suv);
        fftwf_free(sbm);
        fftwf_free(kern);
    }
    
    // Correct for the kernel
    compute_kernel_correction(nPixSide, kernel1D, corr);
    
    // Cleanup
    fftwf_destroy_plan(pF);
    fftwf_destroy_plan(pR);
    fftwf_free(inP);
    
    free(kernel1D);
    
    free(planeStart);
    free(planeStop);
    
    Py_END_ALLOW_THREADS
}


static PyObject *WProjection(PyObject *self, PyObject *args, PyObject *kwds) {
    PyObject *uVec, *vVec, *wVec, *visVec, *wgtVec, *output;
    PyArrayObject *uu=NULL, *vv=NULL, *ww=NULL, *vd=NULL, *wd=NULL, *uvPlane=NULL, *bmPlane=NULL, *kernCorr=NULL;
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
    vd = (PyArrayObject *) PyArray_ContiguousFromObject(visVec,
                                                        PyArray_TYPE((PyArrayObject *) visVec),
                                                        1, 1);
    wd = (PyArrayObject *) PyArray_ContiguousFromObject(wgtVec,
                                                        PyArray_TYPE((PyArrayObject *) visVec),
                                                        1, 1);
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
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input data array to 1-D");
        goto fail;
    }
    if( wd == NULL ) {
        PyErr_Format(PyExc_RuntimeError, "Cannot cast input wgt array to 1-D and the same type as the input data");
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
    /* kernel correction */
    kernCorr = (PyArrayObject *) PyArray_SimpleNew(2, dims, NPY_FLOAT32);
    if( kernCorr == NULL ) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array - kernel correction");
        goto fail;
    }
    PyArray_FILLWBYTE(kernCorr, 0);
    
    // Get pointers to the data we need
    double *u, *v, *w;
    u = (double *) PyArray_DATA(uu);
    v = (double *) PyArray_DATA(vv);
    w = (double *) PyArray_DATA(ww);
    
    // Grid
#define LAUNCH_GRIDDER(IterType) \
    compute_gridding<IterType,Complex32,float>(nVis, nPixSide, uvRes, wRes, \
                                               u, v, w, \
                                               (IterType*) PyArray_DATA(vd), \
                                               (IterType*) PyArray_DATA(wd), \
                                               (Complex32*) PyArray_DATA(uvPlane), \
                                               (Complex32*) PyArray_DATA(bmPlane), \
                                               (float*) PyArray_DATA(kernCorr))
    switch( PyArray_TYPE(vd) ) {
      case( NPY_COMPLEX64  ): LAUNCH_GRIDDER(Complex32); break;
      case( NPY_COMPLEX128 ): LAUNCH_GRIDDER(Complex64); break;
      default: PyErr_Format(PyExc_RuntimeError, "Unsupport input data type"); goto fail;
    }
    
#undef LAUNCH_GRIDDER

    Py_XDECREF(uu);
    Py_XDECREF(vv);
    Py_XDECREF(ww);
    
    output = Py_BuildValue("(OOO)", PyArray_Return(uvPlane), PyArray_Return(bmPlane), PyArray_Return(kernCorr));
    Py_XDECREF(uvPlane);
    Py_XDECREF(bmPlane);
    Py_XDECREF(kernCorr);
    
    return output;
    
fail:
    Py_XDECREF(uu);
    Py_XDECREF(vv);
    Py_XDECREF(ww);
    Py_XDECREF(vd);
    Py_XDECREF(wd);
    Py_XDECREF(uvPlane);
    Py_XDECREF(bmPlane);
    Py_XDECREF(kernCorr);
    
    return NULL;
}

PyDoc_STRVAR(WProjection_doc, \
"w-projection gridder for uv data based on the aipy.img.ImgW class, 'Wide-\n\
field Imaging Problems in Radio Astronomy' (Cornwell et al. 2005), and\n\
'WSClean: an Implementation of a Fast, Generic Wide-Field Imager for Radio\n\
Astronomy' (Offringa et al. 2014).\n\
\n\
Input arguments are:\n\
 * u: 1-D numpy.float64 array of u coordinates\n\
 * v: 1-D numpy.float64 array of v coordinates\n\
 * w: 1-D numpy.float64 array of w coordinates\n\
 * data: 1-D numpy.complex64 or numpy.complex128 array of visibility data\n\
 * wgt: 1-D numpy.complex64 or numpy.complex128 array of visibility weight data\n\
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
 * kernCorr: 2-D numpy.float32 of the image correction for the gridding\n\
             kernel\n\
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
    
    // Version information
    PyModule_AddObject(m, "__version__", PyString_FromString("0.3"));
    
    // LSL FFTW Wisdom
    pModule = PyImport_ImportModule("lsl.common.paths");
    if( pModule != NULL ) {
        pDataPath = PyObject_GetAttrString(pModule, "DATA");
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

#include "Python.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL gofast_ARRAY_API
#include "numpy/arrayobject.h"

#include "readers.h"

/*
  OVRO-LWA triggered voltage buffer dump file data section reader
*/

PyObject *ovro_method = NULL;
PyObject *ovro_size   = NULL;


PyObject *read_ovro(PyObject *self, PyObject *args) {
    PyObject *ph, *buffer, *output;
    PyArrayObject *data;
    long i;
    int ntime, nchan, nstand, npol;
    
    if(!PyArg_ParseTuple(args, "Oiiii", &ph, &ntime, &nchan, &nstand, &npol)) {
        PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
        return NULL;
    }
    
    // Create the output data array
    npy_intp dims[4];
    // 4+4-bit Data -> 6144 samples in the data section as 12 channels, 256 stands, and 2 pols.
    dims[0] = (npy_intp) ntime;
    dims[1] = (npy_intp) nchan;
    dims[2] = (npy_intp) nstand;
    dims[3] = (npy_intp) npol;
    data = (PyArrayObject*) PyArray_ZEROS(4, dims, NPY_COMPLEX64, 0);
    if(data == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array");
        Py_XDECREF(data);
        return NULL;
    }

    // Read from the file
    if( ovro_method == NULL ) {
        ovro_method = Py_BuildValue("s", "read");
        ovro_size = Py_BuildValue("i", ntime*nchan*nstand*npol*1);
    }
    buffer = PyObject_CallMethodObjArgs(ph, ovro_method, ovro_size, NULL);
    if( buffer == NULL ) {
        if( PyObject_HasAttrString(ph, "read") ) {
            PyErr_Format(PyExc_IOError, "An error occured while reading from the file");
        } else {
            PyErr_Format(PyExc_AttributeError, "Object does not have a read() method");
        }
        Py_XDECREF(data);
        return NULL;
    } else if( PyString_GET_SIZE(buffer) != (ntime*nchan*nstand*npol*1) ) {
        PyErr_Format(EOFError, "End of file encountered during filehandle read");
        Py_XDECREF(data);
        Py_XDECREF(buffer);
        return NULL;
    }
    
    Py_BEGIN_ALLOW_THREADS
    
    // Fill the data array
    const float *fp;
    float complex *a;
    a = (float complex *) PyArray_DATA(data);
    for(i=0; i<ntime*nchan*nstand*npol; i++) {
        fp = ovroLUT[ (PyString_AS_STRING(buffer))[i] ];
        *(a + i) = fp[0] + _Complex_I * fp[1];
    }
    
    Py_END_ALLOW_THREADS
    
    output = Py_BuildValue("O", PyArray_Return(data));
    
    Py_XDECREF(buffer);
    Py_XDECREF(data);
    
    return output;
}

char read_ovro_doc[] = PyDoc_STR(\
"Function to read in the data section of a OVRO-LWA triggered voltage buffer\n\
dump file and return it as a 4-D numpy array.\n\
");

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


PyObject *read_ovro_spec(PyObject *self, PyObject *args) {
    PyObject *ph, *buffer, *output;
    PyArrayObject *data;
    long i;
    int nchan, nstand, npol;
    
    if(!PyArg_ParseTuple(args, "Oiii", &ph, &nchan, &nstand, &npol)) {
        PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
        return NULL;
    }
    
    // Create the output data array
    npy_intp dims[3];
    // 4+4-bit Data -> N samples in the data section as nchan channels, nstand stands, and npol pols.
    dims[0] = (npy_intp) nchan;
    dims[1] = (npy_intp) nstand;
    dims[2] = (npy_intp) npol;
    data = (PyArrayObject*) PyArray_ZEROS(3, dims, NPY_COMPLEX64, 0);
    if(data == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array");
        Py_XDECREF(data);
        return NULL;
    }

    // Read from the file
    if( ovro_method == NULL ) {
        ovro_method = Py_BuildValue("s", "read");
        ovro_size = Py_BuildValue("i", nchan*nstand*npol*1);
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
    } else if( PyString_GET_SIZE(buffer) != (nchan*nstand*npol*1) ) {
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
    for(i=0; i<nchan*nstand*npol; i++) {
        fp = ovroLUT[ (PyString_AS_STRING(buffer))[i] ];
        *(a + i) = fp[0] + _Complex_I * fp[1];
    }
    
    Py_END_ALLOW_THREADS
    
    output = Py_BuildValue("O", PyArray_Return(data));
    
    Py_XDECREF(buffer);
    Py_XDECREF(data);
    
    return output;
}

char read_ovro_spec_doc[] = PyDoc_STR(\
"Function to read in a single set of spectra from the data section of a OVRO-LWA\n\
triggered voltage buffer dump file and return it as a 3-D numpy.complex64 array.\n\
");
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
  TBF Reader
*/

#pragma pack(push)
#pragma pack(1)
typedef struct {
    unsigned int syncWord;
    union {
        struct {
            unsigned int frame_count:24;
            unsigned char adp_id:8;
        };
        unsigned int frame_count_word;
    };
    unsigned int second_count;
    signed short int first_chan;
    signed short int unassigned;
} TBFHeader;


typedef struct {
    signed long long timetag;
    unsigned char bytes[6144];
} TBFPayload;


typedef struct {
    TBFHeader header;
    TBFPayload payload;
} TBFFrame;
#pragma pack(pop)


PyObject *tbf_method = NULL;
PyObject *tbf_size   = NULL;


PyObject *read_tbf_ci8(PyObject *self, PyObject *args) {
    PyObject *ph, *buffer, *output, *frame, *fHeader, *fPayload, *temp;
    PyArrayObject *data;
    int i;
    TBFFrame cFrame;
    
    if(!PyArg_ParseTuple(args, "OO", &ph, &frame)) {
        PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
        return NULL;
    }
    
    // Create the output data array
    npy_intp dims[4];
    // 4+4-bit Data -> 6144 samples in the data section as 12 channels, 256 stands, and 2 pols.
    dims[0] = (npy_intp) 12;
    dims[1] = (npy_intp) 256;
    dims[2] = (npy_intp) 2;
    dims[3] = (npy_intp) 2;
    data = (PyArrayObject*) PyArray_ZEROS(4, dims, NPY_INT8, 0);
    if(data == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array");
        Py_XDECREF(data);
        return NULL;
    }

    // Read from the file
    if( tbf_method == NULL ) {
        tbf_method = Py_BuildValue("s", "read");
        tbf_size = Py_BuildValue("i", sizeof(cFrame));
    }
    buffer = PyObject_CallMethodObjArgs(ph, tbf_method, tbf_size, NULL);
    if( buffer == NULL ) {
        if( PyObject_HasAttrString(ph, "read") ) {
            PyErr_Format(PyExc_IOError, "An error occured while reading from the file");
        } else {
            PyErr_Format(PyExc_AttributeError, "Object does not have a read() method");
        }
        Py_XDECREF(data);
        return NULL;
    } else if( PyString_GET_SIZE(buffer) != sizeof(cFrame) ) {
        PyErr_Format(EOFError, "End of file encountered during filehandle read");
        Py_XDECREF(data);
        Py_XDECREF(buffer);
        return NULL;
    }
    memcpy(&cFrame, PyString_AS_STRING(buffer), sizeof(cFrame));
    Py_XDECREF(buffer);
    
    Py_BEGIN_ALLOW_THREADS
    
    // Swap the bits around
    cFrame.header.frame_count_word = __bswap_32(cFrame.header.frame_count_word);
    cFrame.header.second_count = __bswap_32(cFrame.header.second_count);
    cFrame.header.first_chan = __bswap_16(cFrame.header.first_chan);
    cFrame.payload.timetag = __bswap_64(cFrame.payload.timetag);
    
    // Fill the data array
    const signed char *fp;
    signed char *a;
    a = (signed char *) PyArray_DATA(data);
    for(i=0; i<6144; i++) {
        fp = tbfLUT[ cFrame.payload.bytes[i] ];
        *(a + 2*i + 0) = fp[0];
        *(a + 2*i + 1) = fp[1];
    }
    
    Py_END_ALLOW_THREADS
    
    // Validate
    if( !validSync5C(cFrame.header.syncWord) ) {
        PyErr_Format(SyncError, "Mark 5C sync word differs from expected");
        Py_XDECREF(data);
        return NULL;
    }
    
    // Save the data to the frame object
    // 1.  Header
    fHeader = PyObject_GetAttrString(frame, "header");
    
    temp = Py_BuildValue("B", cFrame.header.adp_id);
    PyObject_SetAttrString(fHeader, "adp_id", temp);
    Py_XDECREF(temp);
    
    temp = PyLong_FromUnsignedLong(cFrame.header.frame_count);
    PyObject_SetAttrString(fHeader, "frame_count", temp);
    Py_XDECREF(temp);
    
    temp = PyLong_FromUnsignedLong(cFrame.header.second_count);
    PyObject_SetAttrString(fHeader, "second_count", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("h", cFrame.header.first_chan);
    PyObject_SetAttrString(fHeader, "first_chan", temp);
    Py_XDECREF(temp);
    
    // 2. Data
    fPayload = PyObject_GetAttrString(frame, "payload");
    
    temp = PyLong_FromLongLong(cFrame.payload.timetag);
    PyObject_SetAttrString(fPayload, "timetag", temp);
    Py_XDECREF(temp);
    
    PyObject_SetAttrString(fPayload, "_data_ci8", PyArray_Return(data));
    
    // 3. Frame
    PyObject_SetAttrString(frame, "header", fHeader);
    PyObject_SetAttrString(frame, "payload", fPayload);
    output = Py_BuildValue("O", frame);
    
    Py_XDECREF(fHeader);
    Py_XDECREF(fPayload);
    Py_XDECREF(data);
    
    return output;
}

char read_tbf_ci8_doc[] = PyDoc_STR(\
"Function to read in a single TBW frame (header+data) and store the contents\n\
as a Frame object.\n\
\n\
.. versionadded:: 1.2.0\n\
");

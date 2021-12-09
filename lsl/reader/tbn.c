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
  TBN Reader
*/

#pragma pack(push)
#pragma pack(1)
typedef struct {
    unsigned int sync_word;
    union {
        struct {
            unsigned int frame_count:24;
            unsigned char id:8;
        };
        unsigned int frame_count_word;
    };
    unsigned int tuning_word;
    union {
        unsigned short int tbn_id;
        struct {
            unsigned short int stand:14;
            unsigned char reserved:1;
            unsigned char is_tbw:1;
        };
    };
    unsigned short int gain;
} TBNHeader;


typedef struct {
    unsigned long long timetag;
    unsigned char bytes[1024];
} TBNPayload;


typedef struct {
    TBNHeader header;
    TBNPayload payload;
} TBNFrame;
#pragma pack(pop)


PyObject *tbn_method = NULL;
PyObject *tbn_size   = NULL;


PyObject *read_tbn_cf32(PyObject *self, PyObject *args) {
    PyObject *ph, *buffer, *output, *frame, *fHeader, *fPayload, *temp;
    PyArrayObject *data;
    int i;
    TBNFrame cFrame;
    
    if(!PyArg_ParseTuple(args, "OO", &ph, &frame)) {
        PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
        return NULL;
    }
    
    // Create the output data array
    npy_intp dims[1];
    dims[0] = (npy_intp) 512;
    data = (PyArrayObject*) PyArray_ZEROS(1, dims, NPY_COMPLEX64, 0);
    if(data == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array");
        Py_XDECREF(data);
        return NULL;
    }
    
    // Read from the file
    if( tbn_method == NULL ) {
        tbn_method = Py_BuildValue("s", "read");
        tbn_size = Py_BuildValue("i", sizeof(cFrame));
    }
    buffer = PyObject_CallMethodObjArgs(ph, tbn_method, tbn_size, NULL);
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
    cFrame.header.tuning_word = __bswap_32(cFrame.header.tuning_word);
    cFrame.header.tbn_id = __bswap_16(cFrame.header.tbn_id);
    cFrame.header.gain= __bswap_16(cFrame.header.gain);
    cFrame.payload.timetag = __bswap_64(cFrame.payload.timetag);
    
    // Fill the data array
    float complex *a;
    a = (float complex *) PyArray_DATA(data);
    for(i=0; i<512; i++) {
        *(a + i) = tbnLUT[ cFrame.payload.bytes[2*i+0] ] + _Complex_I * tbnLUT[ cFrame.payload.bytes[2*i+1] ];
    }
    
    Py_END_ALLOW_THREADS
    
    // Validate
    if( !validSync5C(cFrame.header.sync_word) ) {
        PyErr_Format(SyncError, "Mark 5C sync word differs from expected");
        Py_XDECREF(data);
        return NULL;
    }
    
    // Save the data to the frame object
    // 1.  Header
    fHeader = PyObject_GetAttrString(frame, "header");
    
    temp = PyLong_FromUnsignedLong(cFrame.header.frame_count);
    PyObject_SetAttrString(fHeader, "frame_count", temp);
    Py_XDECREF(temp);
    
    temp = PyLong_FromUnsignedLong(cFrame.header.tuning_word);
    PyObject_SetAttrString(fHeader, "tuning_word", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("H", cFrame.header.tbn_id);
    PyObject_SetAttrString(fHeader, "tbn_id", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("H", cFrame.header.gain);
    PyObject_SetAttrString(fHeader, "gain", temp);
    Py_XDECREF(temp);
    
    // 2. Data
    fPayload = PyObject_GetAttrString(frame, "payload");
    
    temp = PyLong_FromUnsignedLongLong(cFrame.payload.timetag);
    PyObject_SetAttrString(fPayload, "timetag", temp);
    Py_XDECREF(temp);
    
    PyObject_SetAttrString(fPayload, "_data", PyArray_Return(data));
    
    // 3. Frame
    PyObject_SetAttrString(frame, "header", fHeader);
    PyObject_SetAttrString(frame, "payload", fPayload);
    output = Py_BuildValue("O", frame);
    
    Py_XDECREF(fHeader);
    Py_XDECREF(fPayload);
    Py_XDECREF(data);
    
    return output;
}

char read_tbn_cf32_doc[] = PyDoc_STR(\
"Function to read in a single TBN frame (header+payload) and store the contents\n\
as a Frame object.\n\
");

PyObject *read_tbn_ci8(PyObject *self, PyObject *args) {
    PyObject *ph, *buffer, *output, *frame, *fHeader, *fPayload, *temp;
    PyArrayObject *data;
    int i;
    TBNFrame cFrame;
    
    if(!PyArg_ParseTuple(args, "OO", &ph, &frame)) {
        PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
        return NULL;
    }
    
    // Create the output data array
    npy_intp dims[2];
    dims[0] = (npy_intp) 512;
    dims[1] = (npy_intp) 2;
    data = (PyArrayObject*) PyArray_ZEROS(2, dims, NPY_INT8, 0);
    if(data == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array");
        Py_XDECREF(data);
        return NULL;
    }
    
    // Read from the file
    if( tbn_method == NULL ) {
        tbn_method = Py_BuildValue("s", "read");
        tbn_size = Py_BuildValue("i", sizeof(cFrame));
    }
    buffer = PyObject_CallMethodObjArgs(ph, tbn_method, tbn_size, NULL);
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
    cFrame.header.tuning_word = __bswap_32(cFrame.header.tuning_word);
    cFrame.header.tbn_id = __bswap_16(cFrame.header.tbn_id);
    cFrame.header.gain= __bswap_16(cFrame.header.gain);
    cFrame.payload.timetag = __bswap_64(cFrame.payload.timetag);
    
    // Fill the data array
    signed char *a;
    a = (signed char *) PyArray_DATA(data);
    for(i=0; i<512; i++) {
        *(a + 2*i + 0) = tbnLUT[ cFrame.payload.bytes[2*i+0] ];
        *(a + 2*i + 1) = tbnLUT[ cFrame.payload.bytes[2*i+1] ];
    }
    
    Py_END_ALLOW_THREADS
    
    // Validate
    if( !validSync5C(cFrame.header.sync_word) ) {
        PyErr_Format(SyncError, "Mark 5C sync word differs from expected");
        Py_XDECREF(data);
        return NULL;
    }
    
    // Save the data to the frame object
    // 1.  Header
    fHeader = PyObject_GetAttrString(frame, "header");
    
    temp = PyLong_FromUnsignedLong(cFrame.header.frame_count);
    PyObject_SetAttrString(fHeader, "frame_count", temp);
    Py_XDECREF(temp);
    
    temp = PyLong_FromUnsignedLong(cFrame.header.tuning_word);
    PyObject_SetAttrString(fHeader, "tuning_word", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("H", cFrame.header.tbn_id);
    PyObject_SetAttrString(fHeader, "tbn_id", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("H", cFrame.header.gain);
    PyObject_SetAttrString(fHeader, "gain", temp);
    Py_XDECREF(temp);
    
    // 2. Data
    fPayload = PyObject_GetAttrString(frame, "payload");
    
    temp = PyLong_FromUnsignedLongLong(cFrame.payload.timetag);
    PyObject_SetAttrString(fPayload, "timetag", temp);
    Py_XDECREF(temp);
    
    PyObject_SetAttrString(fPayload, "_data", PyArray_Return(data));
    
    // 3. Frame
    PyObject_SetAttrString(frame, "header", fHeader);
    PyObject_SetAttrString(frame, "payload", fPayload);
    output = Py_BuildValue("O", frame);
    
    Py_XDECREF(fHeader);
    Py_XDECREF(fPayload);
    Py_XDECREF(data);
    
    return output;
}

char read_tbn_ci8_doc[] = PyDoc_STR(\
"Function to read in a single TBN frame (header+payload) and store the contents\n\
as a Frame object.\n\
");

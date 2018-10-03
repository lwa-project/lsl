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
DRX Reader
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
        /* Also: 
            struct {
                unsigned int frame_count:24;
                unsigned char beam:3;
                unsigned char tune:3;
                unsigned char reserved:1;
                unsigned char pol:1;
            };
        */
        unsigned int frame_count_word;
    };
    unsigned int second_count;
    unsigned short int decimation;
    unsigned short int time_offset;
} DRXHeader;


typedef struct {
    unsigned long long timetag;
    unsigned int tuning_word;
    unsigned int flags;
    unsigned char bytes[4096];
} DRXPayload;


typedef struct {
    DRXHeader header;
    DRXPayload data;
} DRXFrame;
#pragma pack(pop)


PyObject *drx_method = NULL;
PyObject *drx_size   = NULL;


PyObject *readDRX(PyObject *self, PyObject *args) {
    PyObject *ph, *buffer, *output, *frame, *fHeader, *fData, *temp;
    PyArrayObject *data;
    int i;
    DRXFrame cFrame;
    
    if(!PyArg_ParseTuple(args, "OO", &ph, &frame)) {
        PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
        return NULL;
    }
    
    // Create the output data array
    npy_intp dims[1];
    dims[0] = (npy_intp) 4096;
    data = (PyArrayObject*) PyArray_ZEROS(1, dims, NPY_COMPLEX64, 0);
    if(data == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array");
        Py_XDECREF(data);
        return NULL;
    }
    
    // Read from the file
    if( drx_method == NULL ) {
        drx_method = Py_BuildValue("s", "read");
        drx_size = Py_BuildValue("i", sizeof(cFrame));
    }
    buffer = PyObject_CallMethodObjArgs(ph, drx_method, drx_size, NULL);
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
    cFrame.header.decimation = __bswap_16(cFrame.header.decimation);
    cFrame.header.time_offset = __bswap_16(cFrame.header.time_offset);
    cFrame.data.timetag = __bswap_64(cFrame.data.timetag);
    cFrame.data.tuning_word = __bswap_32(cFrame.data.tuning_word);
    cFrame.data.flags = __bswap_32(cFrame.data.flags);
    
    // Fill the data array
    const float *fp;
    float complex *a;
    a = (float complex *) PyArray_DATA(data);
    for(i=0; i<4096; i++) {
        fp = drxLUT[ cFrame.data.bytes[i] ];
        *(a + i) = fp[0] + _Complex_I * fp[1];
    }
    
    Py_END_ALLOW_THREADS
    
    // Validate
    if( !validSync5C(cFrame.header.sync_word) ) {
        PyErr_Format(SyncError, "Mark 5C sync word differs from expected");
        Py_XDECREF(data);
        return NULL;
    }
    
    // Save the data to the frame object
    // 1. Header
    fHeader = PyObject_GetAttrString(frame, "header");
    
    temp = PyLong_FromUnsignedLong(cFrame.header.frame_count);
    PyObject_SetAttrString(fHeader, "frame_count", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("B", cFrame.header.id);
    PyObject_SetAttrString(fHeader, "drx_id", temp);
    Py_XDECREF(temp);
    
    temp = PyLong_FromUnsignedLong(cFrame.header.second_count);
    PyObject_SetAttrString(fHeader, "second_count", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("H", cFrame.header.decimation);
    PyObject_SetAttrString(fHeader, "decimation", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("H", cFrame.header.time_offset);
    PyObject_SetAttrString(fHeader, "time_offset", temp);
    Py_XDECREF(temp);
    
    // 2. Data
    fData = PyObject_GetAttrString(frame, "data");
    
    temp = PyLong_FromUnsignedLongLong(cFrame.data.timetag);
    PyObject_SetAttrString(fData, "timetag", temp);
    Py_XDECREF(temp);
    
    temp = PyLong_FromUnsignedLong(cFrame.data.tuning_word);
    PyObject_SetAttrString(fData, "tuning_word", temp);
    Py_XDECREF(temp);
    
    temp = PyLong_FromUnsignedLong(cFrame.data.flags);
    PyObject_SetAttrString(fData, "flags", temp);
    Py_XDECREF(temp);
    
    PyObject_SetAttrString(fData, "iq", PyArray_Return(data));
    
    // 3. Frame
    PyObject_SetAttrString(frame, "header", fHeader);
    PyObject_SetAttrString(frame, "data", fData);
    output = Py_BuildValue("O", frame);
    
    Py_XDECREF(fHeader);
    Py_XDECREF(fData);
    Py_XDECREF(data);
    
    return output;
}

char readDRX_doc[] = PyDoc_STR(\
"Function to read in a single DRX frame (header+data) and store the contents\n\
as a Frame object.\n\
");

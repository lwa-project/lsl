#include "Python.h"
#include <cmath>
#include <cstdint>

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL gofast_ARRAY_API
#include "numpy/arrayobject.h"

#include "readers.hpp"

/*
  TBX Reader
*/

typedef struct __attribute__((packed)) {
    uint32_t syncWord;
    union {
        struct {
            uint32_t frame_count:24;
            uint8_t  frame_id:8;
        };
        uint32_t frame_count_word;
    };
    uint32_t second_count;
    uint32_t first_chan;
    uint32_t nstand;
    uint32_t nchan;
} TBXHeader;


typedef struct __attribute__((packed)) {
    int64_t timetag;
    uint8_t bytes[8192];       // The maximum possible size
} TBXPayload;


typedef struct __attribute__((packed)) {
    TBXHeader header;
    TBXPayload payload;
} TBXFrame;

/*
  State variables to keep track of the last read_tbx_impl used
*/

static int cached_nstand = 256;
static int cached_nchan = 16;

PyObject *tbx_method = NULL;
PyObject *tbx_size = NULL;

template<int NSTAND, int NCHAN, typename T, NPY_TYPES N>
PyObject *read_tbx_impl(PyObject *self, PyObject *args) {
    PyObject *ph, *buffer, *output, *frame, *fHeader, *fPayload, *temp;
    PyArrayObject *data;
    int i, nstand, nchan;
    TBXFrame cFrame;
    
    static constexpr int frameSize = sizeof(TBXHeader)+8+NCHAN*NSTAND*2*1;
    
    if(!PyArg_ParseTuple(args, "OO", &ph, &frame)) {
        PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
        return NULL;
    }
    
    // Read from the file - header + timestamp + NSTAND-sized buffer
    if( tbx_method == NULL ) {
        tbx_method = Py_BuildValue("s", "read");
    }
    if( tbx_size == NULL ) {
        tbx_size = Py_BuildValue("i", frameSize);
    }
    buffer = PyObject_CallMethodObjArgs(ph, tbx_method, tbx_size, NULL);
    if( buffer == NULL ) {
        if( PyObject_HasAttrString(ph, "read") ) {
            PyErr_Format(PyExc_IOError, "An error occured while reading from the file");
        } else {
            PyErr_Format(PyExc_AttributeError, "Object does not have a read() method");
        }
        return NULL;
    } else if( PyBytes_GET_SIZE(buffer) != frameSize ) {
        PyErr_Format(EOFError, "End of file encountered during filehandle read");
        Py_XDECREF(buffer);
        return NULL;
    }
    memcpy(&cFrame, PyBytes_AS_STRING(buffer), frameSize);
    Py_XDECREF(buffer);
    
    // Determine the number of stands and channels in the frame 
    nstand = __bswap_32(cFrame.header.nstand);
    nchan = __bswap_32(cFrame.header.nstand);
    
    // If nstand is not what we expect, update the cache and retry with correct size
    if( nstand != NSTAND || nchan != NCHAN ) {
        cached_nstand = nstand; // Update for next time
        cached_nchan = nchan;
        Py_XDECREF(tbx_size);   // Update for next time
        tbx_size = NULL;        // Force NULL since Py_XDECREF isn't guaranteed to change tbx_size
        PyObject_CallMethod(ph, "seek", "ii", -frameSize, 1);
        
        switch(nstand) {
            case 64:
                switch(nchan) {
                    case 4:
                        return read_tbx_impl<64, 4, T, N>(self, args);
                    case 8:
                        return read_tbx_impl<64, 8, T, N>(self, args);
                    case 12:
                        return read_tbx_impl<64, 12, T, N>(self, args);
                    case 16:
                        return read_tbx_impl<64, 16, T, N>(self, args);
                    default:
                        PyErr_Format(SyncError, "Mark 5C sync word differs from expected");
                        return NULL;
                }
            case 128:
                switch(nchan) {
                    case 4:
                        return read_tbx_impl<128, 4, T, N>(self, args);
                    case 8:
                        return read_tbx_impl<128, 8, T, N>(self, args);
                    case 12:
                        return read_tbx_impl<128, 12, T, N>(self, args);
                    case 16:
                        return read_tbx_impl<128, 16, T, N>(self, args);
                    default:
                        PyErr_Format(SyncError, "Mark 5C sync word differs from expected");
                        return NULL;
                }
            case 256:
                switch(nchan) {
                    case 4:
                        return read_tbx_impl<256, 4, T, N>(self, args);
                    case 8:
                        return read_tbx_impl<256, 8, T, N>(self, args);
                    case 12:
                        return read_tbx_impl<256, 12, T, N>(self, args);
                    case 16:
                        return read_tbx_impl<256, 16, T, N>(self, args);
                    default:
                        PyErr_Format(SyncError, "Mark 5C sync word differs from expected");
                        return NULL;
                }
            default:
                PyErr_Format(SyncError, "Mark 5C sync word differs from expected");
                return NULL;
        }
    }
    
    // Create the output data array
    npy_intp dims[3];
    // 4+4-bit Data -> 8192 samples in the data section as 16 channels, 256 stands, and 2 pols.
    dims[0] = (npy_intp) NCHAN;
    dims[1] = (npy_intp) NSTAND;
    dims[2] = (npy_intp) 2;
    if(N == NPY_INT8) dims[2] *= 2;
    data = (PyArrayObject*) PyArray_ZEROS(3, dims, N, 0);
    if(data == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array");
        Py_XDECREF(data);
        return NULL;
    }
    
    Py_BEGIN_ALLOW_THREADS
    
    // Swap the bits around
    cFrame.header.frame_count_word = __bswap_32(cFrame.header.frame_count_word);
    cFrame.header.second_count = __bswap_32(cFrame.header.second_count);
    cFrame.header.first_chan = __bswap_16(cFrame.header.first_chan);
    cFrame.payload.timetag = __bswap_64(cFrame.payload.timetag);
    
    // Fill the data array
    const int8_t *fp;
    T *a;
    a = (T *) PyArray_DATA(data);
    for(i=0; i<NCHAN*NSTAND*2*1; i+=4) {
        fp = tbxLUT[ cFrame.payload.bytes[i] ];
        *(a + 2*i + 0) = fp[0];
        *(a + 2*i + 1) = fp[1];
        fp = tbxLUT[ cFrame.payload.bytes[i+1] ];
        *(a + 2*i + 2) = fp[0];
        *(a + 2*i + 3) = fp[1];
        fp = tbxLUT[ cFrame.payload.bytes[i+2] ];
        *(a + 2*i + 4) = fp[0];
        *(a + 2*i + 5) = fp[1];
        fp = tbxLUT[ cFrame.payload.bytes[i+3] ];
        *(a + 2*i + 6) = fp[0];
        *(a + 2*i + 7) = fp[1];
    }
    
    Py_END_ALLOW_THREADS
    
    // Validate
    if( !validSync5C(cFrame.header.syncWord) ) {
        buffer = PyObject_CallMethod(ph, "seek", "ii", -frameSize, 1);
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
    
    temp = Py_BuildValue("h", nstand);
    PyObject_SetAttrString(fHeader, "nstand", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("h", nchan);
    PyObject_SetAttrString(fHeader, "nchan", temp);
    Py_XDECREF(temp);
    
    // 2. Data
    fPayload = PyObject_GetAttrString(frame, "payload");
    
    temp = PyLong_FromLongLong(cFrame.payload.timetag);
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


template<typename T, NPY_TYPES N>
PyObject *read_tbx(PyObject *self, PyObject *args) {
    // Try fast path first with cached nstand
    switch(cached_nstand) {
        case 64:
            switch(cached_nchan) {
                case 4:
                    return read_tbx_impl<64, 4, T, N>(self, args);
                case 8:
                    return read_tbx_impl<64, 8, T, N>(self, args);
                case 12:
                    return read_tbx_impl<64, 12, T, N>(self, args);
                case 16:
                default:
                    return read_tbx_impl<64, 16, T, N>(self, args);
            }
        case 128:
            switch(cached_nchan) {
                case 4:
                    return read_tbx_impl<128, 4, T, N>(self, args);
                case 8:
                    return read_tbx_impl<128, 8, T, N>(self, args);
                case 12:
                    return read_tbx_impl<128, 12, T, N>(self, args);
                case 16:
                default:
                    return read_tbx_impl<128, 16, T, N>(self, args);
            }
        case 256: 
        default:
            switch(cached_nchan) {
                case 4:
                    return read_tbx_impl<256, 4, T, N>(self, args);
                case 8:
                    return read_tbx_impl<256, 8, T, N>(self, args);
                case 12:
                    return read_tbx_impl<256, 12, T, N>(self, args);
                case 16:
                default:
                    return read_tbx_impl<256, 16, T, N>(self, args);
            }
    }
}


PyObject *read_tbx_cf32(PyObject *self, PyObject *args) {
    return read_tbx<float,NPY_COMPLEX64>(self, args);
}

char read_tbx_cf32_doc[] = PyDoc_STR(\
"Function to read in a single TBW frame (header+data) and store the contents\n\
as a Frame object.\n\
\n\
.. versionadded:: 1.2.0\n\
");


PyObject *read_tbx_ci8(PyObject *self, PyObject *args) {
    return read_tbx<int8_t,NPY_INT8>(self, args);
}

char read_tbx_ci8_doc[] = PyDoc_STR(\
"Function to read in a single TBW frame (header+data) and store the contents\n\
as a Frame object.\n\
\n\
.. note::\n\
    This function differs from `read_tbx` in that it returns a\n\
    `lsl.tbx.FramePayload` containing a 4-D numpy.int8 array (channels by\n\
    stands by polarizations by real/complex) rather than a 1-D numpy.complex64\n\
    array.\n\
\n\
.. versionadded:: 2.1.3\n\
");

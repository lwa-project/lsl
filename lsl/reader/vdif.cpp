#include "Python.h"
#include <cmath>
#include <complex>

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL gofast_ARRAY_API
#include "numpy/arrayobject.h"

#include "readers.h"

/*
  VIDF Reader
*/

typedef struct __attribute__((packed)) {
    struct {
        uint32_t seconds_from_epoch:30;
        uint8_t  is_legacy:1;
        uint8_t  is_invalid:1;
    };
    struct {
        uint32_t frame_in_second:24;
        uint16_t refEpoch:6;
        uint8_t  unassigned:2;
    };
    struct {
        uint32_t frame_length:24;
        uint32_t log2_nchan:5;
        uint8_t  version:3;
    };
    struct {
        uint16_t station_id:16;
        uint16_t thread_id:10;
        uint8_t  bits_per_sample_minus_one:5;
        uint8_t  is_complex:1;
    };
} VDIFBasicHeader;

typedef struct __attribute__((packed)) {
    uint32_t extended_data_1;
    uint32_t extended_data_2;
    uint32_t extended_data_3;
    uint32_t extended_data_4;
} VDIFExtendedHeader;


PyObject *vdif_method   = NULL;
PyObject *vdif_size_hdr = NULL;
PyObject *vdif_size_ext = NULL;
PyObject *vdif_size_dat = NULL;


/* 
  Look-up Tables
*/

#define OPTIMAL_2BIT_HIGH	3.3359

static const float HiMag = OPTIMAL_2BIT_HIGH;
static const float FourBit1sigma = 2.95;

static float vdif1LUT[256][8];
static float vdif2LUT[256][4];
static float vdif4LUT[256][2];
static float vdif8LUT[256];

void initVDIFLUTs(void) {
    /*
    These look-up tables come from the VDIFIO library's decode.c file's initluts() function.
    
    Copyright (C) 2013 Walter Brisken
    */
    
    const float lut2level[2] = {-1.0, 1.0};
    const float lut4level[4] = {-HiMag, -1.0, 1.0, HiMag};
    const float lut16level[16] = {-8/FourBit1sigma,-7/FourBit1sigma,-6/FourBit1sigma,-5/FourBit1sigma,-4/FourBit1sigma,
                    -3/FourBit1sigma,-2/FourBit1sigma,-1/FourBit1sigma,0,1/FourBit1sigma,2/FourBit1sigma,
                    3/FourBit1sigma,4/FourBit1sigma,5/FourBit1sigma,6/FourBit1sigma,7/FourBit1sigma};
    int b, i, l;
    
    for(b=0; b<256; b++) {
        /* vdif1LUT */
        for(i=0; i<8; i++) {
            l = (b>>i) & 0x01;
            vdif1LUT[b][i] =  lut2level[l];
        }
        
        /* vdif2LUT */
        for(i=0; i<4; i++) {
            l = (b >> (2*i)) & 0x03;
            vdif2LUT[b][i] = lut4level[l];
        }
        
        /* vdif4LUT */
        for(i=0; i<2; i++) {
            l = (b >> (4*i)) & 0x0F;
            vdif4LUT[b][i] = lut16level[l];
        }
        
        /* vdif8LUT */
        vdif8LUT[b] = (b*2-255)/256.0;
    }
}

/* 
  Data-type specific VDIF Helper Functions
*/

template<int8_t D>
static PyArrayObject * parse_vdif(uint8_t *rawData, uint32_t dataLength, uint32_t samplesPerWord) {
    PyArrayObject *data;
    
    uint32_t nSamples;	
    nSamples = dataLength / 4 * samplesPerWord;			// bytes -> words -> samples
    
    // Create the data holders and output data array
    npy_intp dims[1];
    dims[0] = (npy_intp) nSamples;
    data = (PyArrayObject*) PyArray_ZEROS(1, dims, NPY_FLOAT32, 0);
    if(data == NULL) {
        return data;
    }
    
    // Fill the data
    float *a;
    a = (float *) PyArray_DATA(data);
    
    uint32_t i;
    float *fp;
    for(i=0; i<dataLength; i++) {
        if( D == 8 ) {
            *(a + i) = vdif8LUT[ *(rawData+i) ];
        } else if( D == 4 ) {
            fp = vdif4LUT[ *(rawData+i) ];
            *(a + 2*i + 0) = fp[0];
            *(a + 2*i + 1) = fp[1];
        } else if( D == 2 ) {
            fp = vdif2LUT[ *(rawData+i) ];
            *(a + 4*i + 0) = fp[0];
            *(a + 4*i + 1) = fp[1];
            *(a + 4*i + 2) = fp[2];
            *(a + 4*i + 3) = fp[3];
        } else if( D == 1) {
            fp = vdif1LUT[ *(rawData+i) ];
            *(a + 8*i + 0) = fp[0];
            *(a + 8*i + 1) = fp[1];
            *(a + 8*i + 2) = fp[2];
            *(a + 8*i + 3) = fp[3];
            *(a + 8*i + 4) = fp[4];
            *(a + 8*i + 5) = fp[5];
            *(a + 8*i + 6) = fp[6];
            *(a + 8*i + 7) = fp[7];
        }
    }
    
    // Done
    return data;
}


/*
  VIDF Reader
*/

PyObject *read_vdif(PyObject *self, PyObject *args, PyObject *kwds) {
    PyObject *ph, *buffer, *output, *frame, *fHeader, *fPayload, *temp;
    PyArrayObject *data=NULL;
    uint32_t i;
    float cFreq, sRate;
    cFreq = 0.0;
    sRate = 0.0;
    VDIFBasicHeader bHeader;
    VDIFExtendedHeader eHeader;
    
    char const* kwlist[] = {"fh", "frame", "central_freq", "sample_rate", NULL};
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "OO|ff", const_cast<char **>(kwlist), &ph, &frame, &cFreq, &sRate)) {
        PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
        return NULL;
    }
    
    // Read in the 16 byte common (regular + legacy) header
    if( vdif_method == NULL ) {
        vdif_method = Py_BuildValue("s", "read");
        vdif_size_hdr = Py_BuildValue("i", sizeof(bHeader));
    }
    buffer = PyObject_CallMethodObjArgs(ph, vdif_method, vdif_size_hdr, NULL);
    if( buffer == NULL ) {
        if( PyObject_HasAttrString(ph, "read") ) {
            PyErr_Format(PyExc_IOError, "An error occured while reading from the file");
        } else {
            PyErr_Format(PyExc_AttributeError, "Object does not have a read() method");
        }
        goto fail;
    } else if( PyString_GET_SIZE(buffer) != sizeof(bHeader) ) {
        PyErr_Format(EOFError, "End of file encountered during filehandle read");
        goto fail;
    }
    memcpy(&bHeader, PyString_AS_STRING(buffer), sizeof(bHeader));
    Py_XDECREF(buffer);
    
    // Fix up various bits in the basic header
    uint32_t nChan;
    uint8_t bitsPerSample;
    nChan = 1 << bHeader.log2_nchan;
    bitsPerSample = bHeader.bits_per_sample_minus_one + 1;
    
    // Does this frame look like it is valid?
    if( bHeader.frame_length < sizeof(bHeader)/8 ) {
        buffer = PyObject_CallMethod(ph, "seek", "ii", -sizeof(bHeader), 1);
        PyErr_Format(SyncError, "Frame size is zero, zero-filled frame?");
        goto fail;
    }
    
    
    if( bHeader.is_legacy == 0 ) {
        // Deal with the extra information in standard (non-legacy) headers
        if( vdif_size_ext == NULL ) {
            vdif_size_ext = Py_BuildValue("i", sizeof(eHeader));
        }
        buffer = PyObject_CallMethodObjArgs(ph, vdif_method, vdif_size_ext, NULL);
        if( buffer == NULL ) {
            if( PyObject_HasAttrString(ph, "read") ) {
                PyErr_Format(PyExc_IOError, "An error occured while reading from the file");
            } else {
                PyErr_Format(PyExc_AttributeError, "Object does not have a read() method");
            }
            goto fail;
        } else if( PyString_GET_SIZE(buffer) != sizeof(eHeader) ) {
            PyErr_Format(EOFError, "End of file encountered during filehandle read");
            goto fail;
        }
        memcpy(&eHeader, PyString_AS_STRING(buffer), sizeof(eHeader));
        Py_XDECREF(buffer);
        
    } else {
        // Legacy headers are missing this information
        eHeader.extended_data_1 = 0;
        eHeader.extended_data_2 = 0;
        eHeader.extended_data_3 = 0;
        eHeader.extended_data_4 = 0;
        
    }
    
    // Figure out how much to read in to get the entire data frame
    uint32_t dataLength, samplesPerWord, nSamples;
    dataLength = bHeader.frame_length*8 - 32 + 16*bHeader.is_legacy;	// 8-byte chunks -> bytes - full header + legacy offset
    samplesPerWord = 32 / bitsPerSample;					// dimensionless
    nSamples = dataLength / 4 * samplesPerWord;					// bytes -> words -> samples
    
    // Read in a chunk of the data
    uint8_t *rawData;
    rawData = (uint8_t *) malloc(dataLength);
    if( vdif_size_dat == NULL ) {
        vdif_size_dat = Py_BuildValue("i", dataLength);
    } else if( PyInt_AsLong(vdif_size_dat) != dataLength ) {
        Py_XDECREF(vdif_size_dat);
        vdif_size_dat = Py_BuildValue("i", dataLength);
    }
    buffer = PyObject_CallMethodObjArgs(ph, vdif_method, vdif_size_dat, NULL);
    if( buffer == NULL ) {
        if( PyObject_HasAttrString(ph, "read") ) {
            PyErr_Format(PyExc_IOError, "An error occured while reading from the file");
        } else {
            PyErr_Format(PyExc_AttributeError, "Object does not have a read() method");
        }
        free(rawData);
        goto fail;
    } else if( PyString_GET_SIZE(buffer) != dataLength ) {
        PyErr_Format(EOFError, "End of file encountered during filehandle read");
        free(rawData);
        goto fail;
    }
    memcpy(rawData, PyString_AS_STRING(buffer), sizeof(uint8_t)*dataLength);
    Py_XDECREF(buffer);
    
    // Parse it out
    if( bitsPerSample == 8 ) {
        data = parse_vdif<8>(rawData, dataLength, samplesPerWord);
    } else if( bitsPerSample == 4 ) {
        data = parse_vdif<4>(rawData, dataLength, samplesPerWord);
    } else if( bitsPerSample == 2 ) {
        data = parse_vdif<2>(rawData, dataLength, samplesPerWord);
    } else if( bitsPerSample == 1 ) {
        data = parse_vdif<1>(rawData, dataLength, samplesPerWord);
    } else {
        PyErr_Format(PyExc_RuntimeError, "Cannot parse data with %d bits per sample", bitsPerSample);
        free(rawData);
        goto fail;
    }
    
    // Clean and deal with the unexpected
    free(rawData);
    if(data == NULL) {
        PyErr_Format(PyExc_MemoryError, "Cannot create output array");
        goto fail;
    }
    
    // Deal with complex data as necessary
    if( bHeader.is_complex == 1 ) {
        PyArrayObject *tempArrayComplex;
        
        // Create a new complex64 array to hold the data
        npy_intp dims[1];
        dims[0] = (npy_intp) nSamples/2;
        tempArrayComplex = (PyArrayObject*) PyArray_ZEROS(1, dims, NPY_COMPLEX64, 0);
        if(tempArrayComplex == NULL) {
            PyErr_Format(PyExc_MemoryError, "Cannot create output array");
            Py_XDECREF(tempArrayComplex);
            goto fail;
        }
        
        float *a;
        std::complex<float> *b;
        
        a = (float *) PyArray_DATA(data);
        b = (std::complex<float> *) PyArray_DATA(tempArrayComplex);
        for(i=0; i<nSamples; i+=2) {
            *(b+i/2) = std::complex<float>(*(a+i+0), *(a+i+1));
        }
        
        Py_XDECREF(data);
        
        data = tempArrayComplex;
        nSamples /= 2;
    }
    
    // Reshape to deal with multi-channel data... someday
    if( nChan > 1 ) {
        PyArrayObject *tempArrayMulti;
        
        // Reshape
        npy_intp dims[2];
        dims[0] = (npy_intp) nSamples/nChan;
        dims[1] = (npy_intp) nChan;
        PyArray_Dims padims;
        padims.ptr = &dims[0];
        padims.len = 2;
        
        tempArrayMulti = (PyArrayObject*) PyArray_Newshape(data, &padims, NPY_CORDER);
        Py_XDECREF(data);
        
        // Transpose
        data = (PyArrayObject*) PyArray_Transpose(tempArrayMulti, NULL);
        Py_XDECREF(tempArrayMulti);
        
    }
    
    // Save the data to the frame object
    // 1.  Header
    fHeader = PyObject_GetAttrString(frame, "header");

    temp = Py_BuildValue("B", bHeader.is_invalid);
    PyObject_SetAttrString(fHeader, "is_invalid", temp);
    Py_XDECREF(temp);

    temp = Py_BuildValue("B", bHeader.is_legacy);
    PyObject_SetAttrString(fHeader, "is_legacy", temp);
    Py_XDECREF(temp);

    temp = Py_BuildValue("I", bHeader.seconds_from_epoch);
    PyObject_SetAttrString(fHeader, "seconds_from_epoch", temp);
    Py_XDECREF(temp);

    temp = Py_BuildValue("H", bHeader.refEpoch);
    PyObject_SetAttrString(fHeader, "ref_epoch", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("I", bHeader.frame_in_second);
    PyObject_SetAttrString(fHeader, "frame_in_second", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("B", bHeader.version);
    PyObject_SetAttrString(fHeader, "version", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("I", nChan);
    PyObject_SetAttrString(fHeader, "nchan", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("I", bHeader.frame_length);
    PyObject_SetAttrString(fHeader, "frame_length", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("B", bHeader.is_complex);
    PyObject_SetAttrString(fHeader, "is_complex", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("B", bitsPerSample);
    PyObject_SetAttrString(fHeader, "bits_per_sample", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("H", bHeader.thread_id);
    PyObject_SetAttrString(fHeader, "thread_id", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("H", bHeader.station_id);
    PyObject_SetAttrString(fHeader, "station_id", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("I", eHeader.extended_data_1);
    PyObject_SetAttrString(fHeader, "extended_data_1", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("I", eHeader.extended_data_2);
    PyObject_SetAttrString(fHeader, "extended_data_2", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("I", eHeader.extended_data_3);
    PyObject_SetAttrString(fHeader, "extended_data_3", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("I", eHeader.extended_data_4);
    PyObject_SetAttrString(fHeader, "extended_data_4", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("f", cFreq);
    PyObject_SetAttrString(fHeader, "central_freq", temp);
    Py_XDECREF(temp);
    
    temp = Py_BuildValue("f", sRate);
    PyObject_SetAttrString(fHeader, "sample_rate", temp);
    Py_XDECREF(temp);
    
    // 2. Data
    fPayload = PyObject_GetAttrString(frame, "payload");
    
    PyObject_SetAttrString(fPayload, "_data", PyArray_Return(data));

    // 3. Frame
    PyObject_SetAttrString(frame, "header", fHeader);
    PyObject_SetAttrString(frame, "payload", fPayload);
    output = Py_BuildValue("O", frame);
    
    Py_XDECREF(fHeader);
    Py_XDECREF(fPayload);
    Py_XDECREF(data);
    
    return output;
    
fail:
    Py_XDECREF(buffer);
    Py_XDECREF(data);
    
    return NULL;
}

char read_vdif_doc[] = PyDoc_STR(\
"Function to read in a single VDIF frame (header+payload) and store the contents\n\
as a Frame object.\n\
");

#include "Python.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#if defined(__linux__)
/* Linux */
#include <byteswap.h>
#elif defined(__APPLE__) && defined(__MACH__)
/* OSX */
#include <libkern/OSByteOrder.h>
#define __bswap_16 OSSwapInt16
#define __bswap_32 OSSwapInt32
#define __bswap_64 OSSwapInt64
#endif

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL gofast_ARRAY_API
#include "numpy/arrayobject.h"

#include "readers.h"

/*
  VIDF Reader
*/

#pragma pack(push)
#pragma pack(1)
typedef struct {
	struct {
		unsigned int secondsFromEpoch:30;
		unsigned char isLegacy:1;
		unsigned char isInvalid:1;
	};
	struct {
		unsigned int frameInSecond:24;
		unsigned short int refEpoch:6;
		unsigned char unassigned:2;
	};
	struct {
		unsigned int frameLength:24;
		unsigned int log2nChan:5;
		unsigned char version:3;
	};
	struct {
		unsigned short int stationID:16;
		unsigned short int threadID:10;
		unsigned char bitsPerSampleMinus1:5;
		unsigned char isComplex:1;
	};
} VDIFBasicHeader;

typedef struct {
	unsigned int extendedData1;
	unsigned int extendedData2;
	unsigned int extendedData3;
	unsigned int extendedData4;
} VDIFExtendedHeader;
#pragma pack(pop)


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
	
	for(b = 0; b < 256; b++)
	{
		/* vdif1LUT */
		for(i = 0; i < 8; i++)
		{
			l = (b>>i) & 0x01;
			vdif1LUT[b][i] =  lut2level[l];
		}
		
		/* vdif2LUT */
		for(i = 0; i < 4; i++)
		{
			l = (b >> (2*i)) & 0x03;
			vdif2LUT[b][i] = lut4level[l];
		}
		
		/* vdif4LUT */
		for(i = 0; i < 2; i++)
		{
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

static PyArrayObject * parseVDIF8(unsigned char *rawData, unsigned int dataLength, unsigned int samplesPerWord) {
	PyArrayObject *data;
	
	unsigned int nSamples;	
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
	
	unsigned int i;
	for(i=0; i<dataLength; i++) {
		*(a + i) = vdif8LUT[ *(rawData+i) ];
	}
	
	// Done
	return data;
}

static PyArrayObject * parseVDIF4(unsigned char *rawData, unsigned int dataLength, unsigned int samplesPerWord) {
	PyArrayObject *data;
	
	unsigned int nSamples;	
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
	
	unsigned int i;
	float *fp;
	for(i=0; i<dataLength; i++) {
		fp = vdif4LUT[ *(rawData+i) ];
		*(a + 2*i + 0) = fp[0];
		*(a + 2*i + 1) = fp[1];
	}
	
	// Done
	return data;
}

static PyArrayObject * parseVDIF2(unsigned char *rawData, unsigned int dataLength, unsigned int samplesPerWord) {
	PyArrayObject *data;
	
	unsigned int nSamples;	
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
	
	unsigned int i;
	float *fp;
	for(i=0; i<dataLength; i++) {
		fp = vdif2LUT[ *(rawData+i) ];
		*(a + 4*i + 0) = fp[0];
		*(a + 4*i + 1) = fp[1];
		*(a + 4*i + 2) = fp[2];
		*(a + 4*i + 3) = fp[3];
	}
	
	// Done
	return data;
}


static PyArrayObject * parseVDIF1(unsigned char *rawData, unsigned int dataLength, unsigned int samplesPerWord) {
	PyArrayObject *data;
	
	unsigned int nSamples;	
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
	
	unsigned int i;
	float *fp;
	for(i=0; i<dataLength; i++) {
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
	
	// Done
	return data;
}


/*
  VIDF Reader
*/

PyObject *readVDIF(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *ph, *buffer, *output, *frame, *fHeader, *fData, *temp;
	PyArrayObject *data=NULL;
	unsigned int i;
	float cFreq, sRate;
	cFreq = 0.0;
	sRate = 0.0;
	VDIFBasicHeader bHeader;
	VDIFExtendedHeader eHeader;
	
	static char *kwlist[] = {"fh", "frame", "central_freq", "sample_rate", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OO|ff", kwlist, &ph, &frame, &cFreq, &sRate)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	// Read in the 16 byte common (regular + legacy) header
	buffer = PyObject_CallMethod(ph, "read", "i", sizeof(bHeader));
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
	unsigned int nChan;
	unsigned char bitsPerSample;
	nChan = 1 << bHeader.log2nChan;
	bitsPerSample = bHeader.bitsPerSampleMinus1 + 1;
	
	// Does this frame look like it is valid?
	if( bHeader.frameLength < sizeof(bHeader)/8 ) {
		buffer = PyObject_CallMethod(ph, "seek", "ii", -sizeof(bHeader), 1);
		PyErr_Format(SyncError, "Frame size is zero, zero-filled frame?");
		goto fail;
	}
	
	
	if( bHeader.isLegacy == 0 ) {
		// Deal with the extra information in standard (non-legacy) headers
		buffer = PyObject_CallMethod(ph, "read", "i", sizeof(eHeader));
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
		eHeader.extendedData1 = 0;
		eHeader.extendedData2 = 0;
		eHeader.extendedData3 = 0;
		eHeader.extendedData4 = 0;
		
	}
	
	// Figure out how much to read in to get the entire data frame
	unsigned int dataLength, samplesPerWord, nSamples;
	dataLength = bHeader.frameLength*8 - 32 + 16*bHeader.isLegacy;	// 8-byte chunks -> bytes - full header + legacy offset
	samplesPerWord = 32 / bitsPerSample;					// dimensionless
	nSamples = dataLength / 4 * samplesPerWord;					// bytes -> words -> samples
	
	// Read in a chunk of the data
	unsigned char *rawData;
	rawData = (unsigned char *) malloc(dataLength);
	buffer = PyObject_CallMethod(ph, "read", "i", dataLength);
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
	memcpy(rawData, PyString_AS_STRING(buffer), sizeof(dataLength));
	Py_XDECREF(buffer);
	
	// Parse it out
	if( bitsPerSample == 8 ) {
		data = parseVDIF8(rawData, dataLength, samplesPerWord);
	} else {
		if( bitsPerSample == 4 ) {
			data = parseVDIF4(rawData, dataLength, samplesPerWord);
		} else {
			if( bitsPerSample == 2 ) {
				data = parseVDIF2(rawData, dataLength, samplesPerWord);
			} else {
				if( bitsPerSample == 1 ) {
					data = parseVDIF1(rawData, dataLength, samplesPerWord);
				} else {
					PyErr_Format(PyExc_RuntimeError, "Cannot parse data with %d bits per sample", bitsPerSample);
					free(rawData);
					goto fail;
				}
			}
		}
	}
	
	// Clean and deal with the unexpected
	free(rawData);
	if(data == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		goto fail;
	}
	
	// Deal with complex data as necessary
	if( bHeader.isComplex == 1 ) {
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
		float complex *b;
		
		a = (float *) PyArray_DATA(data);
		b = (float complex *) PyArray_DATA(tempArrayComplex);
		for(i=0; i<nSamples; i+=2) {
			*(b+i/2) = *(a+i+0) + _Complex_I * *(a+i+1);
		}
		
		Py_XDECREF(data);
		
		data = tempArrayComplex;		
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

	temp = Py_BuildValue("B", bHeader.isInvalid);
	PyObject_SetAttrString(fHeader, "is_invalid", temp);
	Py_XDECREF(temp);

	temp = Py_BuildValue("B", bHeader.isLegacy);
	PyObject_SetAttrString(fHeader, "is_legacy", temp);
	Py_XDECREF(temp);

	temp = Py_BuildValue("I", bHeader.secondsFromEpoch);
	PyObject_SetAttrString(fHeader, "seconds_from_epoch", temp);
	Py_XDECREF(temp);

	temp = Py_BuildValue("H", bHeader.refEpoch);
	PyObject_SetAttrString(fHeader, "ref_epoch", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("I", bHeader.frameInSecond);
	PyObject_SetAttrString(fHeader, "frame_in_second", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("B", bHeader.version);
	PyObject_SetAttrString(fHeader, "version", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("I", nChan);
	PyObject_SetAttrString(fHeader, "nchan", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("I", bHeader.frameLength);
	PyObject_SetAttrString(fHeader, "frame_length", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("B", bHeader.isComplex);
	PyObject_SetAttrString(fHeader, "is_complex", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("B", bitsPerSample);
	PyObject_SetAttrString(fHeader, "bits_per_sample", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("H", bHeader.threadID);
	PyObject_SetAttrString(fHeader, "thread_id", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("H", bHeader.stationID);
	PyObject_SetAttrString(fHeader, "station_id", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("I", eHeader.extendedData1);
	PyObject_SetAttrString(fHeader, "extended_data_1", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("I", eHeader.extendedData2);
	PyObject_SetAttrString(fHeader, "extended_data_2", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("I", eHeader.extendedData3);
	PyObject_SetAttrString(fHeader, "extended_data_3", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("I", eHeader.extendedData4);
	PyObject_SetAttrString(fHeader, "extendedData4", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("f", cFreq);
	PyObject_SetAttrString(fHeader, "central_freq", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("f", sRate);
	PyObject_SetAttrString(fHeader, "sample_rate", temp);
	Py_XDECREF(temp);
	
	// 2. Data
	fData = PyObject_GetAttrString(frame, "data");
	
	PyObject_SetAttrString(fData, "data", PyArray_Return(data));

	// 3. Frame
	PyObject_SetAttrString(frame, "header", fHeader);
	PyObject_SetAttrString(frame, "data", fData);
	output = Py_BuildValue("O", frame);
	
	Py_XDECREF(fHeader);
	Py_XDECREF(fData);
	Py_XDECREF(data);
	
	return output;
	
fail:
	Py_XDECREF(buffer);
	Py_XDECREF(data);
	
	return NULL;
}

char readVDIF_doc[] = PyDoc_STR(\
"Function to read in a single VDIF frame (header+data) and store the contents\n\
as a Frame object.  This function serves as a replacement for the pure python\n\
reader lsl.reader.vdif.readFrame.\n\
\n\
In order to use this reader in place of lsl.reader.vdif.readFrame change:\n\
\n\
\t>>> import lsl.reader.vdif as vdif\n\
\t>>> fh = open('some-vdif-file.dat', 'rb')\n\
\t>>> frame = vdif.readFrame(fh)\n\
\n\
to:\n\
\n\
\t>>> import lsl.reader.vdif as vdif\n\
\t>>> from lsl.reader._vdif import readVDIF, SyncError, EOFError\n\
\t>>> fh = open('some-vdif-file.dat', 'rb')\n\
\t>> frame = readVDIF(fh, vdif.Frame())\n\
\n\
In addition, the exceptions checked for in the try...except blocks wrapping the\n\
frame reader need to be changed to 'IOError' since SyncError and EOFError are\n\
are sub-classes of IOError.\n\
");

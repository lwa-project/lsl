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
#else
/* BSD */
#include <sys/endian.h>
#define __bswap_16 __bswap16_var
#define __bswap_32 __bswap32_var
#define __bswap_64 __bswap64_var
#endif

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL gofast_ARRAY_API
#include "numpy/arrayobject.h"

#include "readers.h"

/*
  DRX 8+8-bit Reader
*/

#pragma pack(push)
#pragma pack(1)
typedef struct {
	unsigned int syncWord;
	union {
		struct {
			unsigned int frameCount:24;
			unsigned char id:8;
		};
		/* Also:
			struct {
				unsigned int frameCount:24;
				unsigned char beam:6;
				unsigned char isADP:1;
				unsigned char pol:1;
			};
		*/
		unsigned int frameCountWord;
	};
	unsigned int secondsCount;
	signed short int decimation;
	signed short int timeOffset;
} DRX8Header;


typedef struct {
	signed long long timeTag;
	unsigned int tuningWord;
	union {
		struct {
			unsigned int flags:24;
			unsigned char drx8BW:8;
		};
		unsigned int flagsWord;
	};
	unsigned char bytes[4096];
} DRX8Payload;


typedef struct {
	DRX8Header header;
	DRX8Payload data;
} DRX8Frame;
#pragma pack(pop)


PyObject *readDRX8(PyObject *self, PyObject *args) {
	PyObject *ph, *output, *frame, *fHeader, *fData, *temp;
	PyArrayObject *data;
	int i;
	
	if(!PyArg_ParseTuple(args, "OO", &ph, &frame)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	// Read in a single 4128 byte frame
	FILE *fh = PyFile_AsFile(ph);
	PyFile_IncUseCount((PyFileObject *) ph);
	DRX8Frame cFrame;
	i = fread(&cFrame, sizeof(cFrame), 1, fh);
	if(ferror(fh)) {
		PyFile_DecUseCount((PyFileObject *) ph);
		PyErr_Format(PyExc_IOError, "An error occured while reading from the file");
		return NULL;
	}
	if(feof(fh)) {
		PyFile_DecUseCount((PyFileObject *) ph);
		PyErr_Format(eofError, "End of file encountered during filehandle read");
		return NULL;
	}
	PyFile_DecUseCount((PyFileObject *) ph);
	
	// Validate
	if( !validSync5C(cFrame.header.syncWord) ) {
		PyErr_Format(syncError, "Mark 5C sync word differs from expected");
		return NULL;
	}
	
	// Swap the bits around
	cFrame.header.frameCountWord = __bswap_32(cFrame.header.frameCountWord);
	cFrame.header.secondsCount = __bswap_32(cFrame.header.secondsCount);
	cFrame.header.decimation = __bswap_16(cFrame.header.decimation);
	cFrame.header.timeOffset = __bswap_16(cFrame.header.timeOffset);
	cFrame.data.timeTag = __bswap_64(cFrame.data.timeTag);
	cFrame.data.tuningWord = __bswap_32(cFrame.data.tuningWord);
	cFrame.data.flagsWord = __bswap_32(cFrame.data.flagsWord);
	
	// Create the output data array
	npy_intp dims[1];
	dims[0] = (npy_intp) 2048;
	data = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_COMPLEX64);
	if(data == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		return NULL;
	}
	
	// Fill the data array
	float complex *a;
	a = (float complex *) data->data;
	for(i=0; i<2048; i++) {
		*(a + i) = drx8LUT[ cFrame.data.bytes[2*i+0] ] + imaginary * drx8LUT[ cFrame.data.bytes[2*i+1] ];
	}
	
	// Save the data to the frame object
	// 1. Header
	fHeader = PyObject_GetAttrString(frame, "header");
	
	temp = PyLong_FromUnsignedLong(cFrame.header.frameCount);
	PyObject_SetAttrString(fHeader, "frameCount", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("B", cFrame.header.id);
	PyObject_SetAttrString(fHeader, "drx8ID", temp);
	Py_XDECREF(temp);
	
	temp = PyLong_FromUnsignedLong(cFrame.header.secondsCount);
	PyObject_SetAttrString(fHeader, "secondsCount", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("h", cFrame.header.decimation);
	PyObject_SetAttrString(fHeader, "decimation", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("h", cFrame.header.timeOffset);
	PyObject_SetAttrString(fHeader, "timeOffset", temp);
	Py_XDECREF(temp);
	
	// 2. Data
	fData = PyObject_GetAttrString(frame, "data");
	
	temp = PyLong_FromLongLong(cFrame.data.timeTag);
	PyObject_SetAttrString(fData, "timeTag", temp);
	Py_XDECREF(temp);
	
	temp = PyLong_FromUnsignedLong(cFrame.data.tuningWord);
	PyObject_SetAttrString(fData, "tuningWord", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("B", cFrame.data.drx8BW);
	PyObject_SetAttrString(fData, "drx8BW", temp);
	Py_XDECREF(temp);
	
	temp = PyLong_FromUnsignedLong(cFrame.data.flags);
	PyObject_SetAttrString(fData, "flags", temp);
	Py_XDECREF(temp);
	
	PyObject_SetAttrString(fData, "iq", PyArray_Return(data));
	
	// 3. Frame
	PyObject_SetAttrString(frame, "header", fHeader);
	PyObject_SetAttrString(frame, "data", fData);
	
	Py_XDECREF(fHeader);
	Py_XDECREF(fData);
	Py_XDECREF(data);
	
	output = Py_BuildValue("O", frame);
	return output;
}

char readDRX8_doc[] = PyDoc_STR(\
"Function to read in a single DRX 8+8-bit frame (header+data) and store the contents\n\
as a Frame object.\n\
\n\
.. versionadded:: 1.2.0\n\
");

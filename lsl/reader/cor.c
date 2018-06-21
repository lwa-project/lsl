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
  COR Reader
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
		unsigned int frameCountWord;
	};
	unsigned int secondsCount;
	signed short int firstChan;
	signed short int gain;
} CORHeader;


typedef struct {
	signed long long timeTag;
	signed int nAvg;
	signed short int stand0;
	signed short int stand1;
	float complex vis[COR_NCHAN*4];
} CORPayload;


typedef struct {
	CORHeader header;
	CORPayload data;
} CORFrame;
#pragma pack(pop)


PyObject *readCOR(PyObject *self, PyObject *args) {
	PyObject *ph, *buffer, *output, *frame, *fHeader, *fData, *temp;
	PyArrayObject *data=NULL;
	CORFrame cFrame;
	
	if(!PyArg_ParseTuple(args, "OO", &ph, &frame)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		goto fail;
	}
	
	// Create the output data array
	npy_intp dims[3];
	// 32+32-bit Data -> as COR_NCHAN channels, 2 pols @ stand 0, 2 pols @  stand 1
	dims[0] = COR_NCHAN;
	dims[1] = 2;
	dims[2] = 2;
	data = (PyArrayObject*) PyArray_ZEROS(3, dims, NPY_COMPLEX64, 0);
	if(data == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		goto fail;
	}
	
	// Read from the file
	buffer = PyObject_CallMethodObjArgs(ph, PyString_FromString("read"), PyInt_FromLong(sizeof(cFrame)), NULL);
	if( buffer == NULL ) {
		if( PyObject_HasAttrString(ph, "read") ) {
			PyErr_Format(PyExc_IOError, "An error occured while reading from the file");
		} else {
			PyErr_Format(PyExc_AttributeError, "Object does not have a read() method");
		}
		Py_XDECREF(data);
		return NULL;
	} else if( PyString_GET_SIZE(buffer) != sizeof(cFrame) ) {
		PyErr_Format(eofError, "End of file encountered during filehandle read");
		Py_XDECREF(data);
		Py_XDECREF(buffer);
		return NULL;
	}
	memcpy(&cFrame, PyString_AS_STRING(buffer), sizeof(cFrame));
	Py_XDECREF(buffer);
	
	Py_BEGIN_ALLOW_THREADS
	
	// Swap the bits around
	cFrame.header.frameCountWord = __bswap_32(cFrame.header.frameCountWord);
	cFrame.header.secondsCount = __bswap_32(cFrame.header.secondsCount);
	cFrame.header.firstChan = __bswap_16(cFrame.header.firstChan);
	cFrame.header.gain = __bswap_16(cFrame.header.gain);
	cFrame.data.timeTag = __bswap_64(cFrame.data.timeTag);
	cFrame.data.nAvg = __bswap_32(cFrame.data.nAvg);
	cFrame.data.stand0 = __bswap_16(cFrame.data.stand0);
	cFrame.data.stand1 = __bswap_16(cFrame.data.stand1);
	
	// Fill the data array
	float complex *a;
	a = (float complex *) PyArray_DATA(data);
	/*for(i=0; i<288; i++) {
		*(a + i) = cFrame.data.vis[i];
	}*/
	memcpy(a, &cFrame.data.vis, sizeof(float complex)*COR_NCHAN*4);
	
	Py_END_ALLOW_THREADS
	
	// Validate
	if( !validSync5C(cFrame.header.syncWord) ) {
		PyErr_Format(syncError, "Mark 5C sync word differs from expected");
		goto fail;
	}
	
	// Save the data to the frame object
	// 1.  Header
	fHeader = PyObject_GetAttrString(frame, "header");
	
	temp = Py_BuildValue("B", cFrame.header.id);
	PyObject_SetAttrString(fHeader, "id", temp);
	Py_XDECREF(temp);
	
	temp = PyLong_FromUnsignedLong(cFrame.header.frameCount);
	PyObject_SetAttrString(fHeader, "frameCount", temp);
	Py_XDECREF(temp);
	
	temp = PyLong_FromUnsignedLong(cFrame.header.secondsCount);
	PyObject_SetAttrString(fHeader, "secondsCount", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("H", cFrame.header.firstChan);
	PyObject_SetAttrString(fHeader, "firstChan", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("H", cFrame.header.gain);
	PyObject_SetAttrString(fHeader, "gain", temp);
	Py_XDECREF(temp);
	
	// 2. Data
	fData = PyObject_GetAttrString(frame, "data");
	
	temp = PyLong_FromLongLong(cFrame.data.timeTag);
	PyObject_SetAttrString(fData, "timeTag", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("i", cFrame.data.nAvg);
	PyObject_SetAttrString(fData, "nAvg", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("h", cFrame.data.stand0);
	PyObject_SetAttrString(fData, "stand0", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("h", cFrame.data.stand1);
	PyObject_SetAttrString(fData, "stand1", temp);
	Py_XDECREF(temp);
	
	PyObject_SetAttrString(fData, "vis", PyArray_Return(data));
	
	// 3. Frame
	PyObject_SetAttrString(frame, "header", fHeader);
	PyObject_SetAttrString(frame, "data", fData);
	output = Py_BuildValue("O", frame);
	
	Py_XDECREF(fHeader);
	Py_XDECREF(fData);
	Py_XDECREF(data);
	
	return output;
	
fail:
	Py_XDECREF(data);
	
	return NULL;
}

char readCOR_doc[] = PyDoc_STR(\
"Function to read in a single COR frame (header+data) and store the contents\n\
as a Frame object.\n\
\n\
.. versionchanged:: 1.2.1\n\
\tUpdated readCOR for the switch over to 72 channels, complex64 data, and no\n\
\tdata weights\n\
\n\
.. versionadded:: 1.2.0\n\
");

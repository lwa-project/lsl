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
  TBN Reader
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
	unsigned int tuningWord;
	union {
		unsigned short int tbnID;
		struct {
			unsigned short int stand:14;
			unsigned char reserved:1;
			unsigned char isTBW:1;
		};
	};
	unsigned short int gain;
} TBNHeader;


typedef struct {
	unsigned long long timeTag;
	unsigned char bytes[1024];
} TBNPayload;


typedef struct {
	TBNHeader header;
	TBNPayload data;
} TBNFrame;
#pragma pack(pop)


PyObject *readTBN(PyObject *self, PyObject *args) {
	PyObject *ph, *buffer, *output, *frame, *fHeader, *fData, *temp;
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
	buffer = PyObject_CallMethod(ph, "read", "i", sizeof(cFrame));
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
	cFrame.header.frameCountWord = __bswap_32(cFrame.header.frameCountWord);
	cFrame.header.tuningWord = __bswap_32(cFrame.header.tuningWord);
	cFrame.header.tbnID = __bswap_16(cFrame.header.tbnID);
	cFrame.header.gain= __bswap_16(cFrame.header.gain);
	cFrame.data.timeTag = __bswap_64(cFrame.data.timeTag);
	
	// Fill the data array
	float complex *a;
	a = (float complex *) PyArray_DATA(data);
	for(i=0; i<512; i++) {
		*(a + i) = tbnLUT[ cFrame.data.bytes[2*i+0] ] + _Complex_I * tbnLUT[ cFrame.data.bytes[2*i+1] ];
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
	
	temp = PyLong_FromUnsignedLong(cFrame.header.frameCount);
	PyObject_SetAttrString(fHeader, "frame_count", temp);
	Py_XDECREF(temp);
	
	temp = PyLong_FromUnsignedLong(cFrame.header.tuningWord);
	PyObject_SetAttrString(fHeader, "tuning_word", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("H", cFrame.header.tbnID);
	PyObject_SetAttrString(fHeader, "tbn_id", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("H", cFrame.header.gain);
	PyObject_SetAttrString(fHeader, "gain", temp);
	Py_XDECREF(temp);
	
	// 2. Data
	fData = PyObject_GetAttrString(frame, "data");
	
	temp = PyLong_FromUnsignedLongLong(cFrame.data.timeTag);
	PyObject_SetAttrString(fData, "timetag", temp);
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

char readTBN_doc[] = PyDoc_STR(\
"Function to read in a single TBN frame (header+data) and store the contents\n\
as a Frame object.  This function serves as a replacement for the pure python\n\
reader lsl.reader.tbn.readFrame.\n\
\n\
In order to use this reader in place of lsl.reader.tbn.readFrame change:\n\
\n\
\t>>> import lsl.reader.tbn as tbn\n\
\t>>> fh = open('some-tbn-file.dat', 'rb')\n\
\t>>> frame = tbn.readFrame(fh)\n\
\n\
to:\n\
\n\
\t>>> import lsl.reader.tbn as tbn\n\
\t>>> from lsl.reader._gofast import ReadTBN, SyncError, EOFError\n\
\t>>> fh = open('some-tbn-file.dat', 'rb')\n\
\t>> frame = readTBN(fh, tbn.Frame())\n\
\n\
In addition, the exceptions checked for in the try...except blocks wrapping the\n\
frame reader need to be changed to 'IOError' since SyncError and EOFError are\n\
are sub-classes of IOError.\n\
\n\
.. versionchanged:: 0.4.0\n\
\tThe Go Fast! readers are the default readers used by the :mod:`lsl.reader.tbn`\n\
\tmodule.\n\
");

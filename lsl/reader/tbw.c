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
  TBW Reader (12-bit and 4-bit samples)
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
	union {
		unsigned short int tbwID;
		struct {
			unsigned short int stand:14;
			unsigned char bits:1;
			unsigned char isTBW:1;
		};
	};
	unsigned short int unassigned;
} TBWHeader;


typedef struct {
	unsigned long long timeTag;
	unsigned char bytes[1200];
} TBWPayload;


typedef struct {
	TBWHeader header;
	TBWPayload data;
} TBWFrame;
#pragma pack(pop)


PyObject *readTBW(PyObject *self, PyObject *args) {
	PyObject *ph, *buffer, *output, *frame, *fHeader, *fData, *temp;
	PyArrayObject *data4, *data12;
	int i;
	TBWFrame cFrame;
	
	if(!PyArg_ParseTuple(args, "OO", &ph, &frame)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	// Create the output data arrays
	// 4-bit
	npy_intp dims[2];
	dims[0] = (npy_intp) 2;
	dims[1] = (npy_intp) 1200;
	data4 = (PyArrayObject*) PyArray_ZEROS(2, dims, NPY_INT16, 0);
	// 12-bit
	dims[0] = (npy_intp) 2;
	dims[1] = (npy_intp) 400;
	data12 = (PyArrayObject*) PyArray_ZEROS(2, dims, NPY_INT16, 0);
	if( data4 == NULL || data12 == NULL ) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		goto fail;
	}
	
	// Read from the file
	buffer = PyObject_CallMethod(ph, "read", "i", sizeof(cFrame));
	if( buffer == NULL ) {
		if( PyObject_HasAttrString(ph, "read") ) {
			PyErr_Format(PyExc_IOError, "An error occured while reading from the file");
		} else {
			PyErr_Format(PyExc_AttributeError, "Object does not have a read() method");
		}
		goto fail;
	} else if( PyString_GET_SIZE(buffer) != sizeof(cFrame) ) {
		PyErr_Format(EOFError, "End of file encountered during filehandle read");
		goto fail;
	}
	memcpy(&cFrame, PyString_AS_STRING(buffer), sizeof(cFrame));
	Py_XDECREF(buffer);
	
	Py_BEGIN_ALLOW_THREADS
	
	// Swap the bits around
	cFrame.header.frameCountWord = __bswap_32(cFrame.header.frameCountWord);
	cFrame.header.secondsCount = __bswap_32(cFrame.header.secondsCount);
	cFrame.header.tbwID = __bswap_16(cFrame.header.tbwID);
	cFrame.data.timeTag = __bswap_64(cFrame.data.timeTag);
	
	// Fill the data array
	if(cFrame.header.bits == 0) {
		short int tempR;
		short int *a;
		a = (short int *) PyArray_DATA(data12);
		for(i=0; i<400; i++) {
			tempR = (cFrame.data.bytes[3*i]<<4) | ((cFrame.data.bytes[3*i+1]>>4)&15);
			tempR -= ((tempR&2048)<<1);
			*(a + i) = (short int) tempR;

			tempR = ((cFrame.data.bytes[3*i+1]&15)<<8) | cFrame.data.bytes[3*i+2];
			tempR -= ((tempR&2048)<<1);
			*(a + 400 + i) = (short int) tempR;
		}
	} else {
		const short int *fp;
		short int *a;
		a = (short int *) PyArray_DATA(data4);
		for(i=0; i<1200; i++) {
			fp = tbw4LUT[ cFrame.data.bytes[i] ];
			*(a + i) = (short int) fp[0];
			*(a + 1200 + i) = (short int) fp[1];
		}
	}
	
	Py_END_ALLOW_THREADS
	
	// Validate
	if( !validSync5C(cFrame.header.syncWord) ) {
		PyErr_Format(SyncError, "Mark 5C sync word differs from expected");
		goto fail;
	}
	
	// Save the data to the frame object
	// 1.  Header
	fHeader = PyObject_GetAttrString(frame, "header");
	
	temp = PyLong_FromUnsignedLong(cFrame.header.frameCount);
	PyObject_SetAttrString(fHeader, "frame_count", temp);
	Py_XDECREF(temp);
	
	temp = PyLong_FromUnsignedLong(cFrame.header.secondsCount);
	PyObject_SetAttrString(fHeader, "second_count", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("H", cFrame.header.tbwID);
	PyObject_SetAttrString(fHeader, "tbw_id", temp);
	Py_XDECREF(temp);
	
	// 2. Data
	fData = PyObject_GetAttrString(frame, "data");
	
	temp = PyLong_FromUnsignedLongLong(cFrame.data.timeTag);
	PyObject_SetAttrString(fData, "timetag", temp);
	Py_XDECREF(temp);
	
	if(cFrame.header.bits == 0) {
		PyObject_SetAttrString(fData, "xy", PyArray_Return(data12));
	} else {
		PyObject_SetAttrString(fData, "xy", PyArray_Return(data4));
	}
	
	// 3. Frame
	PyObject_SetAttrString(frame, "header", fHeader);
	PyObject_SetAttrString(frame, "data", fData);
	output = Py_BuildValue("O", frame);
	
	Py_XDECREF(fHeader);
	Py_XDECREF(fData);
	Py_XDECREF(data4);
	Py_XDECREF(data12);
	
	return output;
	
fail:
	Py_XDECREF(data4);
	Py_XDECREF(data12);
	
	return NULL;
}

char readTBW_doc[] = PyDoc_STR(\
"Function to read in a single TBW frame (header+data) and store the contents\n\
as a Frame object.  This function serves as a replacement for the pure python\n\
reader lsl.reader.tbw.readFrame.\n\
\n\
In order to use this reader in place of lsl.reader.tbw.readFrame change:\n\
\n\
\t>>> import lsl.reader.tbw as tbw\n\
\t>>> fh = open('some-tbw-file.dat', 'rb')\n\
\t>>> frame = tbw.readFrame(fh)\n\
\n\
to:\n\
\n\
\t>>> import lsl.reader.tbw as tbw\n\
\t>>> from lsl.reader._gofast import ReadTBW, SyncError, EOFError\n\
\t>>> fh = open('some-tbw-file.dat', 'rb')\n\
\t>>> frame = readTBW(fh, tbw.Frame())\n\
\n\
In addition, the exceptions checked for in the try...except blocks wrapping the\n\
frame reader need to be changed to 'IOError' since SyncError and EOFError are\n\
are sub-classes of IOError.\n\
\n\
.. versionchanged:: 0.4.0\n\
\tThe Go Fast! readers are the default readers used by the :mod:`lsl.reader.tbw`\n\
\tmodule.\n\
");

#include "Python.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "numpy/arrayobject.h"

#define PI 3.1415926535898
#define imaginary _Complex_I


/*
  Minimum function for two values
*/
#define min(a,b) (((a) < (b)) ? (a) : (b))


/*
  Exceptions for the Go Fast! Readers
*/

static PyObject *syncError;
static PyObject *eofError;


/*
  Validate a collection of Mark 5C sync words.  Return 1 is all are
  valid.
*/

int validSync(unsigned char w1, unsigned char w2, unsigned char w3, unsigned char w4) {
	int valid = 0;
	
	if( w1 == 92 && w2 == 222 && w3 == 192 && w4 == 222 ) {
		valid = 1;
	}
	
	return valid;
}


/* 
  Look-up Tables
*/

static short int tbw4LUT[256][2];
static float tbnLUT[256];
static float drxLUT[256][2];

static void initLWALUTs() {
	// Look-up table inialization function from the VDIFIO library
	
	int i,j;
	short int t;
	
	//TBW - 4 bit
	for(i=0; i<256; i++) {
		for(j=0; j<2; j++) {
			t = (i >> 4*(1-j)) & 15;
			tbw4LUT[i][j] = t;
			tbw4LUT[i][j] -= ((t&8)<<1);
		}
	}
	
	// TBN
	for(i=0; i<256; i++) {
		tbnLUT[i] = i;
		tbnLUT[i] -= ((i&128)<<1);
	}
	
	// DRX
	for(i=0; i<256; i++) {
		for(j=0; j<2; j++) {
			t = (i >> 4*(1-j)) & 15;
			drxLUT[i][j] = t;
			drxLUT[i][j] -= ((t&8)<<1);
		}
	}
}

#define OPTIMAL_2BIT_HIGH	3.3359

static const float HiMag = OPTIMAL_2BIT_HIGH;
static const float FourBit1sigma = 2.95;

static float vdif1LUT[256][8];
static float vdif2LUT[256][4];
static float vdif4LUT[256][2];
static float vdif8LUT[256];

static void initVDIFLUTs() {
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
  DR Spectrometer header format
*/

#pragma pack(push)
#pragma pack(1)
typedef struct {
	unsigned int MAGIC1; 			// must always equal 0xC0DEC0DE
	unsigned long long timeTag0; 		// time tag of first frame in ``block''
	unsigned short int timeOffset;	// time offset reported by DP
	unsigned short int decFactor; 	// decimation factor
	unsigned int freqCode[2]; 		// DP frequency codes for each tuning
								//   indexing: 0..1 = Tuning 1..2
	unsigned int fills[4]; 			// fills for each pol/tuning combination
								//   indexing: 0..3 = X0, Y0 X1, Y1
	unsigned char errors[4]; 		// error flag for each pol/tuning combo
								//   indexing: 0..3 = X0, Y0 X1, Y1
	unsigned char beam;				// beam number
	unsigned char stokes_format; 		// ouptut format
	unsigned char spec_version;		// version of the spectrometer data file
	unsigned char flags;			// flag bit-field
	unsigned int nFreqs;			// <Transform Length>
	unsigned int nInts;				// <Integration Count>
	unsigned int satCount[4];		// saturation count for each pol/tuning combo
								//   indexing: 0..3 = X0, Y0 X1, Y1
	unsigned int MAGIC2;			// must always equal 0xED0CED0C
} DRSpecHeader;
#pragma pack(pop)


/*
  TBW Reader (12-bit and 4-bit samples)
*/

static PyObject *readTBW(PyObject *self, PyObject *args) {
	PyObject *ph, *output, *frame, *fHeader, *fData, *temp;
	PyArrayObject *data;
	int i;
	
	if(!PyArg_ParseTuple(args, "OO", &ph, &frame)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Read in a single 1224 byte frame
	FILE *fh = PyFile_AsFile(ph);
	PyFile_IncUseCount((PyFileObject *) ph);
	unsigned char bytes[1224];
	i = fread(bytes, 1, sizeof(bytes), fh);	
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

	// Decode the header
	unsigned char sync1, sync2, sync3, sync4;
	sync1 = bytes[3];
	sync2 = bytes[2];
	sync3 = bytes[1];
	sync4 = bytes[0];
	if( !validSync(sync1, sync2, sync3, sync4) ) {
		PyErr_Format(syncError, "Mark 5C sync word differs from expected");
		return NULL;
	}
	unsigned long int frameCount;
	frameCount = bytes[5]<<16 | bytes[6]<<8 | bytes[7];
	unsigned long int secondsCount;
	secondsCount = bytes[8]<<24 | bytes[9]<<16 | bytes[10]<<8 | bytes[11];
	unsigned short int tbwID;
	tbwID = bytes[12]<<8 | bytes[13];

	int bits = (tbwID>>14)&1;
	unsigned long long timeTag;
	timeTag = ((unsigned long long) bytes[16])<<56 | \
			((unsigned long long) bytes[17])<<48 | \
			((unsigned long long) bytes[18])<<40 | \
			((unsigned long long) bytes[19])<<32 | \
			((unsigned long long) bytes[20])<<24 | \
			((unsigned long long) bytes[21])<<16 | \
			((unsigned long long) bytes[22])<<8 | \
			bytes[23];

	// Create the output data array
	npy_intp dims[2];
	if(bits == 0) {
		// 12-bit Data -> 400 samples in the data section
		dims[0] = 2;
		dims[1] = 400;
		data = (PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_INT16);
		if(data == NULL) {
			PyErr_Format(PyExc_MemoryError, "Cannot create output array");
			Py_XDECREF(data);
			return NULL;
		}
	
		// Fill the data array
		short int tempR;
		short int *a;
		a = (short int *) data->data;
		for(i=0; i<400; i++) {
			tempR = (bytes[24+3*i]<<4) | ((bytes[24+3*i+1]>>4)&15);
			tempR -= ((tempR&2048)<<1);
			*(a + i) = (short int) tempR;

			tempR = ((bytes[24+3*i+1]&15)<<8) | bytes[24+3*i+2];
			tempR -= ((tempR&2048)<<1);
			*(a + 400 + i) = (short int) tempR;
		}
	} else {
		// 4-bit Data -> 1200 samples in the data section
		dims[0] = 2;
		dims[1] = 1200;
		data = (PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_INT16);
		if(data == NULL) {
			PyErr_Format(PyExc_MemoryError, "Cannot create output array");
			Py_XDECREF(data);
			return NULL;
		}
	
		// Fill the data array
		const short int *fp;
		short int *a;
		a = (short int *) data->data;
		for(i=0; i<1200; i++) {
			fp = tbw4LUT[ bytes[i+24] ];
			*(a + i) = (short int) fp[0];
			*(a + 1200 + i) = (short int) fp[1];
		}
		
	}
	
	// Save the data to the frame object
	// 1.  Header
	fHeader = PyObject_GetAttrString(frame, "header");
	
	temp = PyLong_FromUnsignedLong(frameCount);
	PyObject_SetAttrString(fHeader, "frameCount", temp);
	Py_XDECREF(temp);
	
	temp = PyLong_FromUnsignedLong(secondsCount);
	PyObject_SetAttrString(fHeader, "secondsCount", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("H", tbwID);
	PyObject_SetAttrString(fHeader, "tbwID", temp);
	Py_XDECREF(temp);
	
	// 2. Data
	fData = PyObject_GetAttrString(frame, "data");
	
	temp = PyLong_FromUnsignedLongLong(timeTag);
	PyObject_SetAttrString(fData, "timeTag", temp);
	Py_XDECREF(temp);
	
	PyObject_SetAttrString(fData, "xy", PyArray_Return(data));
	
	// 3. Frame
	PyObject_SetAttrString(frame, "header", fHeader);
	PyObject_SetAttrString(frame, "data", fData);
	
	Py_XDECREF(fHeader);
	Py_XDECREF(fData);
	Py_XDECREF(data);
	
	output = Py_BuildValue("O", frame);
	return output;
}

PyDoc_STRVAR(readTBW_doc, \
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
\t>>> from lsl.reader._gofast import ReadTBW, syncError, eofError\n\
\t>>> fh = open('some-tbw-file.dat', 'rb')\n\
\t>>> frame = readTBW(fh, tbw.Frame())\n\
\n\
In addition, the exceptions checked for in the try...except blocks wrapping the\n\
frame reader need to be changed to 'IOError' since syncError and eofError are\n\
are sub-classes of IOError.\n\
\n\
.. versionchanged:: 0.4.0\n\
\tThe Go Fast! readers are the default readers used by the :mod:`lsl.reader.tbw`\n\
\tmodule.\n\
");


/*
  TBN Reader
*/

static PyObject *readTBN(PyObject *self, PyObject *args) {
	PyObject *ph, *output, *frame, *fHeader, *fData, *temp;
	PyArrayObject *data;
	int i;
	
	if(!PyArg_ParseTuple(args, "OO", &ph, &frame)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	// Read in a single 1048 byte frame
	FILE *fh = PyFile_AsFile(ph);
	PyFile_IncUseCount((PyFileObject *) ph);
	unsigned char bytes[1048];
	i = fread(bytes, 1, sizeof(bytes), fh);	
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
	
	// Decode the header
	unsigned char sync1, sync2, sync3, sync4;
	sync1 = bytes[3];
	sync2 = bytes[2];
	sync3 = bytes[1];
	sync4 = bytes[0];
	if( !validSync(sync1, sync2, sync3, sync4) ) {
		PyErr_Format(syncError, "Mark 5C sync word differs from expected");
		return NULL;
	}
	unsigned long int frameCount;
	frameCount = bytes[5]<<16 | bytes[6]<<8 | bytes[7];
	unsigned long int tuningWord;
	tuningWord = bytes[8]<<24 | bytes[9]<<16 | bytes[10]<<8 | bytes[11];
	unsigned short int tbnID;
	tbnID = bytes[12]<<8 | bytes[13];
	unsigned short int gain;
	gain = bytes[14]<<8 | bytes[15];
	
	unsigned long long timeTag;
	timeTag = ((unsigned long long) bytes[16])<<56 | \
			((unsigned long long) bytes[17])<<48 | \
			((unsigned long long) bytes[18])<<40 | \
			((unsigned long long) bytes[19])<<32 | \
			((unsigned long long) bytes[20])<<24 | \
			((unsigned long long) bytes[21])<<16 | \
			((unsigned long long) bytes[22])<<8 | \
			bytes[23];
			
	// Create the output data array
	npy_intp dims[1];
	dims[0] = 512;
	data = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_COMPLEX64);
	if(data == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		return NULL;
	}
	
	// Fill the data array
	float complex *a;
	a = (float complex *) data->data;
	for(i=0; i<512; i++) {
		*(a + i) = tbnLUT[ bytes[24+2*i] ] + imaginary * tbnLUT[ bytes[24+2*i+1] ];
	}
	
	// Save the data to the frame object
	// 1.  Header
	fHeader = PyObject_GetAttrString(frame, "header");
	
	temp = PyLong_FromUnsignedLong(frameCount);
	PyObject_SetAttrString(fHeader, "frameCount", temp);
	Py_XDECREF(temp);
	
	temp = PyLong_FromUnsignedLong(tuningWord);
	PyObject_SetAttrString(fHeader, "tuningWord", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("H", tbnID);
	PyObject_SetAttrString(fHeader, "tbnID", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("H", gain);
	PyObject_SetAttrString(fHeader, "gain", temp);
	Py_XDECREF(temp);
	
	// 2. Data
	fData = PyObject_GetAttrString(frame, "data");
	
	temp = PyLong_FromUnsignedLongLong(timeTag);
	PyObject_SetAttrString(fData, "timeTag", temp);
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

PyDoc_STRVAR(readTBN_doc, \
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
\t>>> from lsl.reader._gofast import ReadTBN, syncError, eofError\n\
\t>>> fh = open('some-tbn-file.dat', 'rb')\n\
\t>> frame = readTBN(fh, tbn.Frame())\n\
\n\
In addition, the exceptions checked for in the try...except blocks wrapping the\n\
frame reader need to be changed to 'IOError' since syncError and eofError are\n\
are sub-classes of IOError.\n\
\n\
.. versionchanged:: 0.4.0\n\
\tThe Go Fast! readers are the default readers used by the :mod:`lsl.reader.tbn`\n\
\tmodule.\n\
");


/*
  DRX Reader
*/

static PyObject *readDRX(PyObject *self, PyObject *args) {
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
	unsigned char bytes[4128];
	i = fread(bytes, 1, sizeof(bytes), fh);	
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
	
	// Decode the header
	unsigned char sync1, sync2, sync3, sync4;
	sync1 = bytes[3];
	sync2 = bytes[2];
	sync3 = bytes[1];
	sync4 = bytes[0];
	if( !validSync(sync1, sync2, sync3, sync4) ) {
		PyErr_Format(syncError, "Mark 5C sync word differs from expected");
		return NULL;
	}
	unsigned char drxID;
	drxID = bytes[4];
	unsigned long int frameCount;
	frameCount = bytes[5]<<16 | bytes[6]<<8 | bytes[7];
	unsigned long int secondsCount;
	secondsCount  = bytes[8]<<24 | bytes[9]<<16 | bytes[10]<<8 | bytes[11];
	unsigned short int decimation;
	decimation = bytes[12]<<8 | bytes[13];
	unsigned short int timeOffset;
	timeOffset = bytes[14]<<8 | bytes[15];

	unsigned long long timeTag;
	timeTag = ((unsigned long long) bytes[16])<<56 | \
			((unsigned long long) bytes[17])<<48 | \
			((unsigned long long) bytes[18])<<40 | \
			((unsigned long long) bytes[19])<<32 | \
			((unsigned long long) bytes[20])<<24 | \
			((unsigned long long) bytes[21])<<16 | \
			((unsigned long long) bytes[22])<<8 | \
			bytes[23];
	unsigned long tuningWord;
	tuningWord = ((unsigned long) bytes[24])<<24 | \
			   ((unsigned long) bytes[25])<<16 | \
			   ((unsigned long) bytes[26])<<8 | \
			   bytes[27];
	unsigned long flags;
	flags  =  ((unsigned long) bytes[28])<<24 | \
			((unsigned long) bytes[29])<<16 | \
			((unsigned long) bytes[30])<<8 | \
			bytes[31];
			
	// Create the output data array
	npy_intp dims[1];
	dims[0] = 4096;
	data = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_COMPLEX64);
	if(data == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		return NULL;
	}
	
	// Fill the data array
	const float *fp;
	float complex *a;
	a = (float complex *) data->data;
	for(i=0; i<4096; i++) {
		fp = drxLUT[ bytes[i+32] ];
		*(a + i) = fp[0] + imaginary * fp[1];
	}
	
	// Save the data to the frame object
	// 1. Header
	fHeader = PyObject_GetAttrString(frame, "header");
	
	temp = PyLong_FromUnsignedLong(frameCount);
	PyObject_SetAttrString(fHeader, "frameCount", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("B", drxID);
	PyObject_SetAttrString(fHeader, "drxID", temp);
	Py_XDECREF(temp);
	
	temp = PyLong_FromUnsignedLong(secondsCount);
	PyObject_SetAttrString(fHeader, "secondsCount", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("H", decimation);
	PyObject_SetAttrString(fHeader, "decimation", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("H", timeOffset);
	PyObject_SetAttrString(fHeader, "timeOffset", temp);
	Py_XDECREF(temp);
	
	// 2. Data
	fData = PyObject_GetAttrString(frame, "data");
	
	temp = PyLong_FromUnsignedLongLong(timeTag);
	PyObject_SetAttrString(fData, "timeTag", temp);
	Py_XDECREF(temp);
	
	temp = PyLong_FromUnsignedLong(tuningWord);
	PyObject_SetAttrString(fData, "tuningWord", temp);
	Py_XDECREF(temp);
	
	temp = PyLong_FromUnsignedLong(flags);
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

PyDoc_STRVAR(readDRX_doc, \
"Function to read in a single DRX frame (header+data) and store the contents\n\
as a Frame object.  This function serves as a replacement for the pure python\n\
reader lsl.reader.drx.readFrame.\n\
\n\
In order to use this reader in place of lsl.reader.drx.readFrame change:\n\
\n\
\t>>> import lsl.reader.tbn as drx\n\
\t>>> fh = open('some-drx-file.dat', 'rb')\n\
\t>>> frame = drx.readFrame(fh)\n\
\n\
to:\n\
\n\
\t>>> import lsl.reader.drx as drx\n\
\t>>> from lsl.reader._gofast import ReadDRX, syncError, eofError\n\
\t>>> fh = open('some-drx-file.dat', 'rb')\n\
\t>>> frame = readDRX(fh, tbn.Frame())\n\
\n\
In addition, the exceptions checked for in the try...except blocks wrapping the\n\
frame reader need to be changed to 'IOError' since syncError and eofError are\n\
are sub-classes of IOError.\n\
\n\
.. versionchanged:: 0.4.0\n\
\tThe Go Fast! readers are the default readers used by the :mod:`lsl.reader.drx`\n\
\tmodule.\n\
");


/*
  DR Spectrometer Reader
*/

static PyObject *readDRSpec(PyObject *self, PyObject *args) {
	PyObject *ph, *output, *frame, *fHeader, *fData, *temp;
	PyObject *tuningWords, *fills, *errors, *saturations;
	PyArrayObject *dataA0, *dataB0, *dataC0, *dataD0;
	PyArrayObject *dataA1, *dataB1, *dataC1, *dataD1;
	
	int i, nSets;
	
	if(!PyArg_ParseTuple(args, "OO", &ph, &frame)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	// Read in a single header
	FILE *fh = PyFile_AsFile(ph);
	PyFile_IncUseCount((PyFileObject *) ph);
	DRSpecHeader header;
	i = fread(&header, sizeof(DRSpecHeader), 1, fh);
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
	
	// Check the header's magic numbers
	if( header.MAGIC1 != 0xC0DEC0DE || header.MAGIC2 != 0xED0CED0C) {
		PyFile_DecUseCount((PyFileObject *) ph);
		PyErr_Format(syncError, "Sync word differs from expected");
		return NULL;
	}
	
	// Get the data format
	i = (int) header.stokes_format;
	nSets = 0;
	if( header.stokes_format < 0x10 ) {
		// Linear
		if( header.stokes_format < 0x09 ) {
			nSets = 2;
		} else {
			if( header.stokes_format == 0x09 ) {
				nSets = 4;
			} else {
				nSets = 8;
			}
		}
	} else {
		// Stokes
		if( header.stokes_format < 0x90 ) {
			nSets = 2;
		} else {
			if( header.stokes_format == 0x90 ) {
				nSets = 4;
			} else {
				nSets = 8;
			}
		}
	}
	
	// Read in the data section
	float data[nSets*header.nFreqs];
	i = fread(&data, sizeof(float), nSets*header.nFreqs, fh);
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
	
	// Create the output data arrays
	npy_intp dims[1];
	dims[0] = (npy_intp) header.nFreqs;
	
	dataA0 = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	if(dataA0 == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array - XX0/I0");
		Py_XDECREF(dataA0);
		return NULL;
	}
	
	dataB0 = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	if(dataB0 == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array - XY0/Q0");
		Py_XDECREF(dataA0);
		Py_XDECREF(dataB0);
		return NULL;
	}
	
	dataC0 = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	if(dataC0 == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array - YX0/U0");
		Py_XDECREF(dataA0);
		Py_XDECREF(dataB0);
		Py_XDECREF(dataC0);
		return NULL;
	}
	
	dataD0 = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	if(dataD0 == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array - YY0/V0");
		Py_XDECREF(dataA0);
		Py_XDECREF(dataB0);
		Py_XDECREF(dataC0);
		Py_XDECREF(dataD0);
		return NULL;
	}
	
	dataA1 = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	if(dataA1 == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array - XX1/I1");
		Py_XDECREF(dataA0);
		Py_XDECREF(dataB0);
		Py_XDECREF(dataC0);
		Py_XDECREF(dataD0);
		Py_XDECREF(dataA1);
		return NULL;
	}
	
	dataB1 = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	if(dataB1 == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array - XY1/Q1");
		Py_XDECREF(dataA0);
		Py_XDECREF(dataB0);
		Py_XDECREF(dataC0);
		Py_XDECREF(dataD0);
		Py_XDECREF(dataA1);
		Py_XDECREF(dataB1);
		return NULL;
	}
	
	dataC1 = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	if(dataC1 == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array - YX1/U1");
		Py_XDECREF(dataA0);
		Py_XDECREF(dataB0);
		Py_XDECREF(dataC0);
		Py_XDECREF(dataD0);
		Py_XDECREF(dataA1);
		Py_XDECREF(dataB1);
		Py_XDECREF(dataC1);
		return NULL;
	}
	
	dataD1 = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	if(dataD1 == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array - YY1/V1");
		Py_XDECREF(dataA0);
		Py_XDECREF(dataB0);
		Py_XDECREF(dataC0);
		Py_XDECREF(dataD0);
		Py_XDECREF(dataA1);
		Py_XDECREF(dataB1);
		Py_XDECREF(dataC1);
		Py_XDECREF(dataD1);
		return NULL;
	}
	
	// Fill the data arrays
	float *a0, *b0, *c0, *d0, *a1, *b1, *c1, *d1;
	a0 = (float *) dataA0->data;
	b0 = (float *) dataB0->data;
	c0 = (float *) dataC0->data;
	d0 = (float *) dataD0->data;
	a1 = (float *) dataA1->data;
	b1 = (float *) dataB1->data;
	c1 = (float *) dataC1->data;
	d1 = (float *) dataD1->data;
	for(i=0; i<header.nFreqs; i++) {
		// Linear
		if( header.stokes_format == 0x01 ) {
			// XX* only
			*(a0 + i) = data[0*header.nFreqs + i] / (float) header.nFreqs / (float) header.fills[0];
			*(a1 + i) = data[1*header.nFreqs + i] / (float) header.nFreqs / (float) header.fills[2];
		} 
		if( header.stokes_format == 0x02 ) {
			// XY* only
			*(b0 + i) = data[0*header.nFreqs + i] / (float) header.nFreqs / (float) min(header.fills[0], header.fills[1]);
			*(b1 + i) = data[1*header.nFreqs + i] / (float) header.nFreqs / (float) min(header.fills[2], header.fills[3]);
		}
		if( header.stokes_format == 0x04 ) {
			// YX* only
			*(c0 + i) = data[0*header.nFreqs + i] / (float) header.nFreqs / (float) min(header.fills[0], header.fills[1]);
			*(c1 + i) = data[1*header.nFreqs + i] / (float) header.nFreqs / (float) min(header.fills[2], header.fills[3]);
		}
		if( header.stokes_format == 0x08 ) {
			// YY* only
			*(d0 + i) = data[0*header.nFreqs + i] / (float) header.nFreqs / (float) header.fills[1];
			*(d1 + i) = data[1*header.nFreqs + i] / (float) header.nFreqs / (float) header.fills[3];
		}
		if( header.stokes_format == 0x09 ) {
			// XX* and YY*
			*(a0 + i) = data[0*header.nFreqs + 2*i + 0] / (float) header.nFreqs / (float) header.fills[0];
			*(d0 + i) = data[0*header.nFreqs + 2*i + 1] / (float) header.nFreqs / (float) header.fills[1];
			*(a1 + i) = data[2*header.nFreqs + 2*i + 0] / (float) header.nFreqs / (float) header.fills[2];
			*(d1 + i) = data[2*header.nFreqs + 2*i + 1] / (float) header.nFreqs / (float) header.fills[3];
		}
		if( header.stokes_format == 0x0f ) {
			// XX*, XY*, YX*, and YY*
			*(a0 + i) = data[0*header.nFreqs + 4*i + 0] / (float) header.nFreqs / (float) header.fills[0];
			*(b0 + i) = data[0*header.nFreqs + 4*i + 1] / (float) header.nFreqs / (float) min(header.fills[0], header.fills[1]);
			*(c0 + i) = data[0*header.nFreqs + 4*i + 2] / (float) header.nFreqs / (float) min(header.fills[0], header.fills[1]);
			*(d0 + i) = data[0*header.nFreqs + 4*i + 3] / (float) header.nFreqs / (float) header.fills[1];
			*(a1 + i) = data[4*header.nFreqs + 4*i + 0] / (float) header.nFreqs / (float) header.fills[2];
			*(b1 + i) = data[4*header.nFreqs + 4*i + 1] / (float) header.nFreqs / (float) min(header.fills[2], header.fills[3]);
			*(c1 + i) = data[4*header.nFreqs + 4*i + 2] / (float) header.nFreqs / (float) min(header.fills[2], header.fills[3]);
			*(d1 + i) = data[4*header.nFreqs + 4*i + 3] / (float) header.nFreqs / (float) header.fills[3];
		}
		
		// Stokes
		if( header.stokes_format == 0x10 ) {
			// I only
			*(a0 + i) = data[0*header.nFreqs + i] / (float) header.nFreqs / (float) min(header.fills[0], header.fills[1]);
			*(a1 + i) = data[1*header.nFreqs + i] / (float) header.nFreqs / (float) min(header.fills[2], header.fills[3]);
		} 
		if( header.stokes_format == 0x20 ) {
			// Q only
			*(b0 + i) = data[0*header.nFreqs + i] / (float) header.nFreqs / (float) min(header.fills[0], header.fills[1]);
			*(b1 + i) = data[1*header.nFreqs + i] / (float) header.nFreqs / (float) min(header.fills[2], header.fills[3]);
		}
		if( header.stokes_format == 0x40 ) {
			// U only
			*(c0 + i) = data[0*header.nFreqs + i] / (float) header.nFreqs / (float) min(header.fills[0], header.fills[1]);
			*(c1 + i) = data[1*header.nFreqs + i] / (float) header.nFreqs / (float) min(header.fills[2], header.fills[3]);
		}
		if( header.stokes_format == 0x80 ) {
			// V only
			*(d0 + i) = data[0*header.nFreqs + i] / (float) header.nFreqs / (float) min(header.fills[0], header.fills[1]);
			*(d1 + i) = data[1*header.nFreqs + i] / (float) header.nFreqs / (float) min(header.fills[2], header.fills[3]);
		}
		if( header.stokes_format == 0x90 ) {
			// I and V
			*(a0 + i) = data[0*header.nFreqs + 2*i + 0] / (float) header.nFreqs / (float) min(header.fills[0], header.fills[1]);
			*(d0 + i) = data[0*header.nFreqs + 2*i + 1] / (float) header.nFreqs / (float) min(header.fills[0], header.fills[1]);
			*(a1 + i) = data[2*header.nFreqs + 2*i + 0] / (float) header.nFreqs / (float) min(header.fills[2], header.fills[3]);
			*(d1 + i) = data[2*header.nFreqs + 2*i + 1] / (float) header.nFreqs / (float) min(header.fills[2], header.fills[3]);
		}
		if( header.stokes_format == 0xf0 ) {
			// I, Q, U, and V
			*(a0 + i) = data[0*header.nFreqs + 4*i + 0] / (float) header.nFreqs / (float) min(header.fills[0], header.fills[1]);
			*(b0 + i) = data[0*header.nFreqs + 4*i + 1] / (float) header.nFreqs / (float) min(header.fills[0], header.fills[1]);
			*(c0 + i) = data[0*header.nFreqs + 4*i + 2] / (float) header.nFreqs / (float) min(header.fills[0], header.fills[1]);
			*(d0 + i) = data[0*header.nFreqs + 4*i + 3] / (float) header.nFreqs / (float) min(header.fills[0], header.fills[1]);
			*(a1 + i) = data[4*header.nFreqs + 4*i + 0] / (float) header.nFreqs / (float) min(header.fills[2], header.fills[3]);
			*(b1 + i) = data[4*header.nFreqs + 4*i + 1] / (float) header.nFreqs / (float) min(header.fills[2], header.fills[3]);
			*(c1 + i) = data[4*header.nFreqs + 4*i + 2] / (float) header.nFreqs / (float) min(header.fills[2], header.fills[3]);
			*(d1 + i) = data[4*header.nFreqs + 4*i + 3] / (float) header.nFreqs / (float) min(header.fills[2], header.fills[3]);
		}
	}
	
	// Fill in the multi-value fields (tuningWords, fills, errors)
	tuningWords = PyList_New(2);
	if(tuningWords == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output list - tuningWords");
		Py_XDECREF(dataA0);
		Py_XDECREF(dataB0);
		Py_XDECREF(dataC0);
		Py_XDECREF(dataD0);
		Py_XDECREF(dataA1);
		Py_XDECREF(dataB1);
		Py_XDECREF(dataC1);
		Py_XDECREF(dataD1);
		Py_XDECREF(tuningWords);
		return NULL;
	}
	for(i=0; i<2; i++) {
		temp = Py_BuildValue("I", header.freqCode[i]);
		PyList_SetItem(tuningWords, i, temp);
	}
	
	fills = PyList_New(4);
	if(fills == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output list - fills");
		Py_XDECREF(dataA0);
		Py_XDECREF(dataB0);
		Py_XDECREF(dataC0);
		Py_XDECREF(dataD0);
		Py_XDECREF(dataA1);
		Py_XDECREF(dataB1);
		Py_XDECREF(dataC1);
		Py_XDECREF(dataD1);
		Py_XDECREF(tuningWords);
		Py_XDECREF(fills);
		return NULL;
	}
	for(i=0; i<4; i++) {
		temp = Py_BuildValue("I", header.fills[i]);
		PyList_SetItem(fills, i, temp);
	}
	
	errors = PyList_New(4);
	if(errors == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output list - flags");
		Py_XDECREF(dataA0);
		Py_XDECREF(dataB0);
		Py_XDECREF(dataC0);
		Py_XDECREF(dataD0);
		Py_XDECREF(dataA1);
		Py_XDECREF(dataB1);
		Py_XDECREF(dataC1);
		Py_XDECREF(dataD1);
		Py_XDECREF(tuningWords);
		Py_XDECREF(fills);
		Py_XDECREF(errors);
		return NULL;
	}
	for(i=0; i<4; i++) {
		temp = Py_BuildValue("B", header.errors[i]);
		PyList_SetItem(errors, i, temp);
	}
	
	saturations = PyList_New(4);
	if(saturations == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output list - saturations");
		Py_XDECREF(dataA0);
		Py_XDECREF(dataB0);
		Py_XDECREF(dataC0);
		Py_XDECREF(dataD0);
		Py_XDECREF(dataA1);
		Py_XDECREF(dataB1);
		Py_XDECREF(dataC1);
		Py_XDECREF(dataD1);
		Py_XDECREF(tuningWords);
		Py_XDECREF(fills);
		Py_XDECREF(errors);
		Py_XDECREF(saturations);
		return NULL;
	}
	for(i=0; i<4; i++) {
		temp = Py_BuildValue("I", header.satCount[i]);
		PyList_SetItem(saturations, i, temp);
	}
	
	// Save the data to the frame object
	// 1. Header
	fHeader = PyObject_GetAttrString(frame, "header");
	
	temp = Py_BuildValue("B", header.beam);
	PyObject_SetAttrString(fHeader, "beam", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("B", header.stokes_format);
	PyObject_SetAttrString(fHeader, "format", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("H", header.decFactor);
	PyObject_SetAttrString(fHeader, "decimation", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("H", header.timeOffset);
	PyObject_SetAttrString(fHeader, "timeOffset", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("I", header.nInts);
	PyObject_SetAttrString(fHeader, "nInts", temp);
	Py_XDECREF(temp);
	
	// 2. Data
	fData = PyObject_GetAttrString(frame, "data");
	
	temp = PyLong_FromUnsignedLongLong(header.timeTag0);
	PyObject_SetAttrString(fData, "timeTag", temp);
	Py_XDECREF(temp);
	
	PyObject_SetAttrString(fData, "tuningWords", tuningWords);
	
	PyObject_SetAttrString(fData, "fills", fills);
	
	PyObject_SetAttrString(fData, "errors", errors);
	
	PyObject_SetAttrString(fData, "saturations", saturations);
	
	// Linear
	if( header.stokes_format & 0x01 ) {
		PyObject_SetAttrString(fData, "XX0", PyArray_Return(dataA0));
		PyObject_SetAttrString(fData, "XX1", PyArray_Return(dataA1));
	}
	if( header.stokes_format & 0x02 ) {
		PyObject_SetAttrString(fData, "XY0", PyArray_Return(dataB0));
		PyObject_SetAttrString(fData, "XY1", PyArray_Return(dataB1));
	}
	if( header.stokes_format & 0x04 ) {
		PyObject_SetAttrString(fData, "YX0", PyArray_Return(dataC0));
		PyObject_SetAttrString(fData, "YX1", PyArray_Return(dataC1));
	}
	if( header.stokes_format & 0x08 ) {
		PyObject_SetAttrString(fData, "YY0", PyArray_Return(dataD0));
		PyObject_SetAttrString(fData, "YY1", PyArray_Return(dataD1));
	}
	
	// Stokes
	if( header.stokes_format & 0x10 ) {
		PyObject_SetAttrString(fData, "I0", PyArray_Return(dataA0));
		PyObject_SetAttrString(fData, "I1", PyArray_Return(dataA1));
	}
	if( header.stokes_format & 0x20 ) {
		PyObject_SetAttrString(fData, "Q0", PyArray_Return(dataB0));
		PyObject_SetAttrString(fData, "Q1", PyArray_Return(dataB1));
	}
	if( header.stokes_format & 0x40 ) {
		PyObject_SetAttrString(fData, "U0", PyArray_Return(dataC0));
		PyObject_SetAttrString(fData, "U1", PyArray_Return(dataC1));
	}
	if( header.stokes_format & 0x80 ) {
		PyObject_SetAttrString(fData, "V0", PyArray_Return(dataD0));
		PyObject_SetAttrString(fData, "V1", PyArray_Return(dataD1));
	}
	
	// 3. Frame
	PyObject_SetAttrString(frame, "header", fHeader);
	PyObject_SetAttrString(frame, "data", fData);
	
	Py_XDECREF(fHeader);
	Py_XDECREF(tuningWords);
	Py_XDECREF(fills);
	Py_XDECREF(errors);
	Py_XDECREF(saturations);
	Py_XDECREF(fData);
	Py_XDECREF(dataA0);
	Py_XDECREF(dataB0);
	Py_XDECREF(dataC0);
	Py_XDECREF(dataD0);
	Py_XDECREF(dataA1);
	Py_XDECREF(dataB1);
	Py_XDECREF(dataC1);
	Py_XDECREF(dataD1);
	
	output = Py_BuildValue("O", frame);
	return output;
}

PyDoc_STRVAR(readDRSpec_doc, \
"Function to read in a DR spectrometer header structure and data and return\n\
a drspec.Frame instance.\n\
\n\
.. note::\n\
\tThis function normalizes the spectra by the number of relevant fills.  For\n\
\tproducts that are a function of more than one primary input, i.e., XY* or\n\
\tI, the minimum fill of X and Y are used for normalization.");


/* 
  Data-type specific VDIF Helper Functions
*/

PyArrayObject * parseVDIF8(unsigned char *rawData, unsigned int dataLength, unsigned int samplesPerWord, unsigned int bitsPerSample) {
	PyArrayObject *data;
	
	unsigned int nSamples;	
	nSamples = dataLength / 4 * samplesPerWord;			// bytes -> words -> samples
	
	// Create the data holders and output data array
	npy_intp dims[1];
	dims[0] = (npy_intp) nSamples;
	data = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	if(data == NULL) {
		return data;
	}
	
	// Fill the data
	float *a;
	
	a = (float *) data->data;
	
	unsigned int i;
	for(i=0; i<dataLength; i++) {
		*(a + i) = vdif8LUT[ *(rawData+i) ];
	}
	
	// Done
	return data;
}

PyArrayObject * parseVDIF4(unsigned char *rawData, unsigned int dataLength, unsigned int samplesPerWord, unsigned int bitsPerSample) {
	PyArrayObject *data;
	
	unsigned int nSamples;	
	nSamples = dataLength / 4 * samplesPerWord;			// bytes -> words -> samples
	
	// Create the data holders and output data array
	npy_intp dims[1];
	dims[0] = (npy_intp) nSamples;
	data = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	if(data == NULL) {
		return data;
	}
	
	// Fill the data
	float *a;
	
	a = (float *) data->data;
	
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

PyArrayObject * parseVDIF2(unsigned char *rawData, unsigned int dataLength, unsigned int samplesPerWord, unsigned int bitsPerSample) {
	PyArrayObject *data;
	
	unsigned int nSamples;	
	nSamples = dataLength / 4 * samplesPerWord;			// bytes -> words -> samples
	
	// Create the data holders and output data array
	npy_intp dims[1];
	dims[0] = (npy_intp) nSamples;
	data = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	if(data == NULL) {
		return data;
	}
	
	// Fill the data
	float *a;
	
	a = (float *) data->data;
	
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


PyArrayObject * parseVDIF1(unsigned char *rawData, unsigned int dataLength, unsigned int samplesPerWord, unsigned int bitsPerSample) {
	PyArrayObject *data;
	
	unsigned int nSamples;	
	nSamples = dataLength / 4 * samplesPerWord;			// bytes -> words -> samples
	
	// Create the data holders and output data array
	npy_intp dims[1];
	dims[0] = (npy_intp) nSamples;
	data = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	if(data == NULL) {
		return data;
	}
	
	// Fill the data
	float *a;
	
	a = (float *) data->data;
	
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

static PyObject *readVDIF(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *ph, *output, *frame, *fHeader, *fData, *temp;
	PyArrayObject *data;
	unsigned int i;
	float cFreq, sRate;
	cFreq = 0.0;
	sRate = 0.0;
	char message[256];
	
	static char *kwlist[] = {"fh", "frame", "centralFreq", "sampleRate", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OO|ff", kwlist, &ph, &frame, &cFreq, &sRate)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Read in the 16 byte common (regular + legacy) header
	FILE *fh = PyFile_AsFile(ph);
	PyFile_IncUseCount((PyFileObject *) ph);
	unsigned char bytes[16];
	i = fread(bytes, 1, sizeof(bytes), fh);	
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

	// Decode the common header
	unsigned char isInvalid, isLegacy;
	unsigned int secondsFromEpoch;
	isInvalid = (bytes[3] >> 30) & 1;
	isLegacy = (bytes[3] >> 29) & 1;
	secondsFromEpoch = ((bytes[3] << 24) | (bytes[2] << 16) | (bytes[1] << 8) | bytes[0]) & 0x7FFFFFFF;
	
	unsigned short int refEpoch;
	unsigned int frameInSecond;
	refEpoch = (bytes[7] & 0x3F);
	frameInSecond = (bytes[6] << 16) | (bytes[5] << 8) | bytes[4];
	
	unsigned char version;
	unsigned int nChan;
	unsigned int frameLength;
	version = (bytes[11] >> 5) & 0x7;
	nChan = 1 << (bytes[11] & 0x1F);
	frameLength = (bytes[10] << 16) | (bytes[9] << 8) | bytes[8];
	
	unsigned char isComplex, bitsPerSample;
	unsigned short int threadID, stationID;
	isComplex = (bytes[15] >> 7) & 1;
	bitsPerSample = ((bytes[15] >> 2) & 0x1F) + 1;
	threadID = ((bytes[15] & 3) <<  8) | bytes[14];
	stationID = (bytes[13] << 8) | bytes[12];
	
	// Does this frame look like it is valid?
	if( frameLength < sizeof(bytes)/8 ) {
		fseek(fh, -sizeof(bytes), SEEK_CUR);
		PyFile_DecUseCount((PyFileObject *) ph);
		PyErr_Format(syncError, "Frame size is zero, zero-filled frame?");
		return NULL;
	}
	
	unsigned int extendedData1, extendedData2, extendedData3, extendedData4;
	if( isLegacy == 0 ) {
		// Deal with the extra information in standard (non-legacy) headers
		PyFile_IncUseCount((PyFileObject *) ph);
		i = fread(bytes, 1, sizeof(bytes), fh);	
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
		
		extendedData1 = (bytes[ 3] << 24) | (bytes[ 2] << 16) | (bytes[ 1] << 8) | bytes[ 0];
		extendedData2 = (bytes[ 7] << 24) | (bytes[ 6] << 16) | (bytes[ 5] << 8) | bytes[ 4];
		extendedData3 = (bytes[11] << 24) | (bytes[10] << 16) | (bytes[ 9] << 8) | bytes[ 8];
		extendedData4 = (bytes[15] << 24) | (bytes[14] << 16) | (bytes[13] << 8) | bytes[12];
		
	} else {
		// Legacy headers are missing this information
		extendedData1 = 0;
		extendedData2 = 0;
		extendedData3 = 0;
		extendedData4 = 0;
		
	}
	
	// Figure out how much to read in to get the entire data frame
	unsigned int dataLength, samplesPerWord, nSamples;
	dataLength = frameLength*8 - 32 + 16*isLegacy;		// 8-byte chunks -> bytes - full header + legacy offset
	samplesPerWord = 32 / bitsPerSample;				// dimensionless
	nSamples = dataLength / 4 * samplesPerWord;			// bytes -> words -> samples
	
	// Read in a chunk of the data
	unsigned char *rawData;
	rawData = (unsigned char *) malloc(dataLength);
	
	PyFile_IncUseCount((PyFileObject *) ph);
	i = fread(rawData, 1, dataLength, fh);
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
	
	// Parse it out
	if( bitsPerSample == 8 ) {
		data = parseVDIF8(rawData, dataLength, samplesPerWord, bitsPerSample);
	} else {
		if( bitsPerSample == 4 ) {
			data = parseVDIF4(rawData, dataLength, samplesPerWord, bitsPerSample);
		} else {
			if( bitsPerSample == 2 ) {
				data = parseVDIF2(rawData, dataLength, samplesPerWord, bitsPerSample);
			} else {
				if( bitsPerSample == 1 ) {
					data = parseVDIF1(rawData, dataLength, samplesPerWord, bitsPerSample);
				} else {
					sprintf(message, "Cannot parse data with %i bits per sample", bitsPerSample);
					PyErr_Format(PyExc_RuntimeError, message);
					free(rawData);
					return NULL;
				}
			}
		}
	}
	
	// Clean and deal with the unexpected
	free(rawData);
	if(data == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		return NULL;
	}
	
	// Deal with complex data as necessary
	if( isComplex == 1 ) {
		PyArrayObject *tempArrayComplex;
		
		// Create a new complex64 array to hold the data
		npy_intp dims[1];
		dims[0] = (npy_intp) nSamples/2;
		tempArrayComplex = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_COMPLEX64);
		if(tempArrayComplex == NULL) {
			PyErr_Format(PyExc_MemoryError, "Cannot create output array");
			Py_XDECREF(tempArrayComplex);
			Py_XDECREF(data);
			return NULL;
		}
		
		float *a;
		float complex *b;
		
		a = (float *) data->data;
		b = (float complex *) tempArrayComplex->data;
		for(i=0; i<nSamples; i+=2) {
			*(b+i/2) = *(a+i+0) + imaginary * *(a+i+1);
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

	temp = Py_BuildValue("B", isInvalid);
	PyObject_SetAttrString(fHeader, "isInvalid", temp);
	Py_XDECREF(temp);

	temp = Py_BuildValue("B", isLegacy);
	PyObject_SetAttrString(fHeader, "isLegacy", temp);
	Py_XDECREF(temp);

	temp = Py_BuildValue("I", secondsFromEpoch);
	PyObject_SetAttrString(fHeader, "secondsFromEpoch", temp);
	Py_XDECREF(temp);

	temp = Py_BuildValue("H", refEpoch);
	PyObject_SetAttrString(fHeader, "refEpoch", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("I", frameInSecond);
	PyObject_SetAttrString(fHeader, "frameInSecond", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("B", version);
	PyObject_SetAttrString(fHeader, "version", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("I", nChan);
	PyObject_SetAttrString(fHeader, "nChan", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("I", frameLength);
	PyObject_SetAttrString(fHeader, "frameLength", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("B", isComplex);
	PyObject_SetAttrString(fHeader, "isComplex", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("B", bitsPerSample);
	PyObject_SetAttrString(fHeader, "bitsPerSample", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("H", threadID);
	PyObject_SetAttrString(fHeader, "threadID", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("H", stationID);
	PyObject_SetAttrString(fHeader, "stationID", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("I", extendedData1);
	PyObject_SetAttrString(fHeader, "extendedData1", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("I", extendedData2);
	PyObject_SetAttrString(fHeader, "extendedData2", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("I", extendedData3);
	PyObject_SetAttrString(fHeader, "extendedData3", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("I", extendedData4);
	PyObject_SetAttrString(fHeader, "extendedData4", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("f", cFreq);
	PyObject_SetAttrString(fHeader, "centralFreq", temp);
	Py_XDECREF(temp);
	
	temp = Py_BuildValue("f", sRate);
	PyObject_SetAttrString(fHeader, "sampleRate", temp);
	Py_XDECREF(temp);
	
	// 2. Data
	fData = PyObject_GetAttrString(frame, "data");
	
	PyObject_SetAttrString(fData, "data", PyArray_Return(data));

	// 3. Frame
	PyObject_SetAttrString(frame, "header", fHeader);
	PyObject_SetAttrString(frame, "data", fData);
	
	Py_XDECREF(fHeader);
	Py_XDECREF(fData);
	Py_XDECREF(data);

	output = Py_BuildValue("O", frame);
	return output;
}

PyDoc_STRVAR(readVDIF_doc, \
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
\t>>> from lsl.reader._vdif import readVDIF, syncError, eofError\n\
\t>>> fh = open('some-vdif-file.dat', 'rb')\n\
\t>> frame = readVDIF(fh, vdif.Frame())\n\
\n\
In addition, the exceptions checked for in the try...except blocks wrapping the\n\
frame reader need to be changed to 'IOError' since syncError and eofError are\n\
are sub-classes of IOError.\n\
");


/*
  Module Setup - Function Definitions and Documentation
*/

static PyMethodDef GoFastMethods[] = {
	{"readTBW",    (PyCFunction) readTBW,    METH_VARARGS,               readTBW_doc   }, 
	{"readTBN",    (PyCFunction) readTBN,    METH_VARARGS,               readTBN_doc   }, 
	{"readDRX",    (PyCFunction) readDRX,    METH_VARARGS,               readDRX_doc   }, 
	{"readDRSpec", (PyCFunction) readDRSpec, METH_VARARGS,               readDRSpec_doc},
	{"readVDIF",   (PyCFunction) readVDIF,   METH_VARARGS|METH_KEYWORDS, readVDIF_doc  }, 
	{NULL, NULL, 0, NULL}
};

PyDoc_STRVAR(GoFast_doc, "Go Fast! (TM) - TBW, TBN, DRX, DR Spectrometer, and VDIF readers written in C");


/*
  Module Setup - Initialization
*/

PyMODINIT_FUNC init_gofast(void) {
	PyObject *m, *dict1, *dict2;
	
	// Initialize the look-up tables
	initLWALUTs();
	initVDIFLUTs();
	
	// Module definitions and functions
	m = Py_InitModule3("_gofast", GoFastMethods, GoFast_doc);
	import_array();

	// Exceptions
	
	//   1.  syncError -> similar to lsl.reader.errors.syncError
	dict1 = (PyObject *) PyDict_New();
	if(dict1 == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create exception dictionary");
		Py_XDECREF(dict1);
		Py_XDECREF(m);
	}
	PyDict_SetItemString(dict1, "__doc__", \
		PyString_FromString("Exception raised when a reader encounters an error with one or more of the four sync. words."));
	syncError = PyErr_NewException("_gofast.syncError", PyExc_IOError, dict1);
	Py_INCREF(syncError);
	PyModule_AddObject(m, "syncError", syncError);
	
	//    2. eofError -> similar to lsl.reader.errors.eofError
	dict2 = (PyObject *) PyDict_New();
	if(dict2 == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create exception dictionary");
		Py_XDECREF(dict1);
		Py_XDECREF(syncError);
		Py_XDECREF(dict2);
		Py_XDECREF(m);
	}
	PyDict_SetItemString(dict2, "__doc__", \
		PyString_FromString("Exception raised when a reader encounters the end-of-file while reading."));
	eofError = PyErr_NewException("_gofast.eofError", PyExc_IOError, dict2);
	Py_INCREF(eofError);
	PyModule_AddObject(m, "eofError", eofError);
	
	// Version and revision information
	PyModule_AddObject(m, "__version__", PyString_FromString("0.7"));
	PyModule_AddObject(m, "__revision__", PyString_FromString("$Rev$"));
	
}

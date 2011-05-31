#include "Python.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/mman.h>
#include "Defines.h"
#include "Disk.h"
#include "FileSystem.h"

#define _CHUNKSIZE_ fs->fileSystemHeaderData->filesystemInfo.chunkSize


/*
  Function to trim the leading and trailing whitespace characters from a string.

  From:
    http://stackoverflow.com/questions/122616/painless-way-to-trim-leading-trailing-whitespace-in-c
*/

char *trim(char *str) {
	size_t len = 0;
	char *frontp = str - 1;
	char *endp = NULL;

	if( str == NULL ) {
		return NULL;
	}

	if( str[0] == '\0' ) {
		return str;
	}

	len = strlen(str);
	endp = str + len;

	/*
	Move the front and back pointers to address
	the first non-whitespace characters from
	each end.
	*/
	while( isspace(*(++frontp)) );
	while( isspace(*(--endp)) && endp != frontp );

	if( str + len - 1 != endp ) {
		*(endp + 1) = '\0';
	} else { 
		if( frontp != str &&  endp == frontp ) {
			*str = '\0';
		}
	}

	/* 
	Shift the string so that it starts at str so
	that if it's dynamically allocated, we can
	still free it on the returned pointer.  Note
	the reuse of endp to mean the front of the
	string buffer now.
	*/
	endp = str;
	if( frontp != str ) {
		while( *frontp ) *endp++ = *frontp++;
		*endp = '\0';
	}

	return str;
}


/*
  Function to query a specified device and (try to) open the LWAFS filesystem
  on the device.  The function takes in a device name string, a malloc'd Disk
  and a malloc'd FileSystem.  If everything is successful, Disk and FileSystem 
  are set and a status code of 1 is returned.  If any step in the process 
  fails, 0 is returned
  
  This function wraps the calls to Disk_IdentifyAll, Disk_GetDiskInfo, and 
  FileSystem_Open, all from the DROS software.
*/

int readyDevice(char* device, Disk* disk, FileSystem* fs) {
	StatusCode out;
	
	out = Disk_IdentifyAll();
	if( out != SUCCESS ) {
		return 0;
	}
	
	out = Disk_GetDiskInfo(device, disk);
	if( out != SUCCESS || disk == NULL ) {
		return 0;
	}
	
	out = FileSystem_Open(fs, disk);
	if( out != SUCCESS || fs == NULL ) {
		return 0;
	}
	
	return 1;
}


/*
  Function to tear down what readyDevice has setup.  Takes the same set of
  inputs and closes the filesystem.  If all went well, a status code of 1 is
  returned, otherwise 0 is returned.
  
  This function wraps the call to FileSystem_Close from the DROS software.
*/

int closeDevice(char* device, Disk* disk, FileSystem* fs) {
	StatusCode out;
	
	out = FileSystem_Close(fs);
	if( out != SUCCESS ) {
		return 0;
	}
	
	return 1;
}


static PyObject *listDevice(PyObject *self, PyObject *args) {
	PyObject *files;
	char *device;
	int index;

	if(!PyArg_ParseTuple(args, "s", &device)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Setup the various DROS variables
	Disk* disk = NULL;
	FileSystem* fs = NULL;
	
	disk = (Disk*) malloc(sizeof(Disk));
	fs = (FileSystem*) malloc(sizeof(FileSystem));
	memset((void*) fs, 0,  sizeof(FileSystem));
	
	// Query and open
	if(! readyDevice(device, disk, fs) ) {
		PyErr_Format(PyExc_IOError, "Cannot ready device");
		free(fs);
		free(disk);
		return NULL;
	}
	
	// List files
	
	files = (PyObject *) PyList_New(fs->fileSystemHeaderData->filesystemInfo.number_files_present);
	if( files == NULL ) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output filename list");
		Py_XDECREF(files);
		return NULL;
	}
	
	index = 1;
	while(index < fs->fileSystemHeaderData->filesystemInfo.number_files_present+1) {
		PyList_SetItem(files, (Py_ssize_t) index-1, Py_BuildValue("s", fs->fileSystemHeaderData->fileInfo[index].name));
		index++;
	}
	
	closeDevice(device, disk, fs);
	free(fs);
	free(disk);
	
	return files;
}

PyDoc_STRVAR(listDevice_doc, \
"Function to return a list of all files on a DRSU specified by device and return them\n\
as a list of filenames.");


static PyObject *getDeviceChunkSize(PyObject *self, PyObject *args) {
	PyLongObject *chunksize;
	char *device;

	if(!PyArg_ParseTuple(args, "s", &device)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}

	// Setup the various DROS variables
	Disk* disk = NULL;
	FileSystem* fs = NULL;
	
	disk = (Disk*) malloc(sizeof(Disk));
	fs = (FileSystem*) malloc(sizeof(FileSystem));
	memset((void*) fs, 0,  sizeof(FileSystem));
	
	// Query and open
	if(! readyDevice(device, disk, fs) ) {
		PyErr_Format(PyExc_IOError, "Cannot ready device");
		free(fs);
		free(disk);
		return NULL;
	}
	
	// Get chunksize
	chunksize = (PyLongObject *) Py_BuildValue("l", _CHUNKSIZE_);
	
	closeDevice(device, disk, fs);
	free(fs);
	free(disk);
	
	return chunksize;
}

PyDoc_STRVAR(getDeviceChunkSize_doc, \
"Function to return a long integer of the chunk size used on the specified DRSU device.");


static PyObject *getFileSize(PyObject *self, PyObject *args) {
	PyLongObject *filesize;
	char *device, *filename;
	int index;

	if(!PyArg_ParseTuple(args, "ss", &device, &filename)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	// Setup the various DROS variables
	Disk* disk = NULL;
	FileSystem* fs = NULL;
	
	disk = (Disk*) malloc(sizeof(Disk));
	fs = (FileSystem*) malloc(sizeof(FileSystem));
	memset((void*) fs, 0,  sizeof(FileSystem));
	
	// Query and open
	if(! readyDevice(device, disk, fs) ) {
		PyErr_Format(PyExc_IOError, "Cannot ready device");
		free(fs);
		free(disk);
		return NULL;
	}
	
	// Check if the file exists
	if(! FileSystem_FileExists(fs, filename) ) {
		PyErr_Format(PyExc_IOError, "File does not exist");
		free(fs);
		free(disk);
		return NULL;
	}
	
	// Get the file index
	index = _FileSystem_GetFileIndex(fs, filename);
	
	// Get the file size
	filesize = (PyLongObject *) Py_BuildValue("l", fs->fileSystemHeaderData->fileInfo[index].size);
	
	closeDevice(device, disk, fs);
	free(fs);
	free(disk);
	
	return filesize;
}

PyDoc_STRVAR(getFileSize_doc, \
"Function to return a long integer of the file size for a particular file on\n\
the specified DRSU device.  An IOError is raised if the file cannot be found\n\
on the device.");


static PyObject *getFileStart(PyObject *self, PyObject *args) {
	PyLongObject *filestart;
	char *device, *filename;
	int index;

	if(!PyArg_ParseTuple(args, "ss", &device, &filename)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	// Setup the various DROS variables
	Disk* disk = NULL;
	FileSystem* fs = NULL;
	
	disk = (Disk*) malloc(sizeof(Disk));
	fs = (FileSystem*) malloc(sizeof(FileSystem));
	memset((void*) fs, 0,  sizeof(FileSystem));
	
	// Query and open
	if(! readyDevice(device, disk, fs) ) {
		PyErr_Format(PyExc_IOError, "Cannot ready device");
		free(fs);
		free(disk);
		return NULL;
	}
	
	// Check if the file exists
	if(! FileSystem_FileExists(fs, filename) ) {
		PyErr_Format(PyExc_IOError, "File does not exist");
		free(fs);
		free(disk);
		return NULL;
	}
	
	// Get the file index
	index = _FileSystem_GetFileIndex(fs, filename);
	
	// Get the file start in bytes from the beginning of the disk
	filestart = (PyLongObject *) Py_BuildValue("l", fs->fileSystemHeaderData->fileInfo[index].startPosition + _CHUNKSIZE_);
	
	closeDevice(device, disk, fs);
	free(fs);
	free(disk);
	
	return filestart;
}

PyDoc_STRVAR(getFileStart_doc, \
"Function to return a long integer of the file start position in bytes \n\
(corrected for the file start pattern) for a particular file on the specified\n\
DRSU device.  An IOError is raised if the file cannot be found on the device.");


static PyObject *getFileType(PyObject *self, PyObject *args) {
	PyStringObject *filetype;
	char *device, *filename;
	int index;

	if(!PyArg_ParseTuple(args, "ss", &device, &filename)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	// Setup the various DROS variables
	Disk* disk = NULL;
	FileSystem* fs = NULL;
	
	disk = (Disk*) malloc(sizeof(Disk));
	fs = (FileSystem*) malloc(sizeof(FileSystem));
	memset((void*) fs, 0,  sizeof(FileSystem));
	
	// Query and open
	if(! readyDevice(device, disk, fs) ) {
		PyErr_Format(PyExc_IOError, "Cannot ready device");
		free(fs);
		free(disk);
		return NULL;
	}
	
	// Check if the file exists
	if(! FileSystem_FileExists(fs, filename) ) {
		PyErr_Format(PyExc_IOError, "File does not exist");
		free(fs);
		free(disk);
		return NULL;
	}
	
	// Get the file index
	index = _FileSystem_GetFileIndex(fs, filename);
	
	// Get the file start in bytes from the beginning of the disk
	filetype = (PyStringObject *) Py_BuildValue("s", trim(fs->fileSystemHeaderData->fileInfo[index].metaData.format));
	
	closeDevice(device, disk, fs);
	free(fs);
	free(disk);
	
	return filetype;
}

PyDoc_STRVAR(getFileType_doc, \
"Function to return a string of the file meta-data type (DEFAULT_DRX, etc.) for\n\
a particular file on the specified DRSU device.  An IOError is raised if the\n\
file cannot be found on the device.");


static PyObject *getFileTime(PyObject *self, PyObject *args) {
	PyLongObject *filemtime;
	char *device, *filename;
	int index;

	if(!PyArg_ParseTuple(args, "ss", &device, &filename)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	// Setup the various DROS variables
	Disk* disk = NULL;
	FileSystem* fs = NULL;
	
	disk = (Disk*) malloc(sizeof(Disk));
	fs = (FileSystem*) malloc(sizeof(FileSystem));
	memset((void*) fs, 0,  sizeof(FileSystem));
	
	// Query and open
	if(! readyDevice(device, disk, fs) ) {
		PyErr_Format(PyExc_IOError, "Cannot ready device");
		free(fs);
		free(disk);
		return NULL;
	}
	
	// Check if the file exists
	if(! FileSystem_FileExists(fs, filename) ) {
		PyErr_Format(PyExc_IOError, "File does not exist");
		free(fs);
		free(disk);
		return NULL;
	}
	
	// Get the file index
	index = _FileSystem_GetFileIndex(fs, filename);
	
	// Get the file creation time
	filemtime = (PyLongObject *) Py_BuildValue("l", fs->fileSystemHeaderData->fileInfo[index].mtime.tv_sec);
	
	closeDevice(device, disk, fs);
	free(fs);
	free(disk);
	
	return filemtime;
}

PyDoc_STRVAR(getFileTime_doc, \
"Function to return a long integer of the file modification time for a\n\
particular file on the specified DRSU device.  An IOError is raised if the\n\
file cannot be found on the device.");


static PyObject *listFiles(PyObject *self, PyObject *args) {
	PyObject *fileClass, *file, *files;
	char *device;
	int index;

	if( !PyArg_ParseTuple(args, "sO", &device, &fileClass) ) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	// Make sure we are passed a class, not an instance
	if( !PyType_Check(fileClass) ) {
		PyErr_Format(PyExc_RuntimeError, "Must pass a File class, not an instance");
		return NULL;
	}

	// Setup the various DROS variables
	Disk* disk = NULL;
	FileSystem* fs = NULL;
	
	disk = (Disk*) malloc(sizeof(Disk));
	fs = (FileSystem*) malloc(sizeof(FileSystem));
	memset((void*) fs, 0,  sizeof(FileSystem));
	
	// Query and open
	if(! readyDevice(device, disk, fs) ) {
		PyErr_Format(PyExc_IOError, "Cannot ready device");
		free(fs);
		free(disk);
		return NULL;
	}
	
	// List files
	files = (PyObject *) PyList_New(fs->fileSystemHeaderData->filesystemInfo.number_files_present);
	if( files == NULL ) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output File list");
		Py_XDECREF(files);
		return NULL;
	}
	
	index = 1;
	while(index < fs->fileSystemHeaderData->filesystemInfo.number_files_present+1) {
		file = PyObject_CallObject(fileClass, Py_BuildValue("(ssl)", device, fs->fileSystemHeaderData->fileInfo[index].name, fs->fileSystemHeaderData->fileInfo[index].size));
		
		// Set file times
		PyObject_SetAttrString(file, "ctime", Py_BuildValue("l", fs->fileSystemHeaderData->fileInfo[index].ctime.tv_sec));
		PyObject_SetAttrString(file, "mtime", Py_BuildValue("l", fs->fileSystemHeaderData->fileInfo[index].mtime.tv_sec));
		PyObject_SetAttrString(file, "atime", Py_BuildValue("l", fs->fileSystemHeaderData->fileInfo[index].atime.tv_sec));
		
		// Set file start/stop/chunk sizes
		PyObject_SetAttrString(file, "start", Py_BuildValue("l", fs->fileSystemHeaderData->fileInfo[index].startPosition + _CHUNKSIZE_));
		PyObject_SetAttrString(file, "stop", Py_BuildValue("l", fs->fileSystemHeaderData->fileInfo[index].stopPosition - _CHUNKSIZE_));
		PyObject_SetAttrString(file, "chunkSize", Py_BuildValue("l", _CHUNKSIZE_));
		
		// Meta data
		PyObject_SetAttrString(file, "mjd", PyLong_FromString(fs->fileSystemHeaderData->fileInfo[index].metaData.StartMJD, NULL, 10));
		PyObject_SetAttrString(file, "mpm", PyLong_FromString(fs->fileSystemHeaderData->fileInfo[index].metaData.StartMPM, NULL, 10));
		PyObject_SetAttrString(file, "mode", Py_BuildValue("s", trim(fs->fileSystemHeaderData->fileInfo[index].metaData.format)));
		
		PyList_SetItem(files, (Py_ssize_t) index-1, file);
		index++;
	}
	
	closeDevice(device, disk, fs);
	free(fs);
	free(disk);
	
	return files;
}

PyDoc_STRVAR(listFiles_doc, \
"Function to list all of the files on a DRSU specified by device and return them\n\
as a list of 'File' instances.  These instances can then be used to provide direct\n\
access to the files stored on the DRSU.\n\
\n\
The use of this function is preferred over listDevice/getFileSize/etc. because it\n\
stores all of the information returned in the other function in the File structure\n\
along with additional information.  It also makes it 'easier' to open the specified\n\
file for reading.");


static PyMethodDef DRSUMethods[] = {
	{"listDevice",         listDevice,         METH_VARARGS, listDevice_doc}, 
	{"getDeviceChunkSize", getDeviceChunkSize, METH_VARARGS, getDeviceChunkSize_doc}, 
	{"getFileSize",        getFileSize,        METH_VARARGS, getFileSize_doc}, 
	{"getFileStart",       getFileStart,       METH_VARARGS, getFileStart_doc}, 
	{"getFileType",        getFileType,        METH_VARARGS, getFileType_doc}, 
	{"getFileTime",        getFileTime,        METH_VARARGS, getFileTime_doc}, 
	{"listFiles",          listFiles,          METH_VARARGS, listFiles_doc}, 
	{NULL,                 NULL,               0,            NULL}
};

PyDoc_STRVAR(DRSUMethods_doc, \
"Module to provide the information necessary for direct access to files stored\n\
on a data recorder storage unit (DRSU).\n\
\n\
.. warning::\n\
\tThis module is currently in an experimental phase.\n\
");


PyMODINIT_FUNC init_drsu(void) {
	PyObject *m;

	// Module definitions and functions
	m = Py_InitModule3("_drsu", DRSUMethods, DRSUMethods_doc);
}

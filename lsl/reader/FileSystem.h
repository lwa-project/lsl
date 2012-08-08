/*
  MCS-DROS: Monitor and Control System - Data Recorder Operating Software
  Copyright (C) 2009-2010  Virginia Tech, Christopher Wolfe <chwolfe2@vt.edu>

  This file is part of MCS-DROS.

  MCS-DROS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  MCS-DROS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with MCS-DROS.  If not, see <http://www.gnu.org/licenses/>.
*/
/*
 * FileSystem.h
 *
 *  Created on: Oct 25, 2009
 *      Author: chwolfe2
 */
#ifndef FILESYSTEM_H_
#define FILESYSTEM_H_
#include "Defines.h"
#include "Disk.h"
#include "Time.h"
#include <aio.h>


typedef int FILEMODE;

#define READ 							0
#define WRITE 							1
#define FILE_HEADER_SIZE 				4096
//#define FILE_HEADER_RESERVED_SPACE      2944
#define FILE_MAX_NAME_LENGTH 		    1023
#define MAX_NUMBER_FILES 			    8191
#define DEFAULT_VOLUME_NAME				"MCS-DRSU-Storage"


typedef struct __FileMetaData{
	char		tag[17];
	char		StartMJD[7];
	char		StartMPM[10];
	char		StopMJD[7];
	char		StopMPM[10];
	char		format[33];
	char		size[16];
	char		usage[16];
	char		complete[4];
} FileMetaData;

typedef struct __FileInfo{
	union {
		struct {
			size_t startPosition;
			size_t stopPosition;
			size_t size;
			size_t isOpen;
			size_t currentPosition;
			size_t tags;
			FILEMODE mode;
			TimeStamp creationDate;
			uid_t uid;
			gid_t gid;
			struct timespec ctime;
			struct timespec atime;
			struct timespec mtime;
			mode_t amode;
			char name[FILE_MAX_NAME_LENGTH+1];
			char reserved[1];
			FileMetaData metaData;
		};
		char dataptr [4096];
	};

} FileInfo;

typedef struct __FileSystemInfo{
	union {
		struct {
			size_t diskSize;				// size of the disk
			size_t chunkSize;				// chunk Size
			size_t numDrives;				// number of drives
			size_t maxFiles;				// max num files
			size_t lastFileIndex;			// last opened file index
			size_t number_files_present;    // number of actual files contained by the fs
			FILEMODE dummy;					// space holder to keep FileInfo and FileSystemInfo same size
			TimeStamp creationDate;			// filesystem creation date
			uid_t uid;
			gid_t gid;
			struct timespec ctime;
			struct timespec atime;
			struct timespec mtime;
			mode_t amode;
			char volname[FILE_MAX_NAME_LENGTH+1]; // volume name
			char reserved[1];
		};
		char dataptr [4096];
	};
} FileSystemInfo;
typedef union __FileSystemBlock{
		FileInfo        fileInfo [MAX_NUMBER_FILES+1];
		FileSystemInfo  filesystemInfo;
}FileSystemBlock;

typedef struct __FileSystem{
	int 				rawDeviceHandle;
	Disk*				diskInfo;
	boolean 			filesystemIsOpen;
	boolean 			filesystemInSync;
	void* 				IMFS_ptr;  				 // pointer
	size_t  			IMFS_Size; 				 // size in bytes needed to hold the filesystem info...
	boolean 			closeIsRegistered; // indicates whether a request to close on exit has been installed
	FileSystemBlock * 	fileSystemHeaderData;
} FileSystem;

typedef struct __MyFile{
	union {
		struct {
			int    			index;
			FILEMODE 		mode;
			boolean 		isOpen;
			boolean 		reof;
			boolean 		lastRPacket;
			boolean 		lastWPacket;
			size_t  		lastRPacketRealSize;
			size_t  		lastWPacketRealSize;
			FileSystem * 	fs;
			struct aiocb 	writeOpInfo;
			struct aiocb 	readOpInfo;
			char * 			buffer;				// a buffer pointer for use with the read only fuse module
			size_t			bufferSize; 		// how big is the buffer
			size_t			buffer1stByte;		// where in the file would the buffer line up with
			size_t			contentSize;		// how many bytes are actually in the buffer (typically buffersize, except EOF)
			size_t			consumed;			// how many bytes have been read already
		};
		char rawdata [1024];
	};
} MyFile;


//#define fs (*((FileSystem*)fsbuf))
StatusCode FileSystem_Create(Disk * disk);
boolean    FileSystem_TestDiskForLWAFS(Disk * disk);
StatusCode FileSystem_Verify(FileSystem * fs);
StatusCode FileSystem_Open(FileSystem* fs, Disk * disk);
StatusCode FileSystem_Close(FileSystem * fs);
StatusCode FileSystem_ListFiles(FileSystem * fs);
size_t     FileSystem_GetFileCount(FileSystem * fs);
void       FileSystem_GetFileMetaData(FileSystem * fs, int index, FileMetaData** data);
StatusCode FileSystem_GetFileFormat(FileSystem * fs, int index,char * buffer);
void       FileSystem_SetFileMetaData(FileSystem * fs, int index, size_t startMPM, size_t startMJD, size_t reference, size_t StopMPM, size_t StopMJD, char * format);
void 	   FileSystem_SetFileMetaDataIsComplete(FileSystem * fs, int index, int complete);

boolean    FileSystem_IsFileCreatable(FileSystem * fs, char* name, size_t size, boolean overwrite);

StatusCode FileSystem_AsynchronousWriteStart(MyFile* file, char * srcPtr, size_t bytes);
StatusCode FileSystem_IsAsynchronousWriteDone(MyFile* file);
StatusCode FileSystem_GetAsynchronousBytesWritten(MyFile* file, size_t* result);

size_t FileSystem_Seek(MyFile* file, size_t offset);
StatusCode FileSystem_AsynchronousReadStart(MyFile* file, char * dstPtr, size_t bytes);
StatusCode FileSystem_IsAsynchronousReadDone(MyFile* file);
StatusCode FileSystem_GetAsynchronousBytesRead(MyFile* file, size_t *result);

StatusCode FileSystem_CreateFile(FileSystem * fs, char* name, size_t size, boolean overwrite, int* fileIndex);
StatusCode FileSystem_OpenFile(FileSystem * fs, MyFile* file,const char* fileName, FILEMODE mode);
StatusCode FileSystem_CloseFile(MyFile* file);
StatusCode FileSystem_DeleteFile(FileSystem * fs, int fileIndex);
boolean    FileSystem_FileExists(FileSystem * fs,const char* fileName);

StatusCode _FileSystem_SynchronousWrite(FileSystem * fs, void* buf, off_t position, size_t count, ssize_t * result);
StatusCode _FileSystem_SynchronousRead(FileSystem * fs, void* buf, off_t position, size_t count, ssize_t * result);
StatusCode _FileSystem_AllocateIMFS(FileSystem * fs);
StatusCode _FileSystem_DeallocateIMFS(FileSystem * fs);

int        _FileSystem_GetFileIndex(FileSystem * fs, const char* name);
void 	   _FileSystem_RegisterActiveFilesystem(FileSystem * fs);
void 	   _FileSystem_UnregisterActiveFilesystem(FileSystem * fs);
void 	   _FiLeSystem_CloseAllActiveSystems(void);
boolean    _FileSystem_CheckPattern(void* ptr, size_t size);
void 	   _FileSystem_CreatePattern(void* ptr, size_t size);
size_t 	   _FileSystem_RangeConflicts(FileSystem * fs, size_t start, size_t stop);
boolean    _FileSystem_RegionOverlaps(size_t starta, size_t stopa, size_t startb, size_t stopb);
boolean    _FileSystem_IsBetween(size_t a, size_t b, size_t c);
void 	   _FileSystem_DumpFileInfo(FileSystem * fs, int index);
size_t FileSystem_GetFreeSpace(FileSystem * fs);

#endif /* FILESYSTEM_H_ */

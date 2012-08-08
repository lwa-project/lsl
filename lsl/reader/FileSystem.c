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
 * FileSystem.c
 *
 *  Created on: Oct 25, 2009
 *      Author: chwolfe2
 */

#include <stdio.h>
#include <stdlib.h>
#include	<sys/stat.h>
#include	<unistd.h>
#include	<fcntl.h>
#include	<errno.h>
#include 	<string.h>
#include    <assert.h>
#include <sys/mman.h>
#include <linux/fs.h>
#include <sys/ioctl.h>
#include <ctype.h>
#include	<aio.h>
#include <time.h>
#include <sys/time.h>

#include "FileSystem.h"
//#include "Globals.h"
#include "Log.h"
#include "ConsoleColors.h"
#define isAligned(offset) ((((size_t)offset) & ~(512l -1)) == ((size_t)offset))

FileSystem* ActiveFileSystems[100];
int 		NumActiveFilesystems=0;

void _FileSystem_RegisterActiveFilesystem(FileSystem * fs){
	Log_Add("[FILE SYSTEM] Registering Active File System %d:<%lu>",NumActiveFilesystems ,(size_t)fs);
	ActiveFileSystems[NumActiveFilesystems++]=fs;
}
void _FileSystem_UnregisterActiveFilesystem(FileSystem * fs){
	int i=0;
	int j=0;
	while (i<NumActiveFilesystems){
		if (ActiveFileSystems[i]==fs){
			Log_Add("[FILE SYSTEM] Unregistering Active File System %d:<%lu>",i ,(size_t)fs);
			for (j=i;j<(NumActiveFilesystems-1);j++){
				Log_Add("[FILE SYSTEM] Copying record %d to %d",j+1,j);
				ActiveFileSystems[j]=ActiveFileSystems[j+1];
			}
			NumActiveFilesystems--;
		}
	}
}
void _FiLeSystem_CloseAllActiveSystems(void){
	while (NumActiveFilesystems>0){
		Log_Add("[FILE SYSTEM] FileSystem_Close(ActiveFileSystems[0]);");
		FileSystem_Close(ActiveFileSystems[0]);
		// automatically done by FileSystem_Close()
		//_FileSystem_UnregisterActiveFilesystem(ActiveFileSystems[0]);
	}
}

StatusCode _FileSystem_DeallocateIMFS(FileSystem * fs){
//(void** fsbuf, size_t fsSize){
	if (fs->IMFS_ptr!=NULL){
		int status = munmap(fs->IMFS_ptr, fs->IMFS_Size);
		if (status == -1) {
			return FAILURE;
		}
		fs->IMFS_ptr=NULL;
		return SUCCESS;
	}
	return FAILURE;
}

StatusCode _FileSystem_AllocateIMFS(FileSystem * fs){
	//(void** fsbuf, size_t fsSize){
	assert (fs!=NULL);
	assert (fs->IMFS_ptr==NULL);
	assert (fs->IMFS_Size>=sizeof(FileSystemInfo));
	fs->IMFS_ptr = mmap(NULL, fs->IMFS_Size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_POPULATE | MAP_LOCKED, -1, 0);
	return ( (fs->IMFS_ptr == MAP_FAILED) ?  FAILURE : SUCCESS);
}


StatusCode _FileSystem_SynchronousRead(FileSystem * fs, void* buf, off_t position, size_t count, ssize_t * result){
	//printf("%lu\n",count);

	assert(fs!=NULL);
	assert(fs->rawDeviceHandle!=-1);
	assert(isAligned(buf));
	assert(isAligned(position));
	assert(isAligned(count));


	ssize_t status;
	size_t bytes = 0;
	//printf("hdl:%d, buf:%lu, pos:%lu, cnt:%lu\n", fs->rawDeviceHandle, (size_t)buf, (size_t)position, count);
	lseek(fs->rawDeviceHandle, position, SEEK_SET);
	while (bytes < count){
		status = read(fs->rawDeviceHandle, buf + bytes, count - bytes);
		if (status == -1) {
			if (errno!=EWOULDBLOCK){
				//printf("status=\'%lu\' or '%s'\n", status, strerror(status));
				//perror("read(): ");

				*result=0;
				return FAILURE;
				/*perror("read(...) <synchronous>");
				close(rawDeviceHandle);
				exit(EXIT_CODE_FAIL);*/
			}
			status=0;
		} else if (status == 0) {
			*result=0;
			return FAILURE;/*
			fprintf(stderr, "Reached EOF\n");
			exit(EXIT_CODE_FAIL);*/
		}
		bytes += status;
	}
	*result= bytes;
	return SUCCESS;
}

StatusCode _FileSystem_SynchronousWrite(FileSystem * fs, void* buf, off_t position, size_t count, ssize_t * result){
	ssize_t status;
	size_t bytes = 0;
	assert(isAligned((size_t)buf));
	assert(isAligned((size_t)position));
	assert(isAligned((size_t)count));
	lseek(fs->rawDeviceHandle, position, SEEK_SET);
	while (bytes < count){
		status = write(fs->rawDeviceHandle, buf + bytes, count - bytes);
		if (status == -1) {
			if (errno!=EWOULDBLOCK){
				*result=0;
				return FAILURE;
				/*perror("write(...) <synchronous>");
				close(rawDeviceHandle);
				exit(EXIT_CODE_FAIL);*/
			}
			status=0;
		}
		bytes += status;
	}
	*result= bytes;
	return SUCCESS;
}


boolean FileSystem_TestDiskForLWAFS(Disk * disk){
	assert(disk!=NULL);
	ssize_t bytesRead=0;
	FileSystem  f;
	FileSystem* fs=&f;
	if (disk->type!=DRIVETYPE_INTERNAL_MULTIDISK_DEVICE)
		return 0;
	memset ((void*)fs,0,sizeof(FileSystem));
	fs->rawDeviceHandle = open(disk->deviceName, O_RDONLY );
	if (fs->rawDeviceHandle == -1) {
		Log_Add("[FILE SYSTEM] TestDiskForLWAFS: Warning: could not open disk to test for LWAFS");
		return 0;
	}else{
		fs->diskInfo=disk;
		fs->IMFS_Size = sizeof(FileSystemBlock);
		if (_FileSystem_AllocateIMFS(fs)==FAILURE){
			Log_Add("[FILE SYSTEM] TestDiskForLWAFS: Warning: Allocation failure");
			close(fs->rawDeviceHandle);
			return 0;
		} else {
			memset(fs->IMFS_ptr,0,fs->IMFS_Size);
		}

		if ((_FileSystem_SynchronousRead(fs, fs->IMFS_ptr,0,fs->IMFS_Size, &bytesRead)==FAILURE) ||	(bytesRead!=fs->IMFS_Size)){ // read the whole filesystem in....
			Log_Add("[FILE SYSTEM] TestDiskForLWAFS: Warning: FileSystemSynchronousRead failure");
			_FileSystem_DeallocateIMFS(fs);
			close(fs->rawDeviceHandle);
			return 0;
		}


		fs->fileSystemHeaderData=(FileSystemBlock *)fs->IMFS_ptr;
		if (  (fs->fileSystemHeaderData->filesystemInfo.chunkSize != Disk_GetArrayChunkSize(disk)) ||
			  (fs->fileSystemHeaderData->filesystemInfo.numDrives != Disk_GetArrayNumDisks(disk))  ||
			  (fs->fileSystemHeaderData->filesystemInfo.diskSize  != Disk_GetArraySize(disk))      ||
			  (fs->fileSystemHeaderData->filesystemInfo.maxFiles  != MAX_NUMBER_FILES)             ||
			  (fs->fileSystemHeaderData->filesystemInfo.number_files_present >= fs->fileSystemHeaderData->filesystemInfo.maxFiles)){
			_FileSystem_DeallocateIMFS(fs);
			close(fs->rawDeviceHandle);
			return 0;
		}
		_FileSystem_DeallocateIMFS(fs);
		close(fs->rawDeviceHandle);
		return 1;
	}
}

StatusCode FileSystem_Create(Disk * disk){
	//size_t numFiles      = MAX_NUMBER_FILES;
	struct timeval dt;

	FileSystem fs;
	memset ((void*)&fs,0,sizeof(FileSystem));
	fs.IMFS_ptr=NULL;
	fs.diskInfo=disk;
	fs.closeIsRegistered=false;
	fs.fileSystemHeaderData=NULL;
	fs.filesystemInSync=false;
	fs.filesystemIsOpen=false;
	fs.rawDeviceHandle=-1;
	fs.IMFS_Size=sizeof(FileSystemBlock);
	if (_FileSystem_AllocateIMFS(&fs)!=SUCCESS){
		Log_Add("[FILE SYSTEM] Error: could not allocate IMFS");
		return FAILURE;
	}
	//printf("----- imfs size = %lu\n",fs.IMFS_Size);fflush(stdout);
	fs.fileSystemHeaderData=(FileSystemBlock *)fs.IMFS_ptr;
	fs.rawDeviceHandle = open(disk->deviceName, O_RDWR | O_DIRECT );
	if (fs.rawDeviceHandle == -1) {
		Log_Add("[FILE SYSTEM] Error: could not open device");
		perror("open(...)");
		return FAILURE;
	}

	size_t diskSize          = Disk_GetArraySize(disk);
	size_t chunkSize         = Disk_GetArrayChunkSize(disk);
	size_t numDrives         = Disk_GetArrayNumDisks(disk);
	size_t maxFiles          = MAX_NUMBER_FILES;
	TimeStamp creationDate   = getTimestamp();

	Log_Add("[FILE SYSTEM] WARNING: Creating filesystem on device %s !!!!",disk->deviceName);
	Log_Add("[FILE SYSTEM] %s: "
			"\n                   [FILE SYSTEM] \tVolume name: %s"
			"\n                   [FILE SYSTEM] \t%18lu chunk size"
			"\n                   [FILE SYSTEM] \t%18lu drives "
			"\n                   [FILE SYSTEM] \t%18lu total bytes "
			"\n                   [FILE SYSTEM] \t%18lu maximum number files"
			"\n                   [FILE SYSTEM] \t%18lu current number files"
			"\n                   [FILE SYSTEM] \t%18lu maximum filename length"
			"\n                   [FILE SYSTEM] \t%18lu bytes reserved for formatting",
			disk->deviceName,
			DEFAULT_VOLUME_NAME,
			chunkSize,
			numDrives,
			diskSize,
			maxFiles,
			0l,
			(unsigned long) FILE_MAX_NAME_LENGTH,
			sizeof(FileSystemBlock)
	);

	// zero the record
	memset(fs.IMFS_ptr,0, fs.IMFS_Size);
	fs.fileSystemHeaderData->filesystemInfo.chunkSize			 = chunkSize;
	fs.fileSystemHeaderData->filesystemInfo.creationDate		 = creationDate;
	fs.fileSystemHeaderData->filesystemInfo.diskSize			 = diskSize;
	fs.fileSystemHeaderData->filesystemInfo.numDrives			 = numDrives;
	fs.fileSystemHeaderData->filesystemInfo.maxFiles			 = maxFiles;
	fs.fileSystemHeaderData->filesystemInfo.number_files_present = 0;
	fs.fileSystemHeaderData->filesystemInfo.lastFileIndex 	     = 0;
	strcpy(fs.fileSystemHeaderData->filesystemInfo.volname,DEFAULT_VOLUME_NAME);
	fs.fileSystemHeaderData->filesystemInfo.uid = getuid();
	fs.fileSystemHeaderData->filesystemInfo.gid = getgid();
	gettimeofday(&dt, NULL);
	fs.fileSystemHeaderData->filesystemInfo.ctime.tv_nsec = 1000l * dt.tv_usec;
	fs.fileSystemHeaderData->filesystemInfo.ctime.tv_sec  = dt.tv_sec;
	memcpy((void*)&(fs.fileSystemHeaderData->filesystemInfo.mtime), (void*)&(fs.fileSystemHeaderData->filesystemInfo.ctime), sizeof(struct timespec));
	memcpy((void*)&(fs.fileSystemHeaderData->filesystemInfo.atime), (void*)&(fs.fileSystemHeaderData->filesystemInfo.ctime), sizeof(struct timespec));
	fs.fileSystemHeaderData->filesystemInfo.amode = S_IFDIR | S_IWGRP | S_IWOTH | S_IWUSR | S_IRGRP | S_IROTH | S_IRUSR;
	// write the whole filesystem structure out...
	ssize_t bytesWritten;
	if ((_FileSystem_SynchronousWrite(&fs, fs.IMFS_ptr, 0, fs.IMFS_Size, &bytesWritten)!=SUCCESS) || bytesWritten!=fs.IMFS_Size){
		Log_Add("[FILE SYSTEM] FATAL error: couldn't write filesystem info to disk. data may be lost.");
		close(fs.rawDeviceHandle);
		munmap(fs.IMFS_ptr,fs.IMFS_Size);
		return FAILURE;
	}
	close(fs.rawDeviceHandle);
	munmap(fs.IMFS_ptr,fs.IMFS_Size);
	return SUCCESS;
}

StatusCode FileSystem_Close(FileSystem * fs){
	StatusCode s=FAILURE;
	Log_Add("[FILE SYSTEM] Closing filesystem.");
	if (fs->rawDeviceHandle!=-1){
		ssize_t bytesWritten;
		if ((_FileSystem_SynchronousWrite(fs,fs->IMFS_ptr,0,fs->IMFS_Size,&bytesWritten)!=SUCCESS) || (bytesWritten!=fs->IMFS_Size)){// write the whole filesystem structure out...
			Log_Add("[FILE SYSTEM] FATAL error: couldn't write filesystem info to disk. data may be lost.");
			s = FAILURE;
		} else {
			s = SUCCESS;
		}
		close(fs->rawDeviceHandle);
		fs->rawDeviceHandle=-1;
		munmap(fs->IMFS_ptr,fs->IMFS_Size);
		fs->IMFS_ptr=NULL;
		_FileSystem_UnregisterActiveFilesystem(fs);
		fs->closeIsRegistered=false;
		fs->fileSystemHeaderData=NULL;
		fs->diskInfo=NULL;
		fs->filesystemInSync=false;
		fs->filesystemIsOpen=false;
	}
	return s;
}
void _FileSystem_DumpFileInfo(FileSystem * fs, int index){
	char buf[FILE_MAX_NAME_LENGTH];
	int i=0;
	while(isprint(fs->fileSystemHeaderData->fileInfo[index].name[i]) && i<FILE_MAX_NAME_LENGTH){
		buf[i]=fs->fileSystemHeaderData->fileInfo[index].name[i];
		i++;
	}
	buf[i]='\0';
	char buf2[100];
	sprintf(buf2,"[FILE SYSTEM] DumpFileInfo: %.2lu/%.2lu/%.4lu %.2lu:%.2lu:%.2lu",
			fs->fileSystemHeaderData->fileInfo[index].creationDate.day,
			fs->fileSystemHeaderData->fileInfo[index].creationDate.month,
			fs->fileSystemHeaderData->fileInfo[index].creationDate.year,
			fs->fileSystemHeaderData->fileInfo[index].creationDate.hour,
			fs->fileSystemHeaderData->fileInfo[index].creationDate.minute,
			fs->fileSystemHeaderData->fileInfo[index].creationDate.second);

	printf( "[FILE SYSTEM] File information for index %d:\n"
			"[FILE SYSTEM]   Name:                %s \n"
			"[FILE SYSTEM]   Creation date:       %s \n"
			"[FILE SYSTEM]   Start position:      %lu\n"
			"[FILE SYSTEM]   Stop position:       %lu\n"
			"[FILE SYSTEM]   Size:                %lu bytes\n"
			"[FILE SYSTEM]   Improperly closed:   %s\n"
			"[FILE SYSTEM]   Last write position: %lu\n"
			"[FILE SYSTEM]   Tags                 %lu\n",
			index,
			buf,
			buf2,
			fs->fileSystemHeaderData->fileInfo[index].startPosition,
			fs->fileSystemHeaderData->fileInfo[index].stopPosition,
			fs->fileSystemHeaderData->fileInfo[index].size,
			((fs->fileSystemHeaderData->fileInfo[index].isOpen)?"yes":"no"),
			fs->fileSystemHeaderData->fileInfo[index].currentPosition,
			fs->fileSystemHeaderData->fileInfo[index].tags);

}
#define SilentVerify 1

StatusCode FileSystem_Verify(FileSystem * fs){
	Log_Add("[FILE SYSTEM] Checking file system integrity...");
	size_t i=0;
	size_t foundcount=0;
	void * patternBuff=mmap(NULL, fs->fileSystemHeaderData->filesystemInfo.chunkSize, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_POPULATE | MAP_LOCKED, -1, 0);
	if ((patternBuff==NULL)||(patternBuff==MAP_FAILED)){
		Log_Add("[FILE SYSTEM] Error: can't allocate space for filesystem verification.");
		perror("mmap(...)");
		close(fs->rawDeviceHandle);
		exit(EXIT_CODE_FAIL);
	}
	if (!SilentVerify) Log_Add("[FILE SYSTEM] Header says there are %lu out of %lu possible files.",fs->fileSystemHeaderData->filesystemInfo.number_files_present,fs->fileSystemHeaderData->filesystemInfo.maxFiles);
	for(i=1;i<fs->fileSystemHeaderData->filesystemInfo.maxFiles;i++){
		if (!SilentVerify) Log_Add("[FILE SYSTEM] Checking index %.8lu of %.8lu\r",i,fs->fileSystemHeaderData->filesystemInfo.maxFiles);

		// is file in use ?
		if (fs->fileSystemHeaderData->fileInfo[i].startPosition!=0){
			foundcount++;
			//printf("Possible file at index %luChecking consistency.",i);

			// check to see that the file starts at a reasonable location
			if (fs->fileSystemHeaderData->fileInfo[i].startPosition<fs->IMFS_Size){
				_FileSystem_DumpFileInfo(fs,i);
				Log_Add("[FILE SYSTEM] File consistency check failed: Start position overlaps filesystem data.");
				Log_Add("[FILE SYSTEM] Deleting File record.");
				memset((void*)&fs->fileSystemHeaderData->fileInfo[i],0,sizeof(FileInfo));
				continue;
			}
			// check to see that the file stops at a reasonable location
			if (fs->fileSystemHeaderData->fileInfo[i].stopPosition < (fs->fileSystemHeaderData->fileInfo[i].startPosition+(2*fs->fileSystemHeaderData->filesystemInfo.chunkSize))){
				_FileSystem_DumpFileInfo(fs,i);
				Log_Add("[FILE SYSTEM] File consistency check failed: End position is less than the minimum.");
				Log_Add("[FILE SYSTEM] Deleting File record.");
				memset((void*)&fs->fileSystemHeaderData->fileInfo[i],0,sizeof(FileInfo));
				continue;
			}
			// check to see that the file's reported size is within the boundaries of the allocated space
			size_t allocatedSize=fs->fileSystemHeaderData->fileInfo[i].stopPosition - fs->fileSystemHeaderData->fileInfo[i].startPosition;
			if (fs->fileSystemHeaderData->fileInfo[i].size > (allocatedSize-2*fs->fileSystemHeaderData->filesystemInfo.chunkSize)){
				_FileSystem_DumpFileInfo(fs,i);
				Log_Add("[FILE SYSTEM] File consistency check failed: Length exceeds allocated space.");
				Log_Add("[FILE SYSTEM] Truncating file to fit allocated space.");
				fs->fileSystemHeaderData->fileInfo[i].size=allocatedSize-2*fs->fileSystemHeaderData->filesystemInfo.chunkSize;
				continue;
			}
			// make sure the file would actually fit on the array
			if (allocatedSize>(fs->fileSystemHeaderData->filesystemInfo.diskSize-fs->IMFS_Size)){
				_FileSystem_DumpFileInfo(fs,i);
				Log_Add("[FILE SYSTEM] File consistency check failed: Length exceeds volume size.");
				Log_Add("[FILE SYSTEM] Deleting File record.");
				memset((void*)&fs->fileSystemHeaderData->fileInfo[i],0,sizeof(FileInfo));
				continue;
			}

			// check for the start tag at the specified location
			//printf("reading pattern at %lu\n",fs->fileSystemHeaderData->fileInfo[i].startPosition);
			ssize_t bytesRead=0;
			StatusCode s;
			s=_FileSystem_SynchronousRead(fs,patternBuff,fs->fileSystemHeaderData->fileInfo[i].startPosition,fs->fileSystemHeaderData->filesystemInfo.chunkSize, &bytesRead);
			if (s!=SUCCESS || bytesRead!=fs->fileSystemHeaderData->filesystemInfo.chunkSize){
				Log_Add("[FILE SYSTEM] File consistency check failed: Can't read start tag.");
				continue;
			}
			if (!_FileSystem_CheckPattern(patternBuff, fs->fileSystemHeaderData->filesystemInfo.chunkSize)){
				_FileSystem_DumpFileInfo(fs,i);
				Log_Add("[FILE SYSTEM] File consistency check failed: File start tag missing or mismatched.");
				Log_Add("[FILE SYSTEM] Deleting File record.");
				memset((void*)&fs->fileSystemHeaderData->fileInfo[i],0,sizeof(FileInfo));
				continue;
			}

			// check for the stop tag at the specified location
			// note that the stop tag is written at the time of creation, so it should always be there...
			// printf("reading pattern at %lu",fs->fileSystemHeaderData->fileInfo[i].stopPosition - fs->fileSystemHeaderData->filesystemInfo.chunkSize);
			s=_FileSystem_SynchronousRead(fs,patternBuff,fs->fileSystemHeaderData->fileInfo[i].stopPosition - fs->fileSystemHeaderData->filesystemInfo.chunkSize, fs->fileSystemHeaderData->filesystemInfo.chunkSize, &bytesRead);
			if (s!=SUCCESS || bytesRead!=fs->fileSystemHeaderData->filesystemInfo.chunkSize){
				Log_Add("[FILE SYSTEM] File consistency check failed: Can't read stop tag.");
				continue;
			}
			if (!_FileSystem_CheckPattern(patternBuff, fs->fileSystemHeaderData->filesystemInfo.chunkSize)){
				_FileSystem_DumpFileInfo(fs,i);
				Log_Add("[FILE SYSTEM] File consistency check failed: File stop tag missing or mismatched.");
				Log_Add("[FILE SYSTEM] Deleting File record.");
				memset((void*)&fs->fileSystemHeaderData->fileInfo[i],0,sizeof(FileInfo));
				continue;
			}

			// last check is to test regions for overlap against all previously checked files
			size_t j=0;
			for (j=1;j<i;j++){
				if (_FileSystem_RegionOverlaps(fs->fileSystemHeaderData->fileInfo[j].startPosition, fs->fileSystemHeaderData->fileInfo[j].stopPosition-1, fs->fileSystemHeaderData->fileInfo[i].startPosition, fs->fileSystemHeaderData->fileInfo[i].stopPosition-1 )){
					_FileSystem_DumpFileInfo(fs,i);
					Log_Add("[FILE SYSTEM] File consistency check failed: file boundaries overlap.");
					Log_Add("[FILE SYSTEM] Conflicting File:");
					_FileSystem_DumpFileInfo(fs,j);
					Log_Add("[FILE SYSTEM] Deleting second file record.");
					memset((void*)&fs->fileSystemHeaderData->fileInfo[i],0,sizeof(FileInfo));
					continue;
				}
			}
			// check to see if the file was closed properly (only applies to write-mode)
			if (fs->fileSystemHeaderData->fileInfo[i].isOpen){
				_FileSystem_DumpFileInfo(fs,i);
				Log_Add("[FILE SYSTEM] File consistency check failed: File was not closed properly. Perhaps the operation was aborted?");
				fs->fileSystemHeaderData->fileInfo[i].isOpen=false;
				// no action necessary, and we should keep going. this is just to notify user ...
			}
			if (!SilentVerify) Log_Add("[FILE SYSTEM] Found file: '%s'",fs->fileSystemHeaderData->fileInfo[i].name);
		}
	}
	if (foundcount!=fs->fileSystemHeaderData->filesystemInfo.number_files_present){
		Log_Add("[FILE SYSTEM] Filesystem consistency check failed: number of files reported does not match the number of files found.");
		fs->fileSystemHeaderData->filesystemInfo.number_files_present=foundcount;
	}
	munmap(patternBuff,fs->fileSystemHeaderData->filesystemInfo.chunkSize);
	Log_Add("[FILE SYSTEM] File check complete.");
	return SUCCESS;
}



StatusCode FileSystem_Open(FileSystem* fs, Disk * disk){
	assert(fs!=NULL);
	assert(disk!=NULL);
	//assert(disk->type==DRIVETYPE_INTERNAL_MULTIDISK_DEVICE);
	//assert(fs->rawDeviceHandle==-1);
	//printf("Opening filesystem.\n");
	/*bzero( (char *)&fs->writeOpInfo, sizeof(struct aiocb) );
	fs->writeOpInfo.aio_fildes=-1;
	bzero( (char *)&fs->readOpInfo, sizeof(struct aiocb) );
	fs->readOpInfo.aio_fildes=-1;
	*/
	fs->rawDeviceHandle = open(disk->deviceName, O_RDWR | O_DIRECT );
	if (fs->rawDeviceHandle == -1) {
		Log_Add("[FILE SYSTEM] Error: could not open filesystem");
		perror("[FILE SYSTEM] open(...)");
		exit(EXIT_CODE_FAIL);
	}else{
		fs->diskInfo=disk;
		// now read header info to determine fs size
		void * buf = mmap(
					NULL,
					sizeof(FileSystemInfo),
					PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_POPULATE | MAP_LOCKED,
					-1,
					0);
		if ((buf!=NULL)&&(buf!=MAP_FAILED)){
			ssize_t bytesRead=0;
			if ((_FileSystem_SynchronousRead(fs, buf,0,sizeof(FileSystemInfo),&bytesRead)==FAILURE) ||
					(bytesRead!=sizeof(FileSystemInfo))){
				//printf("[FILE SYSTEM] bread=%lu\n",bytesRead);
				Log_Add("[FILE SYSTEM] Error: could not read file system info header");

			}
			fs->IMFS_Size = sizeof(FileSystemBlock);
			//printf("----- imfs size = %lu\n",fs->IMFS_Size);
			Log_Add("[FILE SYSTEM] Detected a filesystem size of %lu.",fs->IMFS_Size);
			// we've read the header, so free the resource and prep to load the whole shabang
			munmap(buf,sizeof(FileSystemInfo));
			if (fs->IMFS_Size > (512l*1024l*1024l)){
				Log_Add("[FILE SYSTEM] Error: filesystem reports too large a size, max is 512 MiB, or about 128,000 files");
				exit(EXIT_CODE_FAIL);
			}
			if (fs->IMFS_Size < (32l*4096l)){
				Log_Add("[FILE SYSTEM] Error: filesystem reports too small of a size, min is 128 kiB, or 32 files");
				exit(EXIT_CODE_FAIL);
			}

			if (_FileSystem_AllocateIMFS(fs)==FAILURE){
				Log_Add("[FILE SYSTEM] Allocation failure");
				return FAILURE;
			}

			memset(fs->IMFS_ptr,0,fs->IMFS_Size);
			if ((_FileSystem_SynchronousRead(fs, fs->IMFS_ptr,0,fs->IMFS_Size, &bytesRead)==FAILURE) ||
					(bytesRead!=fs->IMFS_Size)){ // read the whole filesystem in....
				Log_Add("[FILE SYSTEM] FileSystemSynchronousRead failure");
				return FAILURE;
			}
			fs->fileSystemHeaderData=(FileSystemBlock *)fs->IMFS_ptr;
			printf( "                   [FILE SYSTEM] %s: "
					"\n                   [FILE SYSTEM] \tVolume name: %s"
					"\n                   [FILE SYSTEM] \t%13lu chunk size"
					"\n                   [FILE SYSTEM] \t%13lu drives "
					"\n                   [FILE SYSTEM] \t%13lu total bytes "
					"\n                   [FILE SYSTEM] \t%13lu maximum number files"
					"\n                   [FILE SYSTEM] \t%13lu current number files"
					"\n                   [FILE SYSTEM] \t%13lu maximum filename length"
					"\n                   [FILE SYSTEM] \t%13lu bytes reserved for formatting\n",
					fs->diskInfo->deviceName,
					fs->fileSystemHeaderData->filesystemInfo.volname,
					fs->fileSystemHeaderData->filesystemInfo.chunkSize,
					fs->fileSystemHeaderData->filesystemInfo.numDrives,
					fs->fileSystemHeaderData->filesystemInfo.diskSize,
					fs->fileSystemHeaderData->filesystemInfo.maxFiles,
					fs->fileSystemHeaderData->filesystemInfo.number_files_present,
					(unsigned long)FILE_MAX_NAME_LENGTH,
					sizeof(FileSystemBlock)
			);
			if (FileSystem_Verify(fs)==FAILURE){
				Log_Add("[FILE SYSTEM] FileSystemVerify failure");
				return FAILURE;
			}
			ssize_t bytesWritten=0;
			if ((_FileSystem_SynchronousWrite(fs,fs->IMFS_ptr,0,fs->IMFS_Size, &bytesWritten)==FAILURE) ||
					(bytesWritten!=fs->IMFS_Size)){ // write the whole filesystem out...
				Log_Add("[FILE SYSTEM] FileSystemSynchronousWrite failure");
				return FAILURE;
			}
			_FileSystem_RegisterActiveFilesystem(fs);

			fs->closeIsRegistered=true;
			FileSystem_ListFiles(fs);


		} else {
			Log_Add("[FILE SYSTEM] Error: could not allocate space for filesystem header");
			perror("[FILE SYSTEM] mmap(...)");
			exit(EXIT_CODE_FAIL);
		}
	}
	/*
	 *TODO: fixme: no globals!!!
	printf("fs = %lu\n",(size_t) fs);
	printf("globals.fileSystem = %lu \n",(size_t) globals.fileSystem);
	printf("globals.fileSystem->fileSystemHeaderData = %lu \n",(size_t) globals.fileSystem->fileSystemHeaderData);
	printf("globals.fileSystem->fileSystemHeaderData->filesystemInfo.chunkSize = %lu \n",(size_t) globals.fileSystem->fileSystemHeaderData->filesystemInfo.chunkSize);
	 */
	fs->filesystemIsOpen=true;
	fs->filesystemInSync=true;
	return SUCCESS;
}


StatusCode FileSystem_ListFiles(FileSystem * fs){
	assert(fs!=NULL);
	assert(fs->rawDeviceHandle!=-1);
	assert(fs->IMFS_Size!=0);
	assert(fs->IMFS_ptr!=NULL);
	int index=1;
	char buf2[100];
	size_t file_bytes=0;
	size_t f_bytes= sizeof(FileSystemBlock);
	Log_Add("[FILE SYSTEM] Directory listing for volume '%s'",fs->fileSystemHeaderData->filesystemInfo.volname);
	while(index<fs->fileSystemHeaderData->filesystemInfo.number_files_present+1){
		sprintf(buf2,"%.2lu/%.2lu/%.4lu %.2lu:%.2lu:%.2lu",
				fs->fileSystemHeaderData->fileInfo[index].creationDate.day,
				fs->fileSystemHeaderData->fileInfo[index].creationDate.month,
				fs->fileSystemHeaderData->fileInfo[index].creationDate.year,
				fs->fileSystemHeaderData->fileInfo[index].creationDate.hour,
				fs->fileSystemHeaderData->fileInfo[index].creationDate.minute,
				fs->fileSystemHeaderData->fileInfo[index].creationDate.second);
		Log_Add("[FILE SYSTEM] \t%.1s %.20s  %10lu %s",
				((fs->fileSystemHeaderData->fileInfo[index].isOpen) ? "+" : " "),
				buf2,
				fs->fileSystemHeaderData->fileInfo[index].size,
				fs->fileSystemHeaderData->fileInfo[index].name);

		file_bytes+=fs->fileSystemHeaderData->fileInfo[index].size;
		f_bytes+=(fs->fileSystemHeaderData->fileInfo[index].stopPosition - fs->fileSystemHeaderData->fileInfo[index].startPosition)-fs->fileSystemHeaderData->fileInfo[index].size ;
		index++;
	}
	Log_Add("[FILE SYSTEM] Summary:\n                   [FILE SYSTEM]\t%10lu Files\n                   [FILE SYSTEM]\t%18lu Bytes used by files\n                   [FILE SYSTEM]\t%18lu Bytes used for filesystem and formatting",
			fs->fileSystemHeaderData->filesystemInfo.number_files_present,
			file_bytes,
			f_bytes);
	return SUCCESS;
}

size_t FileSystem_GetFreeSpace(FileSystem * fs){
	assert(fs!=NULL);
	assert(fs->rawDeviceHandle!=-1);
	assert(fs->IMFS_Size!=0);
	assert(fs->IMFS_ptr!=NULL);
	int index=1;
	size_t file_bytes=0;
	size_t f_bytes= sizeof(FileSystemBlock);
	while(index<fs->fileSystemHeaderData->filesystemInfo.number_files_present+1){
		file_bytes+=fs->fileSystemHeaderData->fileInfo[index].size;
		f_bytes+=(fs->fileSystemHeaderData->fileInfo[index].stopPosition - fs->fileSystemHeaderData->fileInfo[index].startPosition)-fs->fileSystemHeaderData->fileInfo[index].size ;
		index++;
	}
	return fs->fileSystemHeaderData->filesystemInfo.diskSize - (file_bytes+f_bytes);
}

int _FileSystem_GetFileIndex(FileSystem * fs, const char* name){
	assert(fs->rawDeviceHandle!=-1);
	assert(fs->IMFS_Size!=0);
	assert(fs->IMFS_ptr!=NULL);
	assert(fs->fileSystemHeaderData!=NULL);
	int index=1;
	while(index < fs->fileSystemHeaderData->filesystemInfo.number_files_present+1){
		if(strcmp(name,fs->fileSystemHeaderData->fileInfo[index].name)==0) {
			return  index;
		}
		index++;
	}
	return -1;
}

boolean _FileSystem_IsBetween(size_t a, size_t b, size_t c){
	return (b<=a)&&(a<=c);
}

boolean _FileSystem_RegionOverlaps(size_t starta, size_t stopa, size_t startb, size_t stopb){
	assert(stopa>starta);
	assert(stopb>startb);
	return  _FileSystem_IsBetween(starta, startb, stopb) ||
			_FileSystem_IsBetween(stopa , startb, stopb) ||
			_FileSystem_IsBetween(startb, starta, stopa) ||
			_FileSystem_IsBetween(stopb , starta, stopa);
}

size_t _FileSystem_RangeConflicts(FileSystem * fs, size_t start, size_t stop){
	size_t i=0;
	size_t startf,stopf;
	//printf("Trying range %lu to %lu\n",start,stop);
	for(i=1;i<fs->fileSystemHeaderData->filesystemInfo.number_files_present+1;i++){
		startf=fs->fileSystemHeaderData->fileInfo[i].startPosition;
		stopf=fs->fileSystemHeaderData->fileInfo[i].stopPosition;
		if (startf>0){
			if (_FileSystem_RegionOverlaps(start,stop-1,startf,stopf-1)){
				//printf("Conflicts with file # %d having range %lu to %lu\n",i,startf,stopf);
				return i;
			}
		}
	}
	return 0;
}

void _FileSystem_CreatePattern(void* ptr, size_t size){
	size_t i=0;
	unsigned v=1;
	size_t n=size/(sizeof(unsigned));
	assert (n*(sizeof(unsigned)) == size);
	for(i=0;i<n;i++){
		((unsigned *)ptr)[i] = v;
		v=(v+2)*7;
	}
}

boolean _FileSystem_CheckPattern(void* ptr, size_t size){
	size_t i=0;
	unsigned v=1;
	size_t n=size/(sizeof(unsigned));
	assert (n*(sizeof(unsigned)) == size);
	for(i=0;i<n;i++){
		if(((unsigned *)ptr)[i] != v) return false;
		v=(v+2)*7;
	}
	return true;
}
boolean    FileSystem_IsFileCreatable(FileSystem * fs, char* name, size_t size, boolean overwrite){
	//printf("checking creatability '%s' %lu bytes\n",name,size);
	assert(fs->rawDeviceHandle!=-1);
	int idx = _FileSystem_GetFileIndex(fs, (const char*) name);
	if (idx!=-1){
		if (! overwrite){
			//printf("exists, no overwrite\n",name,size);
			return false;
		} else {
			//printf("exists, overwrite\n",name,size);
			if (fs->fileSystemHeaderData->fileInfo[false].size>size) return true;
		}
	}
	FileInfo fi;
	strncpy(fi.name,name,FILE_MAX_NAME_LENGTH);
	fi.creationDate=getTimestamp();
	fi.currentPosition=0;
	fi.isOpen=0;
	fi.size=size;
	fi.tags=0;

	boolean noconflicts=false;
	size_t i;
	size_t index=fs->fileSystemHeaderData->filesystemInfo.number_files_present+1; // +1 is because the 0th file is actually the filesystem header
	if (index>=fs->fileSystemHeaderData->filesystemInfo.maxFiles){
		//printf("[FILE SYSTEM] too many files\n");
		return false;
	}

	size_t numChunks=((size/fs->fileSystemHeaderData->filesystemInfo.chunkSize)+1);
	size_t sizeOnDisk=(numChunks+2) * fs->fileSystemHeaderData->filesystemInfo.chunkSize;
	if (fs->fileSystemHeaderData->filesystemInfo.number_files_present==0){
		fi.startPosition = fs->IMFS_Size; //start looking at the 1st position after the filesystem block...
		fi.stopPosition  = fi.startPosition + sizeOnDisk;
		//printf("Trying range %lu to %lu\n",fi.startPosition,fi.stopPosition);
		noconflicts=(sizeOnDisk+fs->IMFS_Size < fs->fileSystemHeaderData->filesystemInfo.diskSize);
		//printf("no files noconflicts='%d'\n",noconflicts);
	}else{
		// new code to check space before files
		fi.startPosition = fs->IMFS_Size; //start looking at the 1st position after the filesystem block...
		fi.stopPosition  = fi.startPosition + sizeOnDisk;
		if (!_FileSystem_RangeConflicts(fs,fi.startPosition,fi.stopPosition)){
			noconflicts=true;
		} else {

			for( i=1; i<index; i++){// there are files, so the candidate positions all start at the end of a file :)
				fi.startPosition = fs->fileSystemHeaderData->fileInfo[i].stopPosition; //start looking at the 1st position after this file
				fi.stopPosition  = fi.startPosition + sizeOnDisk;
				if (!_FileSystem_RangeConflicts(fs,fi.startPosition,fi.stopPosition)){
					noconflicts=true;
					break;
				}
			}
		}
		//printf("files noconflicts='%d'\n",noconflicts);
	}
	//return noconflicts;
	if (noconflicts && ((fi.startPosition<fs->fileSystemHeaderData->filesystemInfo.diskSize) &&
			(fi.stopPosition < fs->fileSystemHeaderData->filesystemInfo.diskSize))){
		return true;
	} else {
		return false;
	}
}
StatusCode _WriteStartTag(FileSystem * fs, FileInfo* fi){
	void * patternBuff=mmap(
			NULL,
			fs->fileSystemHeaderData->filesystemInfo.chunkSize,
			PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_POPULATE | MAP_LOCKED,
			-1,
			0);
	if ((patternBuff==NULL)||(patternBuff==MAP_FAILED)){
		Log_Add("[FILE SYSTEM] Error: can't allocate memory for file creation.");
		perror("[FILE SYSTEM] mmap(...)");
		return FAILURE;
	}
	_FileSystem_CreatePattern(patternBuff, fs->fileSystemHeaderData->filesystemInfo.chunkSize);
	//printf("Check  pattern of create: %d",checkPattern(patternBuff, fs.filesystemInfo.chunkSize));
	//printf("writing pattern at %lu",fi.startPosition);
	ssize_t bytesWritten=0;
	if(	_FileSystem_SynchronousWrite(	fs,patternBuff,fi->startPosition,fs->fileSystemHeaderData->filesystemInfo.chunkSize,&bytesWritten) == FAILURE ||
		 bytesWritten != fs->fileSystemHeaderData->filesystemInfo.chunkSize) {
		Log_Add("[FILE SYSTEM] Error: can't write file start tag.");
		perror("[FILE SYSTEM] mmap(...)");
		return FAILURE;
	}
	munmap(patternBuff,fs->fileSystemHeaderData->filesystemInfo.chunkSize);
	return SUCCESS;
}
StatusCode _WriteStopTag(FileSystem * fs, FileInfo* fi){
	void * patternBuff=mmap(
			NULL,
			fs->fileSystemHeaderData->filesystemInfo.chunkSize,
			PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_POPULATE | MAP_LOCKED,
			-1,
			0);
	if ((patternBuff==NULL)||(patternBuff==MAP_FAILED)){
		Log_Add("[FILE SYSTEM] Error: can't allocate memory for file creation.");
		perror("[FILE SYSTEM] mmap(...)");
		return FAILURE;
	}
	_FileSystem_CreatePattern(patternBuff, fs->fileSystemHeaderData->filesystemInfo.chunkSize);
	//printf("Check  pattern of create: %d",checkPattern(patternBuff, fs.filesystemInfo.chunkSize));
	//printf("writing pattern at %lu",fi.startPosition);
	ssize_t bytesWritten=0;
	if(	_FileSystem_SynchronousWrite(fs, patternBuff, fi->stopPosition-fs->fileSystemHeaderData->filesystemInfo.chunkSize,fs->fileSystemHeaderData->filesystemInfo.chunkSize, &bytesWritten)==FAILURE ||
		bytesWritten !=fs->fileSystemHeaderData->filesystemInfo.chunkSize){
		Log_Add("[FILE SYSTEM] Error: can't write file stop tag.");
		perror("[FILE SYSTEM] mmap(...)");
		return FAILURE;
	}
	munmap(patternBuff,fs->fileSystemHeaderData->filesystemInfo.chunkSize);
	return SUCCESS;
}

StatusCode FileSystem_CreateFile(FileSystem * fs, char* name, size_t size, boolean overwrite, int* fileIndex){
	struct timeval dt;

	assert(fs->rawDeviceHandle!=-1);
	assert(fileIndex!=NULL);
	*fileIndex = -1;
	int idx = _FileSystem_GetFileIndex(fs, (const char*) name);
	if (idx!=-1){
		if (! overwrite){
			Log_Add("[FILE SYSTEM] Can't create file, file with name '%s' already exists",name);
			*fileIndex = -1;
			return FAILURE;
		} else {
			FileSystem_DeleteFile(fs,idx);
		}
	}
	Log_Add("[FILE SYSTEM] Creating file '%s' of size %lu bytes",name,size);
	FileInfo fi;
	bzero(&fi,sizeof(FileInfo));
	strncpy(fi.name,name,FILE_MAX_NAME_LENGTH);
	fi.creationDate=getTimestamp();
	fi.currentPosition=0;
	fi.isOpen=0;
	fi.size=size;
	fi.tags=0;
	/* new code for uid,gid, and c,m,a times */
	fi.uid = getuid();
	fi.gid = getgid();
	gettimeofday(&dt, NULL);
	fi.ctime.tv_nsec = 1000l * dt.tv_usec;
	fi.ctime.tv_sec  = dt.tv_sec;
	memcpy((void*)&fi.mtime, (void*)&fi.ctime, sizeof(struct timespec));
	memcpy((void*)&fi.atime, (void*)&fi.ctime, sizeof(struct timespec));
	fi.amode = S_IFREG | S_IWGRP | S_IWOTH | S_IWUSR | S_IRGRP | S_IROTH | S_IRUSR;
	boolean noconflicts=false;
	size_t i=0;
	size_t index=fs->fileSystemHeaderData->filesystemInfo.number_files_present+1; // +1 is because the 0th file is actually the filesystem header
	if (index>=fs->fileSystemHeaderData->filesystemInfo.maxFiles){
		Log_Add("[FILE SYSTEM] Error could not create file, insufficient record space.");
		*fileIndex = -1;
		return FAILURE;
	}

	size_t numChunks=((size/fs->fileSystemHeaderData->filesystemInfo.chunkSize)+1);
	size_t sizeOnDisk=(numChunks+2) * fs->fileSystemHeaderData->filesystemInfo.chunkSize;
	if (fs->fileSystemHeaderData->filesystemInfo.number_files_present==0){
		fi.startPosition = fs->IMFS_Size; //start looking at the 1st position after the filesystem block...
		fi.stopPosition  = fi.startPosition + sizeOnDisk;
		//printf("Trying range %lu to %lu",fi.startPosition,fi.stopPosition);
		noconflicts=(sizeOnDisk+fs->IMFS_Size < fs->fileSystemHeaderData->filesystemInfo.diskSize);
	}else{
		// new code to check space before files
		fi.startPosition = fs->IMFS_Size; //start looking at the 1st position after the filesystem block...
		fi.stopPosition  = fi.startPosition + sizeOnDisk;
		if (!_FileSystem_RangeConflicts(fs,fi.startPosition,fi.stopPosition)){
			noconflicts=true;
		} else {
			for( i=1; i<index; i++){// there are files, so the candidate positions all start at the end of a file :)
				fi.startPosition = fs->fileSystemHeaderData->fileInfo[i].stopPosition; //start looking at the 1st position after this file
				fi.stopPosition  = fi.startPosition + sizeOnDisk;
				if (!_FileSystem_RangeConflicts(fs,fi.startPosition,fi.stopPosition)){
					noconflicts=true;
					break;
				}
			}
		}
	}
	if (noconflicts && ((fi.startPosition<fs->fileSystemHeaderData->filesystemInfo.diskSize) &&
			(fi.stopPosition<fs->fileSystemHeaderData->filesystemInfo.diskSize))){
		fi.currentPosition=fi.startPosition;
		//need to add start and stop tags to make it official
		if (_WriteStartTag(fs,&fi) == FAILURE) {
			*fileIndex = -1;
			return FAILURE;
		}
		if (_WriteStopTag(fs,&fi) == FAILURE) {
			*fileIndex = -1;
			return FAILURE;
		}

		//file is now created, so update filesystem records
		ssize_t bytesWritten=0;
		fs->fileSystemHeaderData->fileInfo[index]=fi;
		fs->fileSystemHeaderData->filesystemInfo.number_files_present++;
		if (_FileSystem_SynchronousWrite(fs,fs->IMFS_ptr,0,fs->IMFS_Size, &bytesWritten)==FAILURE ||
			bytesWritten !=fs->IMFS_Size){
			Log_Add("[FILE SYSTEM] Error: can't synch file system after create file.");
			*fileIndex = -1;
			return FAILURE;
		}
		Log_Add("[FILE SYSTEM] File created successfully");
		*fileIndex = index;
		return SUCCESS;
	} else {
		Log_Add("[FILE SYSTEM] Error could not create file, insufficient contiguous drive space.");
		*fileIndex = -1;
		return FAILURE;
	}
}

StatusCode FileSystem_OpenFile(FileSystem * fs, MyFile* fileptr,const char* fileName, FILEMODE mode){
	assert(fileptr!=NULL);
	assert(fileptr->index==-1);
	assert(fs!=NULL);
	assert(fs->rawDeviceHandle!=-1);
	assert(fs->IMFS_Size!=0);
	assert(fs->IMFS_ptr!=NULL);
	int idx = _FileSystem_GetFileIndex(fs, (const char*) fileName);
	if (idx==-1){
		Log_Add("[FILE SYSTEM] File not found: '%s'",fileName);
		fileptr->index = -1;
		return FAILURE;
	} else {
		//printf("INDEX############################################### %d\n",idx);
		if (fs->fileSystemHeaderData->fileInfo[idx].isOpen){
			Log_Add("[FILE SYSTEM] File already open: '%s'",fileName);
			return FAILURE;
		}
		fileptr->index=-1;
		//printf(HLR("INDEX = %d\n"),idx);
		memset((void *)fileptr,0,sizeof(MyFile));
		//size_t offset= (size_t) fileptr;
		//size_t length= sizeof(MyFile);
		//printf(HLR("CLEARING AT %16lx for %lu bytes\n"),offset,length);
		fileptr->writeOpInfo.aio_fildes  =-1;
		fileptr->readOpInfo.aio_fildes   =-1;
		fileptr->index 				 	 = idx;
		fileptr->fs						 = fs;
		fileptr->lastRPacket 			 = false;
		fileptr->lastRPacketRealSize	 = 0l;
		fileptr->lastWPacket 			 = false;
		fileptr->lastWPacketRealSize	 = 0l;
		fileptr->isOpen					 = true;
		fileptr->mode					 = mode;
		fileptr->reof					 = false;
		fileptr->buffer					 = NULL;
		fileptr->bufferSize				 = 0;
		fileptr->buffer1stByte			 = 0;
		fileptr->contentSize			 = 0;
		fileptr->consumed				 = 0;


		fs->fileSystemHeaderData->fileInfo[fileptr->index].currentPosition = 0;
		fs->fileSystemHeaderData->fileInfo[fileptr->index].isOpen=1;
		fs->fileSystemHeaderData->fileInfo[fileptr->index].mode=mode;
		ssize_t bytesWritten=0;
		StatusCode s = _FileSystem_SynchronousWrite(fs, fs->IMFS_ptr,0,fs->IMFS_Size, &bytesWritten);
		if (bytesWritten!=fs->IMFS_Size || s!=SUCCESS){
			Log_Add("[FILE SYSTEM] Error: can't synch file system after open file.");
			return FAILURE;
		} else {
			return _WriteStartTag(fs, &(fs->fileSystemHeaderData->fileInfo[fileptr->index]));
			/*
			if (FileSystem_Seek(fileptr, 0))
				return FAILURE;
			else
				return SUCCESS;
			*/
		}

	}
}

StatusCode FileSystem_CloseFile(MyFile* file){
	struct timeval dt;
	assert(file!=NULL);
	assert(file->isOpen);
	assert(file->index!=-1);
	assert(file->fs!=NULL);
	assert(file->fs->rawDeviceHandle!=-1);
	assert(file->fs->IMFS_Size!=0);
	assert(file->fs->IMFS_ptr!=NULL);
	assert(file->fs->fileSystemHeaderData->fileInfo[file->index].startPosition!=0);
	assert(file->fs->fileSystemHeaderData->fileInfo[file->index].isOpen!=0);
	gettimeofday(&dt, NULL);
	file->fs->fileSystemHeaderData->fileInfo[file->index].atime.tv_nsec = 1000l * dt.tv_usec;
	file->fs->fileSystemHeaderData->fileInfo[file->index].atime.tv_sec  = dt.tv_sec;
	if (file->mode==WRITE){
		file->fs->fileSystemHeaderData->fileInfo[file->index].size = file->fs->fileSystemHeaderData->fileInfo[file->index].currentPosition;
		memcpy((void*)&file->fs->fileSystemHeaderData->fileInfo[file->index].mtime, (void*)&file->fs->fileSystemHeaderData->fileInfo[file->index].atime, sizeof(struct timespec));
	}
	file->lastRPacket=false;
	file->lastRPacketRealSize=0l;
	file->lastWPacket=false;
	file->lastWPacketRealSize=0l;
	file->mode=READ;
	file->isOpen=false;
	file->reof=true;

	file->fs->fileSystemHeaderData->fileInfo[file->index].currentPosition=0;
	file->fs->fileSystemHeaderData->fileInfo[file->index].isOpen=0;
	file->index=-1;
	StatusCode s;
	ssize_t bytesWritten=0;
	//printf("FileSystem_CloseFile(%s)",file->fs->fileSystemHeaderData->fileInfo[file->index].name);
	s=_FileSystem_SynchronousWrite(file->fs, file->fs->IMFS_ptr,0,file->fs->IMFS_Size, &bytesWritten);
	if (bytesWritten!=file->fs->IMFS_Size){
		Log_Add("[FILE SYSTEM] Error: can't synch file system after close file.");
		return FAILURE;
	} else {
		return s;
	}

}

StatusCode FileSystem_DeleteFile(FileSystem * fs, int fileIndex){
	assert(fs->rawDeviceHandle!=-1);
	assert(fs->IMFS_Size!=0);
	assert(fs->IMFS_ptr!=NULL);
	if (fileIndex<1 || (fileIndex > fs->fileSystemHeaderData->filesystemInfo.number_files_present)){
		Log_Add("[FILE SYSTEM] File not deleted: invalid index specified. '%d'",fileIndex);
		return FAILURE;
	}
	assert(fs->fileSystemHeaderData->fileInfo[fileIndex].startPosition!=0);// why???
	int i=fileIndex;
	while ( i < fs->fileSystemHeaderData->filesystemInfo.number_files_present){
		fs->fileSystemHeaderData->fileInfo[i] = fs->fileSystemHeaderData->fileInfo[i+1];
		i++;
	}
	memset((void*)&fs->fileSystemHeaderData->fileInfo[fs->fileSystemHeaderData->filesystemInfo.number_files_present],0,sizeof(FileInfo));
	fs->fileSystemHeaderData->filesystemInfo.number_files_present--;
	ssize_t bytesWritten=0;
	StatusCode s = _FileSystem_SynchronousWrite(fs, fs->IMFS_ptr, 0 ,fs->IMFS_Size, &bytesWritten);
	if ((bytesWritten != fs->IMFS_Size) || (s != SUCCESS)){
		Log_Add("[FILE SYSTEM] Error: could not synch IMFS after deleting file");
		return FAILURE;
	}
	return SUCCESS;
}

size_t FileSystem_Seek(MyFile* file, size_t offset){
	assert(isAligned(offset));
	assert(file!=NULL);
	assert(file->isOpen);
	assert(file->index!=-1);
	assert(file->fs!=NULL);
	assert(file->fs->rawDeviceHandle!=-1);
	assert(file->fs->IMFS_Size!=0);
	assert(file->fs->IMFS_ptr!=NULL);
	assert(file->fs->fileSystemHeaderData->fileInfo[file->index].isOpen!=0);
	if (offset >= file->fs->fileSystemHeaderData->fileInfo[file->index].size)
		return file->fs->fileSystemHeaderData->fileInfo[file->index].currentPosition;
	file->fs->fileSystemHeaderData->fileInfo[file->index].currentPosition = offset;
	file->lastRPacket=false;
	return offset;
}
StatusCode FileSystem_AsynchronousReadStart(MyFile* file, char * dstPtr, size_t bytes){
	assert(isAligned(dstPtr));
	assert(file!=NULL);
	assert(file->isOpen);
	assert(file->mode==READ);
	assert(file->index!=-1);
	assert(file->fs!=NULL);
	assert(file->fs->rawDeviceHandle!=-1);
	assert(file->fs->IMFS_Size!=0);
	assert(file->fs->IMFS_ptr!=NULL);
	assert(file->fs->fileSystemHeaderData->fileInfo[file->index].isOpen!=0);
	if (file->readOpInfo.aio_fildes==-1) {
		file->readOpInfo.aio_fildes=file->fs->rawDeviceHandle;
		file->readOpInfo.aio_buf=(void*) dstPtr;
		if ((file->fs->fileSystemHeaderData->fileInfo[file->index].currentPosition+bytes)>file->fs->fileSystemHeaderData->fileInfo[file->index].size){// case for the last packet to be read
			if(!file->lastRPacket){
				file->lastRPacket=true;
				file->lastRPacketRealSize=file->fs->fileSystemHeaderData->fileInfo[file->index].size-file->fs->fileSystemHeaderData->fileInfo[file->index].currentPosition;
				file->readOpInfo.aio_nbytes=file->fs->fileSystemHeaderData->filesystemInfo.chunkSize;
			} else {
				Log_Add("[FILE SYSTEM] Error: tried to read a last packet more than once.");
				Log_Add("[FILE SYSTEM] Cpos %15lu",file->fs->fileSystemHeaderData->fileInfo[file->index].currentPosition);
				Log_Add("[FILE SYSTEM] Size %15lu",file->fs->fileSystemHeaderData->fileInfo[file->index].size);
				Log_Add("[FILE SYSTEM] breq %15lu",bytes);
				exit(EXIT_CODE_FAIL);
			}
		} else{
			file->readOpInfo.aio_nbytes=bytes;
		}
		file->readOpInfo.aio_offset =
			file->fs->fileSystemHeaderData->fileInfo[file->index].currentPosition +
			file->fs->fileSystemHeaderData->filesystemInfo.chunkSize +
			file->fs->fileSystemHeaderData->fileInfo[file->index].startPosition;
		if (file->readOpInfo.aio_offset+bytes > (file->fs->fileSystemHeaderData->fileInfo[file->index].stopPosition - file->fs->fileSystemHeaderData->filesystemInfo.chunkSize)){
			Log_Add("[FILE SYSTEM] Error: Attempt to read past the end of the file.");
			bzero( (char *)&file->readOpInfo, sizeof(struct aiocb) );
			file->readOpInfo.aio_fildes=-1;
			return FAILURE;
		} else {
			int status=aio_read(&file->readOpInfo);
			if (status<0) {
				bzero( (char *)&file->readOpInfo, sizeof(struct aiocb) );
				file->readOpInfo.aio_fildes=-1;
				Log_Add("[FILE SYSTEM] Error: Couldn't start asynchronous read.");
				perror("aio_read(...)");
				return NOT_READY;
			} else {
				return SUCCESS;
			}
		}
	} else {
		Log_Add("[FILE SYSTEM] Error: Tried to start a read when a read was already in progress.");
		return FAILURE;
	}
}

StatusCode FileSystem_IsAsynchronousReadDone(MyFile* file){
	assert(file!=NULL);
	assert(file->isOpen);
	assert(file->mode==READ);
	assert(file->index!=-1);
	assert(file->fs!=NULL);
	assert(file->fs->rawDeviceHandle!=-1);
	assert(file->fs->IMFS_Size!=0);
	assert(file->fs->IMFS_ptr!=NULL);
	assert(file->fs->fileSystemHeaderData->fileInfo[file->index].isOpen!=0);
	if (file->readOpInfo.aio_fildes==-1) return FAILURE;
	int status=aio_error(&file->readOpInfo);
	if (status==EINPROGRESS) return NOT_READY;
	return SUCCESS;
}

StatusCode FileSystem_GetAsynchronousBytesRead(MyFile* file, size_t *result){
	assert(file!=NULL);
	assert(file->isOpen);
	assert(file->mode==READ);
	assert(file->index!=-1);
	assert(file->fs!=NULL);
	assert(file->fs->rawDeviceHandle!=-1);
	assert(file->fs->IMFS_Size!=0);
	assert(file->fs->IMFS_ptr!=NULL);
	assert(file->fs->fileSystemHeaderData->fileInfo[file->index].isOpen!=0);
	if (file->readOpInfo.aio_fildes==-1){
		*result=0;
		return FAILURE;
	}

	ssize_t aioResult=aio_return(&file->readOpInfo);
	// check for errors
	if (aioResult<0){
		Log_Add("[FILE SYSTEM] Error: Asynchronous read reported an error.");
		perror("[FILE SYSTEM] aio_return(...)");
		bzero( (char *)&file->readOpInfo, sizeof(struct aiocb) );
		file->readOpInfo.aio_fildes=-1;
		*result=0;
		return FAILURE;
	} else {
		if (!file->lastRPacket){
			if ((size_t)aioResult != file->readOpInfo.aio_nbytes){
				Log_Add("[FILE SYSTEM] Error: Asynchronous read failed to read all bytes.");
			}
			bzero( (char *)&file->readOpInfo, sizeof(struct aiocb) );
			file->readOpInfo.aio_fildes=-1;
			file->fs->fileSystemHeaderData->fileInfo[file->index].currentPosition += (size_t)aioResult;
			*result =  (size_t) aioResult;
			return SUCCESS;
		} else {
			bzero( (char *)&file->readOpInfo, sizeof(struct aiocb) );
			file->readOpInfo.aio_fildes=-1;
			file->fs->fileSystemHeaderData->fileInfo[file->index].currentPosition += file->lastRPacketRealSize;
			file->reof=true;
			*result = file->lastRPacketRealSize;
			return SUCCESS;
		}
	}
}

StatusCode FileSystem_AsynchronousWriteStart(MyFile* file, char * srcPtr, size_t bytes){
	assert(isAligned(srcPtr));
	assert(file!=NULL);
	assert(file->isOpen);
	assert(file->mode==WRITE);
	assert(file->index!=-1);
	assert(file->fs!=NULL);
	assert(file->fs->rawDeviceHandle!=-1);
	assert(file->fs->IMFS_Size!=0);
	assert(file->fs->IMFS_ptr!=NULL);
	assert(file->fs->fileSystemHeaderData->fileInfo[file->index].isOpen!=0);

	if (file->writeOpInfo.aio_fildes==-1) {
		file->writeOpInfo.aio_fildes=file->fs->rawDeviceHandle;
		file->writeOpInfo.aio_buf=(void*) srcPtr;
		if (bytes<file->fs->fileSystemHeaderData->filesystemInfo.chunkSize){// case for the last packet written
			if(!file->lastWPacket){
				file->lastWPacket				= true;
				file->lastWPacketRealSize		= bytes;
				file->writeOpInfo.aio_nbytes	= file->fs->fileSystemHeaderData->filesystemInfo.chunkSize;
			} else {
				Log_Add("[FILE SYSTEM] Error: tried to write a last packet more than once.");
				exit(EXIT_CODE_FAIL);
			}
		} else{
			file->writeOpInfo.aio_nbytes=bytes;
		}
		file->writeOpInfo.aio_offset=
			file->fs->fileSystemHeaderData->fileInfo[file->index].currentPosition +
			file->fs->fileSystemHeaderData->filesystemInfo.chunkSize +
			file->fs->fileSystemHeaderData->fileInfo[file->index].startPosition;
		//printf("SourcePtr %lx\n",(size_t)srcPtr);
		//printf("bytes %lx\n",(size_t)bytes);
		//printf("file->writeOpInfo.aio_offset %lx\n",(size_t)file->writeOpInfo.aio_offset);
		//printf("file->writeOpInfo.aio_nbytes %lx\n",(size_t)file->writeOpInfo.aio_nbytes);
		if (file->writeOpInfo.aio_offset + bytes > (file->fs->fileSystemHeaderData->fileInfo[file->index].stopPosition - file->fs->fileSystemHeaderData->filesystemInfo.chunkSize)){
			Log_Add("[FILE SYSTEM] Error: Attempt to write past the end of the file.");
			Log_Add("\tFile Start: %lu\n\tFile Stop: %lu\n\tCurrent Position: %lu\n\tBytes Tried to Write: %lu\n",
					file->fs->fileSystemHeaderData->fileInfo[file->index].startPosition,
					file->fs->fileSystemHeaderData->fileInfo[file->index].stopPosition,
					file->fs->fileSystemHeaderData->fileInfo[file->index].currentPosition,
					bytes
					);
			//exit(EXIT_CODE_FAIL);
			bzero( (char *)&(file->writeOpInfo), sizeof(struct aiocb) );
			file->writeOpInfo.aio_fildes=-1;
			return FAILURE;
		} else {
			int status=aio_write(&(file->writeOpInfo));
			if (status<0) {
				bzero( (char *)&(file->writeOpInfo), sizeof(struct aiocb) );
				file->writeOpInfo.aio_fildes=-1;
				Log_Add("[FILE SYSTEM] Error: Couldn't start asynchronous write.");
				perror("aio_write(...)");
				return FAILURE;
			} else {
				return SUCCESS;
			}
		}
	} else {
		Log_Add("[FILE SYSTEM] Error: Tried to start a write when a write was already in progress.");
		return FAILURE;
	}
}



StatusCode FileSystem_IsAsynchronousWriteDone(MyFile* file){
	assert(file!=NULL);
	assert(file->isOpen);
	assert(file->mode==WRITE);
	assert(file->index!=-1);
	assert(file->fs!=NULL);
	assert(file->fs->rawDeviceHandle!=-1);
	assert(file->fs->IMFS_Size!=0);
	assert(file->fs->IMFS_ptr!=NULL);
	assert(file->fs->fileSystemHeaderData->fileInfo[file->index].isOpen!=0);
	if (file->writeOpInfo.aio_fildes==-1) return FAILURE;
	int status=aio_error(&(file->writeOpInfo));
	//printf("aio_error %s",strerror(status));
	if (status==EINPROGRESS) return NOT_READY;
	//else printf("aio_error %s",strerror(status));
	return SUCCESS;
}

StatusCode FileSystem_GetAsynchronousBytesWritten(MyFile* file, size_t* result){
	assert(file!=NULL);
	assert(file->isOpen);
	assert(file->mode==WRITE);
	assert(file->index!=-1);
	assert(file->fs!=NULL);
	assert(file->fs->rawDeviceHandle!=-1);
	assert(file->fs->IMFS_Size!=0);
	assert(file->fs->IMFS_ptr!=NULL);
	assert(file->fs->fileSystemHeaderData->fileInfo[file->index].isOpen!=0);
	if (file->writeOpInfo.aio_fildes==-1) {
		*result =  0;
		return FAILURE;
	}
	ssize_t aioResult=aio_return(&(file->writeOpInfo));
	// check for errors
	if (aioResult<0){
		Log_Add("[FILE SYSTEM] Error: Asynchronous write reported an error. %ld", aioResult);
		perror("[FILE SYSTEM] aio_return(...)");
		bzero( (char *)&file->writeOpInfo, sizeof(struct aiocb) );
		file->writeOpInfo.aio_fildes=-1;
		*result =  0;
		return FAILURE;
	} else {
		if (!file->lastWPacket){
			if ((size_t)aioResult != file->writeOpInfo.aio_nbytes){
				Log_Add("[FILE SYSTEM] Error: Asynchronous write failed to write all bytes.");
			}
			bzero( (char *)&file->writeOpInfo, sizeof(struct aiocb) );
			file->writeOpInfo.aio_fildes=-1;
			file->fs->fileSystemHeaderData->fileInfo[file->index].currentPosition+=(size_t)aioResult;
			*result = (size_t)aioResult;
			return SUCCESS;
		} else {
			bzero( (char *)&file->writeOpInfo, sizeof(struct aiocb) );
			file->writeOpInfo.aio_fildes=-1;
			file->fs->fileSystemHeaderData->fileInfo[file->index].currentPosition+=file->lastWPacketRealSize;
			*result = file->lastWPacketRealSize;
			return SUCCESS;
		}
	}
}
boolean    FileSystem_FileExists(FileSystem * fs,const char* fileName){
	return (_FileSystem_GetFileIndex(fs,fileName)!=-1);
}
size_t FileSystem_GetFileCount(FileSystem * fs){
	assert (fs!=NULL);
	assert (fs->fileSystemHeaderData!=NULL);
	return fs->fileSystemHeaderData->filesystemInfo.number_files_present;
}
void FileSystem_GetFileMetaData(FileSystem * fs, int index, FileMetaData** data){
	if (!fs || index<0 || index>fs->fileSystemHeaderData->filesystemInfo.number_files_present || !data){
		Log_Add("[FILE SYSTEM] Error: tried to get metadata, but no filesystem or bad index or bad result pointer.");
		return;
	}
	size_t size=fs->fileSystemHeaderData->fileInfo[index].size;
	size_t AllocatedSize =
		fs->fileSystemHeaderData->fileInfo[index].stopPosition -
		fs->fileSystemHeaderData->fileInfo[index].startPosition;

	sprintf(fs->fileSystemHeaderData->fileInfo[index].metaData.size,"%015lu",size);
	sprintf(fs->fileSystemHeaderData->fileInfo[index].metaData.usage,"%015lu",AllocatedSize);
	*data = &(fs->fileSystemHeaderData->fileInfo[index].metaData);
	Log_Add("[FILE SYSTEM] Get File Metadata for '%s'.",fs->fileSystemHeaderData->fileInfo[index].name);
}
StatusCode FileSystem_GetFileFormat(FileSystem * fs, int index,char * buffer){
	FileMetaData* fmd = NULL;
	FileSystem_GetFileMetaData(fs, index, &fmd);
	if (fmd==NULL)
		return FAILURE;
	else {
		strncpy(buffer,fmd->format,33);
		buffer[32]='\0';
		return SUCCESS;
	}
}
void FileSystem_SetFileMetaData(FileSystem * fs, int index, size_t StartMPM, size_t StartMJD, size_t reference, size_t StopMPM, size_t StopMJD, char * format){
	if (!fs || index<0 || index>fs->fileSystemHeaderData->filesystemInfo.number_files_present){
		Log_Add("[FILE SYSTEM] Error: specified metadata, but no filesystem or bad index.");
		return;
	}
	sprintf(fs->fileSystemHeaderData->fileInfo[index].metaData.tag,"%06lu_%09lu",StartMJD,reference);
	sprintf(fs->fileSystemHeaderData->fileInfo[index].metaData.StartMJD,"%06lu",StartMJD);
	sprintf(fs->fileSystemHeaderData->fileInfo[index].metaData.StartMPM,"%09lu",StartMPM);
	sprintf(fs->fileSystemHeaderData->fileInfo[index].metaData.StopMJD,"%06lu",StopMJD);
	sprintf(fs->fileSystemHeaderData->fileInfo[index].metaData.StopMPM,"%09lu",StopMPM);
	sprintf(fs->fileSystemHeaderData->fileInfo[index].metaData.format,"%32s",format);
	sprintf(fs->fileSystemHeaderData->fileInfo[index].metaData.size,"000000000000000");
	sprintf(fs->fileSystemHeaderData->fileInfo[index].metaData.usage,"000000000000000");
	sprintf(fs->fileSystemHeaderData->fileInfo[index].metaData.complete,"NO ");
	Log_Add("[FILE SYSTEM] Set File Metadata for '%s'.",fs->fileSystemHeaderData->fileInfo[index].name);
}
void FileSystem_SetFileMetaDataIsComplete(FileSystem * fs, int index, int complete){
	if (!fs || index<0 || index>fs->fileSystemHeaderData->filesystemInfo.number_files_present){
		Log_Add("[FILE SYSTEM] Error: specified metadata, but no filesystem or bad index.");
		return;
	}
	if (complete){
		sprintf(fs->fileSystemHeaderData->fileInfo[index].metaData.complete,"YES");
	} else {
		sprintf(fs->fileSystemHeaderData->fileInfo[index].metaData.complete,"NO ");
	}
	Log_Add("[FILE SYSTEM] Set Complete Flag for '%s'.",fs->fileSystemHeaderData->fileInfo[index].name);
}




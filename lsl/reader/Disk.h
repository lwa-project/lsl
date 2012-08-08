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
 * Disk.h
 *
 *  Created on: Oct 18, 2009
 *      Author: chwolfe2
 */

#ifndef DISK_H_
#define DISK_H_
#include <stdlib.h>
#include "Defines.h"

typedef int DriveType;
#define DRIVETYPE_UNKNOWN					0
#define DRIVETYPE_SYSTEM_BLOCK_DEVICE		1
#define DRIVETYPE_SYSTEM_PARTITION			2
#define DRIVETYPE_SYSTEM_OPTICAL			3
#define DRIVETYPE_SYSTEM_FLOPPY				4
#define DRIVETYPE_INTERNAL_BLOCK_DEVICE		5
#define DRIVETYPE_INTERNAL_MULTIDISK_DEVICE	6
#define DRIVETYPE_EXTERNAL_BLOCK_DEVICE		7
#define DRIVETYPE_EXTERNAL_PARTITION		8
#define DRIVETYPE_EXTERNAL_MULTIDISK_DEVICE	9

#define MAX_NUM_DEVICES_TO_ENUMERATE 50
#define MAX_DEVICE_NAME_LENGTH				64
#define EXT_FILE_MAX_NAME_LENGTH			128

typedef struct __Disk{
	char  		deviceName[MAX_DEVICE_NAME_LENGTH];
	boolean		locked;
	DriveType	type;
	size_t 		bytesUsed;
	size_t 		bytesTotal;
	boolean		mounted;
	dev_t		devNumber;
	//char		mountpoint[MAX_DEVICE_NAME_LENGTH];
	//size_t		numberOfChildren;
	//char  		parent[MAX_DEVICE_NAME_LENGTH];
} Disk;
#define DiskSize sizeof(Disk)

StatusCode Disk_IdentifyAll(void);
StatusCode Disk_GetDiskInfo(char *  deviceName, Disk* disk);
StatusCode Disk_CreateArray(Disk disksToUse[], int count);
StatusCode Disk_Format(char* deviceName);
StatusCode Disk_MountByName(char* deviceName);
StatusCode Disk_Mount(Disk* disk);
StatusCode Disk_UnmountByName(char* deviceName);
StatusCode Disk_Unmount(Disk* disk);
StatusCode Disk_UpdateFreeSpace(Disk* disk);

int DiskGetType(int index);
size_t DiskGetFreeSpace(int index);
int DiskFindIndexByName(char * deviceName);
size_t Disk_GetArrayNumDisks(Disk* disk);
size_t Disk_GetArraySize(Disk* disk);
size_t Disk_GetArrayChunkSize(Disk* disk);

extern const char *DriveTypeString[];

int    DiskGetUsableCount(void);
char * DiskGetUsableName(int index);
size_t DiskGetUsableFreeSpace(int index);
size_t DiskGetUsableTotalSpace(int index);

#endif /* DISK_H_ */

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
 * Disk.c
 *
 *  Created on: Oct 18, 2009
 *      Author: chwolfe2
 */


#include "string.h"
#include "Disk.h"
#include "HostInterface.h"
#include "FileSystem.h"
#include "Log.h"
//#include "Globals.h"
#include <stdio.h>
#include <assert.h>
#include <errno.h>
#include <sys/mount.h>
#include <sys/statvfs.h>

// define some constants to aid in device/partition type recognition
//const char* SysBlockDevs[]={"sda"};
//const int NumSysBlockDevs = 1;
//const char* SysPartitions[]={"sda1","sda2","sda5"};
//const int NumSysPartitions = 3;
//const char* IntMultiDiskPartitions[]={"sdb", "sdc", "sdd", "sde", "sdf"};
//const int NumIntMultiDiskPartitions = 5;
//(DEBUGGING))const int NumIntMultiDiskPartitions = 5;
const char* IntMultiDiskDevs[]={"md0", "md1", "md2", "md3", "md4", "md5"};
const int NumIntMultiDiskDevs = 6;
const char * DriveTypeString[]={
	"DRIVETYPE_UNKNOWN",
	"DRIVETYPE_SYSTEM_BLOCK_DEVICE",
	"DRIVETYPE_SYSTEM_PARTITION",
	"DRIVETYPE_SYSTEM_OPTICAL",
	"DRIVETYPE_SYSTEM_FLOPPY",
	"DRIVETYPE_INTERNAL_BLOCK_DEVICE",
	"DRIVETYPE_INTERNAL_MULTIDISK_DEVICE",
	"DRIVETYPE_EXTERNAL_BLOCK_DEVICE",
	"DRIVETYPE_EXTERNAL_PARTITION",
	"DRIVETYPE_EXTERNAL_MULTIDISK_DEVICE"
};



// record to hold information regarding detected devices
Disk 		devices[MAX_NUM_DEVICES_TO_ENUMERATE];
size_t 		devicesCount=0;
boolean 	deviceRecordLocked=0;

int DiskGetType(int index){
	if (index<0 || index>devicesCount) return DRIVETYPE_UNKNOWN;
	return devices[index].type;
}
size_t DiskGetFreeSpace(int index){
	if (index<0 || index>devicesCount) return 0;
	switch (DiskGetType(index)){
		case DRIVETYPE_INTERNAL_MULTIDISK_DEVICE:
			//return FileSystem_GetFreeSpace(globals.fileSystem);
			//TODO: fixme no globals!!!!!
			return 0;
			break;
		case DRIVETYPE_EXTERNAL_PARTITION:
			return devices[index].bytesTotal - devices[index].bytesUsed;
			break;
		default:
			return 0;
	}

}
int DiskFindIndexByName(char * deviceName){
	int i;
	for (i=0;i<devicesCount;i++){
		if (strcmp(deviceName,devices[i].deviceName)==0)	return i;
		//printf("%s is not %s\n",deviceName,devices[i].deviceName);
	}
	return -1;
}


boolean matchStringInStringList(char* string, const char** matchList, int matchListSize ){
	int i;
	for (i=0;i<matchListSize;i++){
		if (strcmp(string,matchList[i])==0)
			return true;
	}
	return false;
}
StatusCode DiskLockDeviceInfo(){
	if (deviceRecordLocked) return FAILURE;
	deviceRecordLocked=true;
	return SUCCESS;
}
StatusCode DiskUnlockDeviceInfo(){
	deviceRecordLocked=false;
	return SUCCESS;
}
StatusCode DiskLockDevice(char * deviceName){
	int i=DiskFindIndexByName(deviceName);
	if (i<0) return FAILURE;
	if (devices[i].locked) return FAILURE;
	devices[i].locked=true;
	return SUCCESS;
}
StatusCode DiskUnlockDevice(char * deviceName){
	int i=DiskFindIndexByName(deviceName);
	if (i<0) return FAILURE;
	devices[i].locked=false;
	return SUCCESS;
}
char * getDevShortName(char*deviceName){
	return &deviceName[5];
}
/**********************************************************************************
new code to prevent CV12 firmware drives from being used as DRSUs (C.N.W. 27MAY11)
**********************************************************************************/
boolean partitionIsOnCV12FirmwareDrive(const char* linuxPartitionDeviceIdentifier){
	char cmd[2048];
	sprintf(cmd,"hdparm -I %s | grep CV12 | wc -l",linuxPartitionDeviceIdentifier);
	size_t linecount;
	ShellGetNumber(cmd,&linecount);
	if (linecount==1)
		return TRUE;
	else
		return FALSE;
}
boolean arrayContainsCV12FirmwareDrives(const char* id){
	char cmd[2048];
	sprintf(cmd,"for x in `cat /proc/mdstat | grep %s | cut -d \" \" -f 5- | sed -e 's/\\[[0-9]\\]//g'`; do hdparm -I /dev/$x | grep CV12; done | wc -l",id);
	size_t linecount;
	ShellGetNumber(cmd,&linecount);
	if (linecount==1)
		return TRUE;
	else
		return FALSE;
}
DriveType detectMultidiskFileSystemType(const char* id){
	char fullDeviceName[1024];	bzero((void*)fullDeviceName, 1024);
	char cmd[2048];				bzero((void*)cmd, 2048);
	char typeName[1024];		bzero((void*)typeName, 1024);
	sprintf(fullDeviceName,"/dev/%s",id);
	sprintf(cmd,"blkid -s TYPE -o value %s",fullDeviceName);
	ShellGetString(cmd,typeName,1024);
	if ((strncmp(typeName, "ext3", 4)==0) || (strncmp(typeName, "ext3", 4)==0)){
		return DRIVETYPE_UNKNOWN;
	} else if (strncmp(typeName, "ext2", 4)==0){
		return DRIVETYPE_EXTERNAL_MULTIDISK_DEVICE;
	} else if (strlen(typeName)==0){
		if (arrayContainsCV12FirmwareDrives(id)){
			return DRIVETYPE_UNKNOWN;
		} else {
			// note that there is no need to actually contain an LWAFS filesystem since we
			// can initialize the drives later
			return DRIVETYPE_INTERNAL_MULTIDISK_DEVICE;
		}
	} else {
		return DRIVETYPE_UNKNOWN;
	}

}
/**********************************************************************************
 end new code
**********************************************************************************/

StatusCode Disk_IdentifyAll(){
	if (deviceRecordLocked) return FAILURE;
	char buffer[2048];
	int i=0;
	int j=0;
	if (issueShellCommand("cat /proc/partitions  | awk '{if (NR>2)print $4}'", buffer, 2048) == -1) {
		return FAILURE;
	}
	if (strlen(buffer)==0){
		return FAILURE;
	}
	// now process the output to determine which devices of what type are present
	devicesCount=0;
	char * token =strtok(buffer," \n");
	char cmdTemp [1024];
	// pass 1 finds instantly recognizable drive device IDs and assigns type/size information
	while (token!=NULL){
		sprintf(devices[devicesCount].deviceName,"/dev/%s",token);
		sprintf(cmdTemp,"sfdisk -s \"/dev/%s\"",token);
		ShellGetNumber(cmdTemp,&devices[devicesCount].bytesTotal);
		devices[devicesCount].bytesTotal *= 1024l;
		sprintf(cmdTemp,"stat \"/dev/%s\" | grep Device | nawk '{print $2}' | sed -e's/.*\\///' | tr -d 'd'",token);
		ShellGetNumber(cmdTemp,&devices[devicesCount].devNumber);



		if((strncmp(token, "sd", 2)==0) || (strncmp(token, "hd", 2)==0) || (strncmp(token, "loop", 4)==0)){
			devices[devicesCount].type = DRIVETYPE_SYSTEM_BLOCK_DEVICE;
			//printf("00FOUND: %s setting to DRIVETYPE_SYSTEM_BLOCK_DEVICE\n",token);
		}
		/********************************************************************************
			New drive detection semantics. previously assumed md0 .. md6 were internal storage.
			unfortunately, managing linux device identifiers is not as straightforward as it should be.
			therefore, we must detect the filesystems (ext2/3/4 vs. LWAFS) directly
		********************************************************************************/
		else if(strncmp(token, "md", 2)==0){
			// detectMultidiskFileSystemType() will return one of :
			//    A) DRIVETYPE_INTERNAL_MULTIDISK_DEVICE
			//    B) DRIVETYPE_EXTERNAL_MULTIDISK_DEVICE
			//    C) DRIVETYPE_UNKNOWN
			devices[devicesCount].type = detectMultidiskFileSystemType(token);
		/********************************************************************************
		 * Old detection semantics follow:
		else if(matchStringInStringList(token, IntMultiDiskDevs, NumIntMultiDiskDevs)) {
			devices[devicesCount].type = DRIVETYPE_INTERNAL_MULTIDISK_DEVICE;
			//printf("01FOUND: %s setting to DRIVETYPE_INTERNAL_MULTIDISK_DEVICE\n",token);
		} else if(strncmp(token, "md", 2)==0){
			//by this point, all of the legal DRSU names have been accounted for, and remaining should be external arrays
			devices[devicesCount].type = DRIVETYPE_EXTERNAL_MULTIDISK_DEVICE;
			//printf("02FOUND: %s setting to DRIVETYPE_EXTERNAL_MULTIDISK_DEVICE\n",token);
		********************************************************************************/


		}

/*
		if (matchStringInStringList(token, SysBlockDevs, NumSysBlockDevs)) {
			devices[devicesCount].type = DRIVETYPE_SYSTEM_BLOCK_DEVICE;
			devices[devicesCount].mounted=false;
			devices[devicesCount].locked=true;
		} else if(matchStringInStringList(token, SysPartitions, NumSysPartitions)) {
			devices[devicesCount].type = DRIVETYPE_SYSTEM_PARTITION;
			devices[devicesCount].mounted=true;
		} else if(matchStringInStringList(token, IntMultiDiskPartitions, NumIntMultiDiskPartitions)) {
			devices[devicesCount].type = DRIVETYPE_INTERNAL_BLOCK_DEVICE;
			devices[devicesCount].mounted=false;
			devices[devicesCount].locked=true;
		} else if(matchStringInStringList(token, IntMultiDiskDevs, NumIntMultiDiskDevs)) {
			devices[devicesCount].type = DRIVETYPE_INTERNAL_MULTIDISK_DEVICE;
		}*/

		else if(strncmp(token, "fd", 2)==0) {
			devices[devicesCount].type = DRIVETYPE_SYSTEM_FLOPPY;
			//printf("03FOUND: %s setting to DRIVETYPE_SYSTEM_FLOPPY\n",token);

		} else if(strncmp(token, "sr", 2)==0) {
			devices[devicesCount].type = DRIVETYPE_SYSTEM_OPTICAL;
			//printf("04FOUND: %s setting to DRIVETYPE_SYSTEM_OPTICAL\n",token);

		} else {
			// at this point, options are DRIVETYPE_UNKNOWN, DRIVETYPE_INTERNAL_BLOCK_DEVICE, DRIVETYPE_EXTERNAL_BLOCK_DEVICE, DRIVETYPE_EXTERNAL_PARTITION
			//DONE cnw, 22AUG10 // TODO: advanced type determination
			devices[devicesCount].type = DRIVETYPE_UNKNOWN;
			//printf("05FOUND: %s setting to DRIVETYPE_UNKNOWN\n",token);


		}
		devicesCount++;
		token = strtok(NULL," \n"); //continue tokenization
	}
	/*
	// TODO: occupied size determination
			sprintf(cmdTemp,"",token);
			ShellGetNumber(cmdTemp,&devices[devicesCount].bytesTotal);
			devices[devicesCount].bytesUsed=999999999999l;*/
	//printf("\nPASS 2\n");
	// pass 2 looks at multi-disk arrays and assigns the proper type to array member devices
	for (i=0;i<devicesCount; i++){
		if (devices[i].type == DRIVETYPE_INTERNAL_MULTIDISK_DEVICE){
			//printf("MULTI_DISK_INTERNAL_CHECKOMATIC %s\n",devices[i].deviceName);
			sprintf(cmdTemp,
					"cat /proc/mdstat | awk '/^%s/,/chunks$/ {for(i=4;i<NF;i++){print $(i+1)}}' | sed 's/\\[[0-9]*\\]//'"
					,getDevShortName(devices[i].deviceName));
			//printf("\n>>>>>>> %s\n",cmdTemp);
			if (issueShellCommand(cmdTemp, buffer, 2048) == -1) {
				return FAILURE;
			} else {
				//printf("BUFFER= '%s'\n",buffer);
				// now buffer contains a list of devices belonging to our multi-disk device
				token =strtok(buffer," \n");
				while (token!=NULL){
					for(j=0;j<devicesCount; j++){
						//printf("1 Test : %s == %s\n",token,getDevShortName(devices[j].deviceName));
						if (strcmp(token,getDevShortName(devices[j].deviceName))==0){
							devices[j].type = DRIVETYPE_INTERNAL_BLOCK_DEVICE;
							devices[devicesCount].mounted=false;
							devices[devicesCount].locked=true;
							//printf("06FOUND: %s setting to DRIVETYPE_INTERNAL_BLOCK_DEVICE\n",devices[j].deviceName);

						}
					}
					token = strtok(NULL," \n"); //continue tokenization
				}
			}
		} else if (devices[i].type == DRIVETYPE_EXTERNAL_MULTIDISK_DEVICE){
			//printf("MULTI_DISK_EXTERNAL_CHECKOMATIC %s\n",devices[i].deviceName);
			sprintf(cmdTemp,
					"cat /proc/mdstat | awk '/^%s/,/chunks$/ {for(i=4;i<NF;i++){print $(i+1)}}' | sed 's/\\[[0-9]*\\]//'"
					,getDevShortName(devices[i].deviceName));
			//printf("\n>>>>>>> %s\n",cmdTemp);
			if (issueShellCommand(cmdTemp, buffer, 2048) == -1) {
				return FAILURE;
			} else {
				//printf("BUFFER= '%s'\n",buffer);
				// now buffer contains a list of devices belonging to our multi-disk device
				token =strtok(buffer," \n");
				while (token!=NULL){
					for(j=0;j<devicesCount; j++){
						//printf("2 Test : %s == %s\n",token,getDevShortName(devices[j].deviceName));
						if (strcmp(token,getDevShortName(devices[j].deviceName))==0){
							devices[j].type = DRIVETYPE_EXTERNAL_BLOCK_DEVICE;
							//printf("07FOUND: %s setting to DRIVETYPE_EXTERNAL_BLOCK_DEVICE\n",devices[j].deviceName);
						}
					}
					token = strtok(NULL," \n"); //continue tokenization
				}
			}

		}
	}

	//for (i=0;i<devicesCount; i++){Log_Add("[DISK MANAGER] Disk_IdentifyAll(): Device found: '%s' of type '%s'", devices[i].deviceName, DriveTypeString[devices[i].type]);}
	//printf("\nPASS 3\n");
	// pass 3 looks for block device parenting and considers childless block devices to be legal partitions
	int childFound=0;
	size_t length1,length2;
	for (i=0;i<devicesCount; i++){

		if (devices[i].type == DRIVETYPE_SYSTEM_BLOCK_DEVICE){
			//printf("Test %s for children\n",devices[i].deviceName);
			childFound=0;
			length1=strlen(devices[i].deviceName);
			for (j=0;j<devicesCount; j++){
				if (j!=i){
					length2=strlen(devices[j].deviceName);
					if (length2>length1){
						if (strncmp(devices[i].deviceName,devices[j].deviceName,length1)==0){
							childFound=1;
						}
					}
				}
			}
			if (!childFound) {
				//could also be an internal partition, but we'll remove those in next pass
				devices[i].type = DRIVETYPE_EXTERNAL_PARTITION;
				//printf("08UPDATE: %s setting to DRIVETYPE_EXTERNAL_PARTITION\n",devices[i].deviceName);

			} else {
				devices[i].type = DRIVETYPE_SYSTEM_BLOCK_DEVICE; //both external and internal;
				// ideally we would differentiate between the two, but since neither is usable, we don't
				// differentiate further.
				devices[devicesCount].mounted=false;
				devices[devicesCount].locked=true;

			}
		}
	}

	//for (i=0;i<devicesCount; i++){Log_Add("[DISK MANAGER] Disk_IdentifyAll(): Device found: '%s' of type '%s'", devices[i].deviceName, DriveTypeString[devices[i].type]);}
	//printf("\nPASS 4\n");
	// pass 4: treat external arrays as external partitions
	for (i=0;i<devicesCount; i++){
		if (devices[i].type == DRIVETYPE_EXTERNAL_MULTIDISK_DEVICE){
			devices[i].type = DRIVETYPE_EXTERNAL_PARTITION;
		}

	}

	//for (i=0;i<devicesCount; i++){Log_Add("[DISK MANAGER] Disk_IdentifyAll(): Device found: '%s' of type '%s'", devices[i].deviceName, DriveTypeString[devices[i].type]);}
	//printf("\nPASS 5a\n");
	// pass 5 removes system partitions from partitions currently identified as DRIVETYPE_EXTERNAL_PARTITION
	// check that
	//  A) partition does not contain root mount point
	//  B) partition contains ext2 file system (ext3/4 to be supported soon)

	// A:
	// get the device which contains the root mountpoint
	if (issueShellCommand("mount | grep 'on / ' | nawk '{print $1}'", buffer, 2048) == -1) {
		return FAILURE;
	}
	//buffer now contains the device name of the root mountpoint
	for (i=0;i<devicesCount; i++){
		if (devices[i].type == DRIVETYPE_EXTERNAL_PARTITION){
			//printf("Test : %s == %s\n",devices[i].deviceName,buffer);
			if (strcmp(devices[i].deviceName,buffer)==0){
				devices[i].type=DRIVETYPE_SYSTEM_PARTITION;
				devices[devicesCount].mounted=true;
				//printf("09UPDATE: %s setting to DRIVETYPE_SYSTEM_PARTITION\n",devices[i].deviceName);

			}
		}
	}

	//for (i=0;i<devicesCount; i++){Log_Add("[DISK MANAGER] Disk_IdentifyAll(): Device found: '%s' of type '%s'", devices[i].deviceName, DriveTypeString[devices[i].type]);}
	//printf("\nPASS 5b\n");
	// B:
	// get a list of partitions with ext2 file systems
	if (issueShellCommand("blkid | grep ext2 | nawk '{print $1}' | tr -d ':' | tr '\\n' ' '", buffer, 2048) == -1) {
		return FAILURE;
	}
	//buffer now contains a list of partitions with ext2 file systems
	for (i=0;i<devicesCount; i++){
		if (devices[i].type == DRIVETYPE_EXTERNAL_PARTITION){
			//printf("Test : %s in '%s'\n",devices[i].deviceName,buffer);
			if (!strstr(buffer,devices[i].deviceName)){
				devices[i].type=DRIVETYPE_SYSTEM_PARTITION;
			}
		}
	}

	//printf("\nPASS 6\n");
	// pass 6: check drive free space
	size_t blocksUsed;
	for (i=0;i<devicesCount; i++){
		blocksUsed=0;
		switch(devices[i].type){
			case DRIVETYPE_UNKNOWN					: //fall through
			case DRIVETYPE_SYSTEM_BLOCK_DEVICE		: //fall through
			case DRIVETYPE_SYSTEM_PARTITION			: //fall through
			case DRIVETYPE_SYSTEM_OPTICAL			: //fall through
			case DRIVETYPE_SYSTEM_FLOPPY			: //fall through
			case DRIVETYPE_INTERNAL_BLOCK_DEVICE	: //fall through
			case DRIVETYPE_EXTERNAL_BLOCK_DEVICE	: //fall through
			case DRIVETYPE_EXTERNAL_MULTIDISK_DEVICE: //fall through (note that external multidisk devices should already be re-classed as external partitions)
			case DRIVETYPE_INTERNAL_MULTIDISK_DEVICE:
			default :
				devices[i].bytesUsed=devices[i].bytesTotal;
				break;
			case DRIVETYPE_EXTERNAL_PARTITION		:
				sprintf(cmdTemp,"df | grep \"%s \" | nawk '{print $3}'",devices[i].deviceName);
				if (ShellGetNumber(cmdTemp, &blocksUsed)!=SUCCESS)
					Log_Add("[Disk] Warning: freespace determination for '%s' may be inaccurate",devices[i].deviceName);
				devices[i].bytesUsed=(blocksUsed * 1024);
				break;
		}
	}

	// pass 7: print discovered drives
	for (i=0; i<devicesCount; i++){
		Log_Add("[DISK MANAGER] Device found: '%s' of type '%s'; Total size: %lu bytes", devices[i].deviceName, DriveTypeString[devices[i].type], devices[i].bytesTotal);
	}
	return (devicesCount) ? SUCCESS : FAILURE;
}
StatusCode Disk_GetDiskInfo(char *  deviceName, Disk* disk){
	int i = DiskFindIndexByName(deviceName);
	if (i!=-1){
		*disk=devices[i];
		return SUCCESS;
	}

	// if we got here, then we didn't find the drive
	return FAILURE;

}
StatusCode Disk_Format(char* deviceName){
	int i = DiskFindIndexByName(deviceName);
	char buffer[2048];
	if (i!=-1){
		if (devices[i].mounted || devices[i].locked) return FAILURE;
		if ((devices[i].type == DRIVETYPE_EXTERNAL_PARTITION) || ((devices[i].type == DRIVETYPE_INTERNAL_MULTIDISK_DEVICE))){
			if (devices[i].type == DRIVETYPE_EXTERNAL_PARTITION){
				if (issueShellCommand("mke2fs %s -q", buffer, 2048) != -1) return SUCCESS;
			} else {
				return FileSystem_Create(&devices[i]);
			}
		}
	}
	// if we got here, then we didn't format the drive, or an error occurred
	return FAILURE;
}

StatusCode Disk_MountByName(char* deviceName){
	int i = DiskFindIndexByName(deviceName);
	char buffer[2048];
	char mountpoint [1024];
	char cmdTemp [1024];
	//printf("mount()\n");
	if (i!=-1){
		if (devices[i].mounted) return SUCCESS;
		if (devices[i].locked) return FAILURE;
		//printf("not mounted locked\n");
		if ((devices[i].type == DRIVETYPE_EXTERNAL_PARTITION) || ((devices[i].type == DRIVETYPE_INTERNAL_MULTIDISK_DEVICE))){
			//printf("good type\n");
			if (devices[i].type == DRIVETYPE_EXTERNAL_PARTITION){
				sprintf(mountpoint,"/LWA_EXT/%s",deviceName);
				sprintf(cmdTemp,"mkdir -p /%s",mountpoint);
				if (issueShellCommand(cmdTemp, buffer, 2048) == -1) {
					printf("[DISK MANAGER] Error: mount failed: Could not create mountpoint.\n");
					return FAILURE;
				}
				if (mount(deviceName,mountpoint,"ext2",0,NULL)==-1){
					switch (errno) {
					case EBUSY:
						//source is already mounted. Or, it cannot be remounted read-only, because it still holds files open for writing. Or, it cannot be mounted on target because target is still busy (it is the working directory of some task, the mount point of another device, has open files, etc.). Or, it could not be unmounted because it is busy.
						if (!umount(mountpoint) && !mount(deviceName,mountpoint,"ext2",0,NULL)){
							// unmount then mount succeeded
							return SUCCESS;
						}else {
							// unmount or mount failed
						}
					case EACCES:
						//A component of a path was not searchable. (See also path_resolution(2).) Or, mounting a read-only filesystem was attempted without giving the MS_RDONLY flag. Or, the block device source is located on a filesystem mounted with the MS_NODEV option.
					case EPERM:
						//The caller does not have the required privileges.
						//Conforming to

					case EAGAIN:
						//A call to umount2() specifying MNT_EXPIRE successfully marked an unbusy file system as expired.
					case EFAULT:
						//One of the pointer arguments points outside the user address space.
					case EINVAL:
						//source had an invalid superblock. Or, a remount (MS_REMOUNT) was attempted, but source was not already mounted on target. Or, a move (MS_MOVE) was attempted, but source was not a mount point, or was '/'. Or, an unmount was attempted, but target was not a mount point. Or, umount2() was called with MNT_EXPIRE and either MNT_DETACH or MNT_FORCE.
					case ELOOP:
						//Too many link encountered during pathname resolution. Or, a move was attempted, while target is a descendant of source.
					case EMFILE:
						//(In case no block device is required:) Table of dummy devices is full.
					case ENAMETOOLONG:
						//A pathname was longer than MAXPATHLEN.
					case ENODEV:
						//filesystemtype not configured in the kernel.
					case ENOENT:
						//A pathname was empty or had a nonexistent component.
					case ENOMEM:
						//The kernel could not allocate a free page to copy filenames or data into.
					case ENOTBLK:
						//source is not a block device (and a device was required).
					case ENOTDIR:
						//The second argument, or a prefix of the first argument, is not a directory.
					case ENXIO:
						//The major number of the block device source is out of range.
					default:
						printf("[DISK MANAGER] Error: mount failed: %s\n",strerror(errno));
						return FAILURE;
					}
				} else {
					//initial mount succeeded
					return SUCCESS;
				}
				/*sprintf(cmdTemp,"mkdir -p /LWA_EXT/%s",deviceName);
				if (issueShellCommand(cmdTemp, buffer, 2048) == -1) return FAILURE;
				//printf("mkdir good\n");
				sprintf(cmdTemp,"mount %s -t ext2 /LWA_EXT/%s",deviceName,deviceName);
				if (issueShellCommand(cmdTemp, buffer, 2048) != -1) {
					devices[i].mounted=true;
					return SUCCESS;
				}*/
				//printf("mount not good\n");
			} else {
				devices[i].mounted=true;
				return SUCCESS;
			}
		}
	}
	// if we got here, then we didn't mount the drive, or an error occurred
	return FAILURE;
}

StatusCode Disk_Mount(Disk* disk){
	assert(disk);
	StatusCode sc = Disk_MountByName(disk->deviceName);
	int i = DiskFindIndexByName(disk->deviceName);
	if (i!=-1){
		*disk=devices[i];
		return sc;
	} else {
		bzero (disk,sizeof(Disk));
		return FAILURE;
	}
}




StatusCode Disk_UnmountByName(char* deviceName){
	int i = DiskFindIndexByName(deviceName);
	char buffer[2048];
	char cmdTemp [1024];
	if (i!=-1){
		if (!devices[i].mounted || devices[i].locked) return FAILURE;
		if ((devices[i].type == DRIVETYPE_EXTERNAL_PARTITION) || ((devices[i].type == DRIVETYPE_INTERNAL_MULTIDISK_DEVICE))){
			if (devices[i].type == DRIVETYPE_EXTERNAL_PARTITION){
				sprintf(cmdTemp,"umount %s",deviceName);
				if (issueShellCommand(cmdTemp, buffer, 2048) != -1) return SUCCESS;
			} else {
				devices[i].mounted=false;
				return SUCCESS;
			}
		}
	}
	// if we got here, then we didn't unmount the drive, or an error occurred
	return FAILURE;
}
StatusCode Disk_Unmount(Disk* disk){
	assert(disk);
	StatusCode sc = Disk_UnmountByName(disk->deviceName);
	int i = DiskFindIndexByName(disk->deviceName);
	if (i!=-1){
		*disk=devices[i];
		return sc;
	} else {
		bzero (disk,sizeof(Disk));
		return FAILURE;
	}
}

size_t Disk_GetArrayNumDisks(Disk* disk){
	assert(disk!=NULL);
	if (disk->type==DRIVETYPE_INTERNAL_MULTIDISK_DEVICE){
		char cmdTemp[1024];
		size_t disks=0;
		sprintf(cmdTemp,"mdadm --detail %s | grep \"Raid Devices\" | awk '{print $4}'",disk->deviceName);
		ShellGetNumber(cmdTemp,&disks);
		if (ShellGetNumber(cmdTemp,&disks)==FAILURE) {
			Log_Add("[DISK MANAGER] Error: couldn't get RAID array number of disks.\n");
			return 0l;
		}
		return disks;
	} else {
		return 1;
	}
}

size_t Disk_GetArraySize(Disk* disk){
	assert(disk!=NULL);
	return disk->bytesTotal;
	/*
	char cmdTemp[1024];
	size_t size=0;
	sprintf(cmdTemp,"sfdisk -s | grep \"%s\" | nawk '{print $2 \" * 1024\" }' | bc",disk->deviceName);
	if (ShellGetNumber(cmdTemp,&size)==FAILURE) {
		printf("Error: couldn't get RAID array size.\n");
		return 0l;
	}
	return size*1024;
	*/
}

size_t Disk_GetArrayChunkSize(Disk* disk){
	assert(disk!=NULL);
	if (disk->type==DRIVETYPE_INTERNAL_MULTIDISK_DEVICE){
		char cmdTemp[1024];
		size_t size=0;

		sprintf(cmdTemp,"mdadm --detail %s | grep \"Chunk Size\" | awk '{print $4}' | sed 's/K/\\*1024/' | bc",disk->deviceName);
		if (ShellGetNumber(cmdTemp,&size)==FAILURE) {
			Log_Add("[DISK MANAGER] Error: couldn't get RAID array size.\n");
			return 0l;
		}
		return size;
	} else {
		return 4096;
	}
}


int    DiskGetUsableCount(){
	int i,count=0;
	if (devicesCount<=0) return 0;
	for(i=0; i<devicesCount; i++){
		if (devices[i].type==DRIVETYPE_EXTERNAL_PARTITION){
			count++;
		}
	}
	return count;
}

char * DiskGetUsableName(int index){
	int i,count=0;
	if (devicesCount<=0) return NULL;
	for(i=0; i<devicesCount; i++){
		if (count==index && devices[i].type==DRIVETYPE_EXTERNAL_PARTITION){
			return devices[i].deviceName;
		}
		if (devices[i].type==DRIVETYPE_EXTERNAL_PARTITION){
			count++;
		}
	}
	return NULL;
}

size_t DiskGetUsableFreeSpace(int index){
	int i,count=0;
	if (devicesCount<=0) return 0;
	for(i=0; i<devicesCount; i++){
		if (count==index && devices[i].type==DRIVETYPE_EXTERNAL_PARTITION){
			Disk_UpdateFreeSpace(&devices[i]);
			return devices[i].bytesTotal - devices[i].bytesUsed;
		}
		if (devices[i].type==DRIVETYPE_EXTERNAL_PARTITION){
			count++;
		}
	}
	return 0;
}

size_t DiskGetUsableTotalSpace(int index){
	int i,count=0;
	if (devicesCount<=0) return 0;
	for(i=0; i<devicesCount; i++){
		if (count==index && devices[i].type==DRIVETYPE_EXTERNAL_PARTITION){
			return devices[i].bytesTotal;
		}
		if (devices[i].type==DRIVETYPE_EXTERNAL_PARTITION){
			count++;
		}
	}
	return 0;
}

StatusCode Disk_UpdateFreeSpace(Disk* disk){
	if (disk==NULL)  return FAILURE;
	int i = DiskFindIndexByName(disk->deviceName);
	if (i==-1) return FAILURE;
	char cmdTemp[8192];
	size_t blocksUsed=0;
	switch(devices[i].type){
		case DRIVETYPE_UNKNOWN					: //fall through
		case DRIVETYPE_SYSTEM_BLOCK_DEVICE		: //fall through
		case DRIVETYPE_SYSTEM_PARTITION			: //fall through
		case DRIVETYPE_SYSTEM_OPTICAL			: //fall through
		case DRIVETYPE_SYSTEM_FLOPPY			: //fall through
		case DRIVETYPE_INTERNAL_BLOCK_DEVICE	: //fall through
		case DRIVETYPE_EXTERNAL_BLOCK_DEVICE	: //fall through
		case DRIVETYPE_EXTERNAL_MULTIDISK_DEVICE: //fall through (note that external multidisk devices should already be re-classed as external partitions)
		default :
			break;
		case DRIVETYPE_INTERNAL_MULTIDISK_DEVICE:
			break;
		case DRIVETYPE_EXTERNAL_PARTITION		:
			sprintf(cmdTemp,"df | grep \"%s \" | nawk '{print $3}'",devices[i].deviceName);
			if (ShellGetNumber(cmdTemp, &blocksUsed)!=SUCCESS)
				Log_Add("[Disk] Warning: freespace determination for '%s' may be inaccurate",devices[i].deviceName);
			devices[i].bytesUsed=(blocksUsed * 1024);
			break;
	}
	return SUCCESS;

}





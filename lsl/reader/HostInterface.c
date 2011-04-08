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
 * HostInterface.c
 *
 *  Created on: Aug 13, 2009
 *      Author: chwolfe2
 */

#include "HostInterface.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>
int issueShellCommand(char * command, char* result, size_t resultSizeLimit){
	if (command==NULL || result==NULL || resultSizeLimit==0) return -1;
	FILE *pipedOutput = popen(command, "r");
	size_t count=0;
	if (pipedOutput){
		while (!feof(pipedOutput)){
			result[count]=fgetc ( pipedOutput );
			if (result[count]==EOF) {
				result[count]=0;
				break;
			}
			count++;
			if (count == resultSizeLimit-1){
				while (!feof(pipedOutput)) fgetc ( pipedOutput ); // dump any remaining characters
			}
		}
	}
	result[count]='\0';
	return pclose(pipedOutput);
}

StatusCode ShellGetNumber(char * command, size_t* number){
	char buffer[2048];
		if (issueShellCommand(command, buffer, 2048) == -1) {
			return FAILURE;
		}
		if (strlen(buffer)==0){
			return FAILURE;
		}
		*number=strtoul(buffer,NULL,10);
		return SUCCESS;
}

StatusCode ShellGetString(char * command, char* result, size_t resultSizeLimit){
	if (issueShellCommand(command, result, resultSizeLimit) == -1)
		return FAILURE;
	else
		return SUCCESS;
}


StatusCode  LaunchProcess(const char* command, int keepOutput, size_t resultSizeLimit, char * result, ProcessHandle** ph){
	if (command==NULL || ph==NULL) return FAILURE;
	if ( keepOutput && (resultSizeLimit==0 || result==NULL)) return FAILURE;
	*ph=(ProcessHandle*) malloc(sizeof(ProcessHandle));
	if (*ph==NULL) return FAILURE;
	(*ph)->fHandle			= popen(command, "r");
	(*ph)->resultSize		= 0;
	(*ph)->resultBuffer		= result;
	(*ph)->resultSizeMax	= resultSizeLimit;
	(*ph)->keepOutput 		= keepOutput;
	if ((*ph)->fHandle){
		return SUCCESS;
	} else {
		return FAILURE;
	}
}

StatusCode  CheckProcessDone(ProcessHandle* ph){
	assert(ph);
	int eofDetected=0;
	char ch;
	if (!ph) {
		printf("CheckProcessDone: bad pointer 'ph'.\n");
		return FAILURE;
	}
	if (!ph->fHandle){
		printf("CheckProcessDone: bad pointer 'ph->fHandle'.\n");
		return FAILURE;
	} else {
		if (!feof(ph->fHandle)){
			if (ph->keepOutput && (ph->resultSize < ph->resultSizeMax-1)){
				ph->resultBuffer[ph->resultSize]=fgetc ( ph->fHandle );
				if (ph->resultBuffer[ph->resultSize]==EOF) {
					ph->resultBuffer[ph->resultSize]=0;
					eofDetected=1;
				} else {
					//printf("%c",ph->resultBuffer[ph->resultSize]);
					ph->resultSize++;
				}
			} else {
				ch = fgetc ( ph->fHandle );
				ph->resultSize++;
				if (ch==EOF) eofDetected=1;
			}
		} else {
			eofDetected=1;
		}
	}
	if (eofDetected){
		if (ph->keepOutput){
			if (ph->resultSize < ph->resultSizeMax-1){
				ph->resultBuffer[ph->resultSize]=0;
			} else {
				ph->resultBuffer[ph->resultSizeMax-1]=0;
			}
		}
		return SUCCESS;
	} else {
		return NOT_READY;
	}

}
StatusCode  ReleaseProcessHandle(ProcessHandle* ph){
	assert(ph);
	if (ph->fHandle)
		pclose(ph->fHandle);
	free(ph);
	return SUCCESS;

}

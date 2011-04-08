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
 * HostInterface.h
 *
 *  Created on: Aug 13, 2009
 *      Author: chwolfe2
 */

#ifndef HOSTINTERFACE_H_
#define HOSTINTERFACE_H_
#include <stdlib.h>
#include <stdio.h>

#include "Defines.h"
typedef struct __ProcessHandle{
	FILE*  fHandle;
	char*  resultBuffer;
	size_t resultSizeMax;
	size_t resultSize;
	int    keepOutput; // do we record what the process dumps to stdout?
}ProcessHandle;

int issueShellCommand(char * command, char* result, size_t resultSizeLimit);
StatusCode ShellGetNumber(char * command, size_t* number);
StatusCode ShellGetString(char * command, char* result, size_t resultSizeLimit);

StatusCode  LaunchProcess(const char* command, int keepOutput, size_t resultSizeLimit, char * result, ProcessHandle** ph);
StatusCode  ReleaseProcessHandle(ProcessHandle* ph);
StatusCode  CheckProcessDone(ProcessHandle* ph);


#endif /* HOSTINTERFACE_H_ */

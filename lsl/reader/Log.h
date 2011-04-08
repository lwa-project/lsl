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
 * Log.h
 *
 *  Created on: Oct 29, 2009
 *      Author: chwolfe2
 */

#ifndef LOG_H_
#define LOG_H_

#include "Defines.h"
#include <stdio.h>
#include <stdarg.h>

#define LOGFILENAME "/LWA/database/logDb.gdbm"

StatusCode Log_Initialize(boolean clear);
StatusCode Log_Add(char* logentryFormat, ...);
StatusCode Log_GetCount(int * count);
StatusCode Log_GetEntry(int whichEntry, char* result);
StatusCode Log_DumpToFile(FILE * stream);

void Log_Close();


#endif /* LOG_H_ */

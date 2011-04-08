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
 * Defines.h
 *
 *  Created on: Oct 18, 2009
 *      Author: chwolfe2
 */

#ifndef DEFINES_H_
#define DEFINES_H_
#pragma pack(1)

typedef int boolean;
#define true 1
#define false 0
#define TRUE 1
#define FALSE 0

#define forever 1
#define never 0

typedef int StatusCode;
#define SUCCESS 0
#define FAILURE 1
#define NOT_READY 2

#define statusString(s) ((s==SUCCESS) ? "SUCCESS" : ((s==FAILURE) ? "FAILURE" : "NOT_READY"))

#define EXIT_CODE_FAIL           0x00
#define EXIT_CODE_FLUSH_DATA     0x01
#define EXIT_CODE_FLUSH_LOG      0x08
#define EXIT_CODE_FLUSH_SCHEDULE 0x04
#define EXIT_CODE_FLUSH_CONFIG   0x02
#define EXIT_CODE_SHUTDOWN       0x10
#define EXIT_CODE_REBOOT         0x20


#endif /* DEFINES_H_ */

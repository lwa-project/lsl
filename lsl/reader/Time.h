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
 * Time.h
 *
 *  Created on: Oct 25, 2009
 *      Author: chwolfe2
 */

#ifndef TIME_H_
#define TIME_H_
#include <stdlib.h>
#include "Defines.h"
#define MILLISECONDS_PER_DAY (1000l*60l*60l*24l)
typedef struct __timestamp {
	size_t MPM;
	size_t MJD;
	size_t hour;
	size_t minute;
	size_t second;
	size_t millisecond;
	size_t year;
	size_t month;
	size_t day;
} TimeStamp;

TimeStamp getTimestamp(void);
size_t    getMPM(void);
size_t    getMJD(void);
size_t    getElapsedMilliseconds(TimeStamp start, TimeStamp stop);
boolean   before (size_t MJDa, size_t MPMa, size_t MJDb, size_t MPMb); // a is before b
boolean   after (size_t MJDa, size_t MPMa, size_t MJDb, size_t MPMb);  // a is after b
boolean   equal (size_t MJDa, size_t MPMa, size_t MJDb, size_t MPMb);  // a is equal to b
void      addTime(size_t MJD, size_t MPM, size_t milliseconds, size_t* resultMJD, size_t* resultMPM);
void      subTime(size_t MJD, size_t MPM, size_t milliseconds, size_t* resultMJD, size_t* resultMPM);
char* HumanTime(void);
int compareDates(size_t mjda, size_t mpma, size_t mjdb, size_t mpmb);

#endif /* TIME_H_ */

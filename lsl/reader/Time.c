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
 * Time.c
 *
 *  Created on: Oct 25, 2009
 *      Author: chwolfe2
 */

#include "Time.h"
#include <time.h>
#include <sys/time.h>
#include <stdio.h>

size_t getMPM(){
	struct timeval dt;
	gettimeofday(&dt, NULL);
	struct tm timeStruct;
	gmtime_r(&dt.tv_sec, &timeStruct);
	size_t hour        = timeStruct.tm_hour;
	size_t minute      = timeStruct.tm_min;
	size_t second      = timeStruct.tm_sec;
	size_t millisecond = dt.tv_usec / 1000;
	// compute MPM
    return (hour*3600 + minute*60 + second)*1000 + millisecond;
}

size_t getMJD(){
	struct timeval dt;
	gettimeofday(&dt, NULL);
	struct tm timeStruct;
	gmtime_r(&dt.tv_sec, &timeStruct);
	size_t year        = timeStruct.tm_year +1900;
	size_t month       = timeStruct.tm_mon +1;
	size_t day         = timeStruct.tm_mday;
	// compute MJD
    // adapted from http://paste.lisp.org/display/73536
    // can check result using http://www.csgnetwork.com/julianmodifdateconv.html
    size_t a = (14 - month) / 12;
    size_t y = year + 4800 - a;
    size_t m = month + (12 * a) - 3;
    size_t p = day + (((153 * m) + 2) / 5) + (365 * y);
    size_t q = (y / 4) - (y / 100) + (y / 400) - 32045;
    return ( (p+q) - 2400000.5);

}

TimeStamp getTimestamp(){
	TimeStamp rv;
	struct timeval dt;
	gettimeofday(&dt, NULL);
	struct tm timeStruct;
	gmtime_r(&dt.tv_sec, &timeStruct);
	rv.hour        = timeStruct.tm_hour;
	rv.minute      = timeStruct.tm_min;
	rv.second      = timeStruct.tm_sec;
	rv.millisecond = dt.tv_usec / 1000;
	rv.year        = timeStruct.tm_year +1900;
	rv.month       = timeStruct.tm_mon +1;
	rv.day         = timeStruct.tm_mday;
	// compute MPM
	rv.MPM =  (rv.hour*3600 + rv.minute*60 + rv.second)*1000 + rv.millisecond;
	size_t a = (14 - rv.month) / 12;
	size_t y = rv.year + 4800 - a;
	size_t m = rv.month + (12 * a) - 3;
	size_t p = rv.day + (((153 * m) + 2) / 5) + (365 * y);
	size_t q = (y / 4) - (y / 100) + (y / 400) - 32045;
	// compute MJD
	rv.MJD = ( (p+q) - 2400000.5);
	return rv;
}
size_t getElapsedMilliseconds(TimeStamp start, TimeStamp stop){
	if (start.MJD==stop.MJD){
		return (stop.MPM-start.MPM);
	} else {
		return ((stop.MJD-start.MJD)*1000l*60l*60l*20l)+(stop.MPM-start.MPM);
	}
}
boolean before (size_t MJDa, size_t MPMa, size_t MJDb, size_t MPMb){ // a is before b
	if (MJDa<MJDb) return true;
	if (MJDa>MJDb) return false;
	if (MJDa==MJDb){
		return (MPMa<MPMb);
	}
	return false;
}
boolean after (size_t MJDa, size_t MPMa, size_t MJDb, size_t MPMb){  // a is after b
	if (MJDa>MJDb) return true;
	if (MJDa<MJDb) return false;
	if (MJDa==MJDb){
		return (MPMa>MPMb);
	}
	return false;

}
boolean equal (size_t MJDa, size_t MPMa, size_t MJDb, size_t MPMb){  // a is equal to b
	return ((MPMa==MPMb) && (MJDa==MJDb));

}
void      addTime(size_t MJD, size_t MPM, size_t milliseconds, size_t* resultMJD, size_t* resultMPM){
	*resultMJD = MJD;
	*resultMPM = MPM + milliseconds;
	while (*resultMPM > MILLISECONDS_PER_DAY){
		(*resultMJD)++;
		*resultMPM -= MILLISECONDS_PER_DAY;
	}
}
void      subTime(size_t MJD, size_t MPM, size_t milliseconds, size_t* resultMJD, size_t* resultMPM){
	*resultMJD = MJD;
	*resultMPM = MPM;
	while (milliseconds > MILLISECONDS_PER_DAY){
		(*resultMJD)--;
		milliseconds -= MILLISECONDS_PER_DAY;
	}
	if (milliseconds > *resultMPM){
		(*resultMJD)--;
		*resultMPM+=*resultMPM;
	}
	*resultMPM -= milliseconds;
	//printf("%lu:%lu - %lu = %lu:%lu\n",MJD,MPM,milliseconds,*resultMJD,*resultMPM);
	//exit(EXIT_CODE_FAIL);
}
char HumanTimeBuffer[1024];
char* mnames[]={
		"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"
};
char* HumanTime(){
	TimeStamp ts=getTimestamp();
	sprintf(HumanTimeBuffer, "%02lu %3s %04lu %02lu:%02lu:%02lu:%03lu", ts.day,mnames[ts.month-1],ts.year,ts.hour,ts.minute,ts.second,ts.millisecond);
	return HumanTimeBuffer;
}
int compareDates(size_t mjda, size_t mpma, size_t mjdb, size_t mpmb){
	if (mjda<mjdb) return -1;
	if (mjda>mjdb) return  1;
	if (mpma<mpmb) return -1;
	if (mpma>mpmb) return  1;
	return 0;
}


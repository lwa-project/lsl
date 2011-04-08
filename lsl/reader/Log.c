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
 * Log.c
 *
 *  Created on: Oct 29, 2009
 *      Author: chwolfe2
 */

#include "Persistence.h"
#include "Log.h"
#include "Time.h"


Database logDb;
boolean isLogOpen=false;


StatusCode Log_Initialize(boolean clear){
	if (!isLogOpen){
		if (!clear){
			StatusCode sc = Persistence_Open(&logDb,LOGFILENAME);
			if (sc!=SUCCESS){
				sc = Persistence_Create(LOGFILENAME,NULL);
				if (sc!=SUCCESS){
					return sc;
				} else {
					sc=Persistence_Open(&logDb,LOGFILENAME);
					if (sc==SUCCESS){
						atexit(Log_Close);
						isLogOpen=true;
					}
					sc = Persistence_List_Create(&logDb,"logitems");
					return sc;
				}
			} else {
				Log_DumpToFile(stdout);
				atexit(Log_Close);
				isLogOpen=true;
				return sc;
			}
		} else {
			StatusCode sc = Persistence_Create(LOGFILENAME,NULL);
			if (sc!=SUCCESS){
				return sc;
			} else {
				sc=Persistence_Open(&logDb,LOGFILENAME);
				if (sc==SUCCESS){
					atexit(Log_Close);
					isLogOpen=true;
				}
				sc = Persistence_List_Create(&logDb,"logitems");
				return sc;
			}
		}
	} else {
		return FAILURE;
	}
}

StatusCode Log_Add(char* logentryFormat, ...){
	char logentry[8173];
	va_list args;
	va_start (args, logentryFormat);
	vsprintf (logentry,logentryFormat, args);
	//vprintf(logentryFormat, args);fflush(stdout);
	va_end (args);
	char buffer[8192];
	sprintf(buffer,"%06lu-%09lu : %s",getMJD(),getMPM(),logentry);
	//printf("%s\n",buffer);
	if (!isLogOpen) return FAILURE;
	return Persistence_List_AddItem(&logDb,"logitems",POSITION_FRONT,buffer);

}
StatusCode Log_GetCount(int * count){
	if (!isLogOpen) return FAILURE;
	return Persistence_List_GetCount(&logDb,"logitems",count);
}
StatusCode Log_GetEntry(int whichEntry, char* result){
	if (!isLogOpen) return FAILURE;
	int count;
	StatusCode sc = Persistence_List_GetCount(&logDb,"logitems",&count);
	if (sc!=SUCCESS) return sc;
	if ((whichEntry >= count)||(whichEntry <= -3)){
		sprintf(result,"<INVALID LOG INDEX SPECIFIED>");
		return FAILURE;
	} else {
		return Persistence_List_GetItem(&logDb,"logitems",whichEntry,result);
	}

}
StatusCode Log_DumpToFile(FILE * stream){
	char buffer[8192];
	if (stream==NULL) return FAILURE;
	int count;
	StatusCode sc = Persistence_List_GetCount(&logDb,"logitems",&count);
	if (sc!=SUCCESS) return sc;
	fprintf(stream, "   --------------------- Log Contents Dump ---------------------\n");
	int i=0;
	for (i=0;i<count;i++){
		sc = Persistence_List_GetItem(&logDb,"logitems",i,buffer);
		if (sc!=SUCCESS) return sc;
		fprintf(stream,"[*** LOG DUMP ***] %s\n",buffer);
	}
	fprintf(stream, "   ------------------------ End of Log  ------------------------\n");
	return SUCCESS;
}

void Log_Close(){
	if (isLogOpen){
		isLogOpen=false;
		Persistence_Close(&logDb);
	}
}

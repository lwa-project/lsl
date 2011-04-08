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
 * Persistence.c
 *
 *  Created on: Oct 28, 2009
 *      Author: chwolfe2
 */

#include "Persistence.h"
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include "Log.h"

char FirstNonWS(char* str){
	while(str && (*str)){
		if ((*str==' ') || (*str=='\t') || (*str=='\n')){
			str++;
		} else {
			return *str;
		}
	}
	return '\0';
}
void fatalFunc(void){
	printf("HIT FATAL FUNCTION\n");
	exit(-1);
}
StatusCode Persistence_Create(					 char * fileName, char * defaultsFileName ){
	char buffer[65536];
	Database database;
	bzero(&database,sizeof(Database));
	database = gdbm_open(fileName, 8192, GDBM_NEWDB , 0666, fatalFunc);
	if (gdbm_errno != GDBM_NO_ERROR){
		Log_Add("[DATABASE] gdbm reported error: %d:'%s'",gdbm_errno,gdbm_strerror(gdbm_errno));
		exit(EXIT_CODE_FAIL);
	} else {
		Log_Add("[DATABASE] Database File Created: '%s'",fileName);
	}
	if (defaultsFileName==NULL){
		gdbm_sync(database);
		gdbm_close(database);
		return SUCCESS;
	}
	Log_Add("[DATABASE] Opening default values file: '%s'", defaultsFileName);
	FILE * fp = fopen (defaultsFileName, "r" );
	if (fp==NULL){
		Log_Add("[DATABASE] Error: defaults file '%s' not found",fileName);
		return FAILURE;
	}
	Log_Add("[DATABASE] Now filling database with default values.");
	char VariableName[40];
	char VariableValue[65535];
	char * token;
	datum key,value;
	int result;
	size_t lineCount=0;
	size_t recordCount=0;
	//int i=0;
	//int l=0;
	while (fgets(buffer,65536,fp)){
		//l=strlen(buffer);
		/*for (i=0;i<l;i++){
			if (buffer[i]=='\n') buffer[i] = ' ';
		}*/
		lineCount++;
		//printf("Processing '%s'\n",buffer);
		if (FirstNonWS(buffer)!='#'){
			token = strtok(buffer," \t\n");
			if (token!=NULL && strlen(token)>0){
				VariableName[0]='\0';
				strcpy(VariableName, token);
				//printf("VariableName '%s'\n",VariableName);
				token = strtok(NULL,"\t\n");
				if (token!=NULL && strlen(token)>0){
					VariableValue[0]='\0';
					strcpy(VariableValue, token);
					//printf("VariableValue '%s'\n",VariableValue);
					key.dptr=VariableName;
					key.dsize=strlen(VariableName)+1;
					value.dptr=VariableValue;
					value.dsize=strlen(VariableValue)+1;
					result = gdbm_store(database,key,value,GDBM_INSERT);
					recordCount++;
					if (result != 0 ){
						Log_Add("[DATABASE] gdbm insert reported error: %d:'%s'",gdbm_errno,gdbm_strerror(gdbm_errno));
						fclose(fp);
						gdbm_sync(database);
						gdbm_close(database);
						return FAILURE;
					}
				}
			}
		}
	}
	fclose(fp);
	Persistence_Close(&database);
	Log_Add("[DATABASE] Finished initializing database. \tProcessed %lu lines \t Added %lu records",lineCount,recordCount);
	return SUCCESS;

}
void Persistence_ListAll(Database* database){
	datum key,value;
	key = gdbm_firstkey(*database);
	while (key.dptr){
		value=gdbm_fetch(*database,key);
		if (value.dptr!=NULL){
			Log_Add("[DATABASE] DB: '%s' => '%s'", key.dptr, value.dptr);
			free(value.dptr);
		} else {
			Log_Add("[DATABASE] DB: '%s' => <<< VALUE MISSING >>>", key.dptr);
		}
		key = gdbm_nextkey(*database,key);
	}
}
StatusCode Persistence_Open(					Database* database, char * fileName){
	bzero(database,sizeof(Database));
	*database = gdbm_open(fileName, 262144, GDBM_WRITER , 0666, fatalFunc);
	if (gdbm_errno != GDBM_NO_ERROR){
		Log_Add("[DATABASE] gdbm reported warning/error: %d:'%s'",gdbm_errno,gdbm_strerror(gdbm_errno));
		return FAILURE ;
	} else {
		Log_Add("[DATABASE] Database File Opened: '%s'",fileName);
		//printf("AT OPEN DB, ADDRESS is %lx\n",(size_t) database);
		//int i = 0;
		//printf("DATABASE FILE DATA IN BYTES\n");
		//for(i=0;i<(10*sizeof(int));i++){
		//	printf("%.2hx ",(((char*)database)[i]&0xff));
		//	if ((i%10) == 9) printf("\n");
		//}
		//printf("---------------------------\n");
		//printf("DATABASE FILE DATA IN ints\n");
		//for(i=0;i<10;i++){
		//	printf("%20d  (%20x)\n",((int*)database)[i],((int*)database)[i]);
		//}
		//printf("---------------------------\n");

		//Persistence_ListAll(database);
		return SUCCESS;
	}
}
StatusCode Persistence_Close(					Database* database){
	//printf("AT CLOSE DB, ADDRESS is %lx\n",(size_t) database);
	//int i = 0;
	//printf("DATABASE FILE DATA IN BYTES\n");
	//for(i=0;i<(10*sizeof(int));i++){
	//	printf("%.2hx ",(((char*)database)[i]&0xff));
	//	if ((i%10) == 9) printf("\n");
	//}
	//printf("---------------------------\n");
	//printf("DATABASE FILE DATA IN ints\n");
	//for(i=0;i<10;i++){
	//	printf("%20d  (%20x)\n",((int*)database)[i],((int*)database)[i]);
	//}
	//printf("---------------------------\n");

	gdbm_sync(*database);
	gdbm_reorganize(*database);
	gdbm_close(*database);
	return SUCCESS;
}
StatusCode Persistence_DeleteVariable(			Database* database, char* variableName){
	datum key;
	key.dptr=variableName;
	key.dsize= strlen(variableName)+1;
	int res=gdbm_delete(*database,key);
	if (res!=0){
		return FAILURE;
	} else
		return SUCCESS;
}
void Unescape(char* escaped, char* dest){
	int i=0;
	char tmp[4];
	assert(escaped!=NULL);
	assert(dest!=NULL);
	while(*escaped){
		if (*escaped=='\\'){
			for( i=0;i<3;i++) tmp[i]=escaped[i];
			*escaped = (char) atoi (tmp);
			dest++;
			escaped++;
		} else {
			*dest = *escaped;
			dest++;
			escaped++;
		}
	}
	*dest='\0';

}
void Escape(char* unescaped, char* dest){
	assert(unescaped!=NULL);
	assert(dest!=NULL);
	char * s=unescaped;
	char * d=dest;
	while(*s){
		if ((!isprint(*s)) || (*s=='\\')){
			*d='\\';
			d++;
			char tmp[4];
			sprintf(tmp,"%03i",(int)*s);
			strcpy(d,tmp);
			d+=strlen(tmp);
			s++;
		} else {
			*d=*s;
			d++;
			s++;
		}
	}
	*d='\0';
}

StatusCode Persistence_Get_String(				Database* database, char* variableName, char*   result){
	assert(variableName!=NULL);
	datum key;
	key.dptr=variableName;
	key.dsize= strlen(variableName)+1;

	datum value;
	if (gdbm_exists(*database, key)){
		value = gdbm_fetch(*database, key);

		if (value.dptr!=NULL){
			Unescape(value.dptr,result);
			free(value.dptr);
			return SUCCESS;
		} else
			return FAILURE;
	} else {
		return FAILURE;
	}
}

StatusCode Persistence_Get_UnsignedLong(		Database* database, char* variableName, size_t* result){

	//printf("Persistence_Get_UnsignedLong (db, '%s', %lu)", variableName, (size_t)result);
	char buffer[65536];
	StatusCode s = Persistence_Get_String(database, variableName, buffer);
	if (s==SUCCESS){
		*result = strtoul(buffer, NULL, 10);
		return SUCCESS;
	} else
		return FAILURE;
}

StatusCode Persistence_Get_Int(					Database* database, char* variableName, int*    result){
	char buffer[65536];
	StatusCode s = Persistence_Get_String(database, variableName, buffer);
	if (s==SUCCESS){
		*result = atoi(buffer);
		return SUCCESS;
	} else
		return FAILURE;
}

StatusCode Persistence_Set_String(				Database* database, char* variableName, char*   value){
	char buffer[65536];
	Escape(value,buffer);
	datum key;
	key.dptr=variableName;
	key.dsize= strlen(variableName)+1;
	datum val;
	val.dptr=buffer;
	val.dsize =strlen(buffer)+1;
	int res=gdbm_store(*database,key,val,GDBM_REPLACE);
	if (res!=0){
		return FAILURE;
	} else
		return SUCCESS;

}

StatusCode Persistence_Set_UnsignedLong(		Database* database, char* variableName, size_t value){
	char buffer[65536];
	sprintf(buffer,"%lu",value);
	return Persistence_Set_String(database,variableName,buffer);
}

StatusCode Persistence_Set_Int(					Database* database, char* variableName, int    value){
	char buffer[65536];
	sprintf(buffer,"%d",value);
	return Persistence_Set_String(database,variableName,buffer);
}

char * mkArrName(char * variableName, int index){
	static char buffer[512];
	sprintf(buffer,"%s[%d]",variableName, index);
	return buffer;
}

StatusCode Persistence_GetArray_String(			Database* database, char* variableName, int index, char*   result){
	return Persistence_Get_String(database, mkArrName(variableName,index), result);
}

StatusCode Persistence_GetArray_UnsignedLong(	Database* database, char* variableName, int index, size_t* result){
	return Persistence_Get_UnsignedLong(database, mkArrName(variableName,index), result);
}

StatusCode Persistence_GetArray_Int(			Database* database, char* variableName, int index, int*    result){
	return Persistence_Get_Int(database, mkArrName(variableName,index), result);
}

StatusCode Persistence_SetArray_String(			Database* database, char* variableName, int index, char*   value){
	return Persistence_Set_String(database, mkArrName(variableName,index), value);
}

StatusCode Persistence_SetArray_UnsignedLong(	Database* database, char* variableName, int index, size_t value){
	return Persistence_Set_UnsignedLong(database, mkArrName(variableName,index), value);
}

StatusCode Persistence_SetArray_Int(			Database* database, char* variableName, int index, int    value){
	return Persistence_Set_Int(database, mkArrName(variableName,index), value);
}

StatusCode Persistence_List_Exists(Database* database, char* listName){
	size_t ignored;
	char buffer[65536];
	sprintf(buffer,"%s.count",listName);
	return Persistence_Get_UnsignedLong(database,buffer,&ignored);
}

StatusCode Persistence_List_Create(Database* database, char* listName){
	if (Persistence_List_Exists(database,listName)!=SUCCESS){
	char buffer[65536];
		sprintf(buffer,"%s.count",listName);
		return Persistence_Set_UnsignedLong(database,buffer,0);
	}
	return FAILURE;
}

StatusCode Persistence_List_AddItem(Database* database, char* listName, int position, char * value){
	StatusCode sc;
	int count=0;
	char copytemp[65536];
	bzero(copytemp,65536);
	sc = Persistence_List_GetCount(database,listName,&count);
	if (sc!=SUCCESS) return sc;
	if (position == POSITION_FRONT) 	position = 0;
	if (position == POSITION_BACK) 		position = count;
	if (position < 0) position=0;
	if (position > count) position=count;

	size_t i=0;
	if (position < count){
		for(i=count; i>position; i--){
			sc = Persistence_GetArray_String(database,listName,i-1,copytemp);
			if (sc!=SUCCESS) return sc;
			sc = Persistence_SetArray_String(database,listName,i,copytemp);
			if (sc!=SUCCESS) return sc;
		}
	}
	sc = Persistence_SetArray_String(database,listName,position,value);
	if (sc!=SUCCESS) return sc;
	return Persistence_List_SetCount(database,listName,count+1);
}

StatusCode Persistence_List_RemoveItem(Database* database, char* listName, int position){
	StatusCode sc;
	int count;
	char buffer[1024];
	char copytemp[65536];
	sprintf(buffer,"%s.count",listName);
	sc = Persistence_List_GetCount(database,listName,&count);
	if (sc!=SUCCESS) return sc;
	if (position == POSITION_FRONT) 	position = 0;
	if (position == POSITION_BACK) 		position = count-1;
	if (position < 0) position=0;
	if (position > count) position=count;
	size_t i=0;
	if (position < count){
		for(i=position; i<count-1; i++){
			sc = Persistence_GetArray_String(database, listName, i+1, copytemp);
			if (sc!=SUCCESS) return sc;
			sc = Persistence_SetArray_String(database,listName,i,copytemp);
			if (sc!=SUCCESS) return sc;
		}
		sprintf(buffer,"%s[%d]",listName,count-1);
		sc = Persistence_DeleteVariable(database,buffer);
		if (sc!=SUCCESS) return sc;
		return Persistence_List_SetCount(database,listName,count-1);
	} else {
		return FAILURE;
	}
}
StatusCode Persistence_List_RemoveItemBinary(Database* database, char* listName, int position){
	StatusCode sc;
	int count;
	char buffer[1024];
	char copytemp[65536];
	sprintf(buffer,"%s.count",listName);
	sc = Persistence_List_GetCount(database,listName,&count);

	if (sc!=SUCCESS) return sc;
	if (position == POSITION_FRONT) 	position = 0;
	if (position == POSITION_BACK) 		position = count-1;
	if (position < 0) position=0;
	if (position > count) position=count;

	int i=0;
	size_t bytes;
	if (position < count){
		for(i=position; i<count-1; i++){

			sc = Persistence_GetArray_Binary(database, listName, i+1, copytemp, &bytes);
			if (sc!=SUCCESS) return sc;
			sc = Persistence_SetArray_Binary(database,listName,i,copytemp, bytes);
			if (sc!=SUCCESS) return sc;
		}

		sprintf(buffer,"%s[%d]",listName,count-1);
		sc = Persistence_DeleteVariable(database,buffer);
		if (sc!=SUCCESS) return sc;


		return Persistence_List_SetCount(database,listName,count-1);
	} else {
		return FAILURE;
	}
}
StatusCode Persistence_List_AddItemBinary(Database* database, char* listName, int position, char * value, size_t size){
	StatusCode sc;
	int count;
	char copytemp[65536];
	sc = Persistence_List_GetCount(database,listName,&count);
	if (sc!=SUCCESS) return sc;
	if (position == POSITION_FRONT) 	position = 0;
	if (position == POSITION_BACK) 		position = count;
	if (position < 0) position=0;
	if (position > count) position=count;

	size_t i=0;
	size_t bytes;
	if (position < count){
		for(i=count; i>position; i--){
			sc = Persistence_GetArray_Binary(database,listName,i-1,copytemp, &bytes);
			if (sc!=SUCCESS) return sc;
			sc = Persistence_SetArray_Binary(database,listName,i,copytemp, bytes);
			if (sc!=SUCCESS) return sc;
		}
	}
	sc = Persistence_SetArray_Binary(database,listName,position,value,size);
	if (sc!=SUCCESS) return sc;
	return Persistence_List_SetCount(database,listName,count+1);
}

StatusCode Persistence_List_GetItemBinary(Database* database, char* listName, int position, char* result, size_t* size){
	StatusCode sc;
	int count;
	sc = Persistence_List_GetCount(database,listName,&count);
	if (sc!=SUCCESS) return sc;
	if (position == POSITION_FRONT) 	position = 0;
	if (position == POSITION_BACK) 		position = count-1;
	if (position < 0) position=0;
	if (position > count) position=count;
	if (position < count){
		return Persistence_GetArray_Binary(database, listName, position, result, size);
	} else {
		return FAILURE;
	}
}

StatusCode Persistence_List_GetItem(Database* database, char* listName, int position, char* result){
	StatusCode sc;
	int count;
	sc = Persistence_List_GetCount(database,listName,&count);
	if (sc!=SUCCESS) return sc;
	if (position == POSITION_FRONT) 	position = 0;
	if (position == POSITION_BACK) 		position = count-1;
	if (position < 0) position=0;
	if (position > count) position=count;

	if (position < count){
		return Persistence_GetArray_String(database, listName, position, result);
	} else {
		return FAILURE;
	}
}

StatusCode Persistence_List_GetCount(Database* database, char* listName, int* count){
	char buffer[1024];
	sprintf(buffer,"%s.count",listName);
	return Persistence_Get_Int(database,buffer,count);
}

StatusCode Persistence_List_SetCount(Database* database, char* listName, int count){
	char buffer[1024];
	sprintf(buffer,"%s.count",listName);
	//printf("setting %s.count to %d\n",listName, count);
	return Persistence_Set_Int(database,buffer,count);
}
StatusCode Persistence_List_Delete(Database* database, char* listName){
	StatusCode sc;
	int count;
	char buffer[1024];
	sprintf(buffer,"%s.count",listName);
	sc = Persistence_List_GetCount(database,listName,&count);
	if (sc!=SUCCESS) return sc;
	size_t i=0;
	for(i=0; i<count; i++){
		sprintf(buffer,"%s[%lu]",listName,i);
		sc = Persistence_DeleteVariable(database,buffer);
		if (sc!=SUCCESS) return sc;
	}
	sprintf(buffer,"%s.count",listName);
	return Persistence_DeleteVariable(database,buffer);
}



StatusCode Persistence_Get_Binary(				Database* database, char* variableName, char*   result, size_t* datalen){
	assert(result!=NULL);
	assert(variableName!=NULL);
	datum key;
	key.dptr=variableName;
	key.dsize= strlen(variableName)+1;
	datum value;
	value = gdbm_fetch(*database, key);
	if (value.dptr!=NULL){
		memcpy(result, value.dptr, value.dsize);
		*datalen = value.dsize;
		free(value.dptr);
		return SUCCESS;
	} else
		return FAILURE;
}
StatusCode Persistence_Set_Binary(				Database* database, char* variableName, char*   value,  size_t datalen){
	datum key;
	datum val;
	key.dptr  = variableName;
	key.dsize = strlen(variableName)+1;
	val.dptr  = value;
	val.dsize = datalen;
	int res=gdbm_store(*database,key,val,GDBM_REPLACE);
	if (res!=0){
		return FAILURE;
	} else{
		return SUCCESS;
	}
}
StatusCode Persistence_GetArray_Binary(			Database* database, char* variableName, int index, char*   result, size_t* datalen){
	return Persistence_Get_Binary(database, mkArrName(variableName,index), result, datalen);
}
StatusCode Persistence_SetArray_Binary(			Database* database, char* variableName, int index, char*   value, size_t datalen){
	return Persistence_Set_Binary(database, mkArrName(variableName,index), value, datalen);
}



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
 * Persistence.h
 *
 *  Created on: Oct 28, 2009
 *      Author: chwolfe2
 */

#ifndef PERSISTENCE_H_
#define PERSISTENCE_H_
#include <gdbm.h>
#include <stdlib.h>
#include "Defines.h"

typedef GDBM_FILE Database;

#define POSITION_FRONT -1
#define POSITION_BACK  -2

StatusCode Persistence_Create(					 char * fileName, char * defaultsFileName );
StatusCode Persistence_Open(					Database* database, char * fileName);
StatusCode Persistence_Close(					Database* database);
StatusCode Persistence_Get_String(				Database* database, char* variableName, char*   result);
StatusCode Persistence_Get_UnsignedLong(		Database* database, char* variableName, size_t* result);
StatusCode Persistence_Get_Int(					Database* database, char* variableName, int*    result);
StatusCode Persistence_Get_Binary(				Database* database, char* variableName, char*   result, size_t* datalen);
StatusCode Persistence_Set_String(				Database* database, char* variableName, char*   value);
StatusCode Persistence_Set_UnsignedLong(		Database* database, char* variableName, size_t value);
StatusCode Persistence_Set_Int(					Database* database, char* variableName, int    value);
StatusCode Persistence_Set_Binary(				Database* database, char* variableName, char*   value,  size_t datalen);
StatusCode Persistence_GetArray_String(			Database* database, char* variableName, int index, char*   result);
StatusCode Persistence_GetArray_UnsignedLong(	Database* database, char* variableName, int index, size_t* result);
StatusCode Persistence_GetArray_Int(			Database* database, char* variableName, int index, int*    result);
StatusCode Persistence_GetArray_Binary(			Database* database, char* variableName, int index, char*   result, size_t* datalen);
StatusCode Persistence_SetArray_String(			Database* database, char* variableName, int index, char*   value);
StatusCode Persistence_SetArray_UnsignedLong(	Database* database, char* variableName, int index, size_t value);
StatusCode Persistence_SetArray_Int(			Database* database, char* variableName, int index, int   value);
StatusCode Persistence_SetArray_Binary(			Database* database, char* variableName, int index, char*   value, size_t datalen);
StatusCode Persistence_DeleteVariable(			Database* database, char* variableName);


StatusCode Persistence_List_Create(Database* database, char* listName);
StatusCode Persistence_List_Exists(Database* database, char* listName);
StatusCode Persistence_List_Delete(Database* database, char* listName);
StatusCode Persistence_List_GetCount(Database* database, char* listName, int* count);
StatusCode Persistence_List_SetCount(Database* database, char* listName, int count);
StatusCode Persistence_List_AddItem(Database* database, char* listName, int position, char * value);
StatusCode Persistence_List_GetItem(Database* database, char* listName, int position, char* result);
StatusCode Persistence_List_AddItemBinary(Database* database, char* listName, int position, char * value, size_t size);
StatusCode Persistence_List_GetItemBinary(Database* database, char* listName, int position, char* result, size_t* size);
StatusCode Persistence_List_RemoveItem(Database* database, char* listName, int position);
StatusCode Persistence_List_RemoveItemBinary(Database* database, char* listName, int position);
#endif /* PERSISTENCE_H_ */

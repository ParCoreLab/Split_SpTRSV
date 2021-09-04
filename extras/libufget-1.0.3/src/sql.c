/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright (C) Martin Koehler, 2015
 */

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#include <sqlite3.h>

#include "libufget.h"
#include "uf_internal.h"



/*-----------------------------------------------------------------------------
 *  Execute SQL statement 
 *-----------------------------------------------------------------------------*/
int sql_execute(sqlite3 *db, const char* sql, ...)
{
    char *err, *tmp;
    va_list ap;
    va_start(ap, sql);
    tmp = sqlite3_vmprintf(sql, ap);
#ifdef DEEP_DEBUG
    printf("statement: %s\n", tmp);
#endif 
    va_end(ap);
    int rc = sqlite3_exec(db, tmp, NULL, NULL, &err);
    if(rc != SQLITE_OK) {
        if (err != NULL) {
            fprintf(stdout, "sql_execute() : Error %i : %s\n", rc, err);
            sqlite3_free(err);
        }
    }
    sqlite3_free(tmp);
    return rc;
}


/*-----------------------------------------------------------------------------
 *  Get Table
 *-----------------------------------------------------------------------------*/
int sql_get_table(sqlite3 *db, char ***resultp, int *nrow, int *ncol, const char *sql, ...)
{
    char *err, *tmp;
    va_list ap;
    va_start(ap, sql);
    tmp = sqlite3_vmprintf(sql, ap);
    va_end(ap);
    int rc = sqlite3_get_table(db, tmp, resultp, nrow, ncol, &err);
    if(rc != SQLITE_OK) {
        if (err != NULL) {
            fprintf(stdout, "sql_get_table() : Error %i : %s\n", rc, err);
            sqlite3_free(err);
        }
    }
    sqlite3_free(tmp);
    return rc;
}


/*-----------------------------------------------------------------------------
 *  Check if table exists
 *-----------------------------------------------------------------------------*/
int sql_table_exists(sqlite3 * db, const char *table_name) 
{
    char **result; 
    int nrow, ncol; 
    int rc; 

    rc = sql_get_table(db, &result, &nrow, &ncol, "SELECT name FROM sqlite_master WHERE type='table' AND name='%s'", table_name);
    sqlite3_free_table(result); 

    if ( rc != SQLITE_OK ) {
        return 0; 
    }
    if ( nrow == 1 ) 
        return 1; 
    else 
        return 0; 
}


/*-----------------------------------------------------------------------------
 *  Sql Get Integer
 *-----------------------------------------------------------------------------*/
int sql_get_long(sqlite3 * db, const char * sql, ...) 
{
    char *err, *tmp;
    va_list ap;
    int nrow, ncol; 
    char **result; 
    long ret; 

    va_start(ap, sql);
    tmp = sqlite3_vmprintf(sql, ap);
    va_end(ap);
    
    int rc = sqlite3_get_table(db, tmp, &result, &nrow, &ncol, &err);
    if(rc != SQLITE_OK) {
        if (err != NULL) {
            fprintf(stdout, "sql_get_long() : Error %i : %s\n", rc, err);
            sqlite3_free(err);
        }
        ret = 0; 
    } else {
        if (nrow == 0 ) 
            ret = 0; 
        else {
            ret = atol(result[1]); 
        }
    }
    sqlite3_free_table(result); 
    sqlite3_free(tmp);
    return ret;
}



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

#include <curl/curl.h>
#include <sqlite3.h>
#include <matio.h>
#include <archive.h>
#include <archive_entry.h>

#include "libufget.h"
#include "uf_internal.h"
#include "io/io.h"


void  uf_query_free(uf_query_t * query){
    if ( query == NULL) return ; 
    sqlite3_finalize(query->stmt); 
    free(query); 
    return; 
}

uf_query_t *uf_query_sql(uf_collection_t *col, const char *sql)
{
    uf_query_t *query; 
    int rc; 
    char *isql; 
    
    if ( col == NULL) return NULL; 
    
    query = (uf_query_t *) malloc(sizeof(uf_query_t) * (1)); 
    if ( query == NULL ) return NULL; 
    query->done = 0; 
    query->col = col; 

    isql = sqlite3_mprintf("SELECT * FROM matrices WHERE %s;", sql); 
    rc = sqlite3_prepare_v2(col->db, isql, -1, &(query->stmt), NULL); 
    sqlite3_free(isql); 
    if ( rc != SQLITE_OK ) {
        fprintf(stderr, "Failed to setup the sql query. Error: %s\n",  sqlite3_errmsg(col->db));
        sqlite3_finalize(query->stmt); 
        free(query); 
        return NULL; 
    }
    return query; 
}

uf_matrix_t * uf_query_next(uf_query_t *query) 
{
    uf_matrix_t *mat; 
    int rc; 

    if ( query == NULL ) return NULL; 
    if ( query->done ==1 ) return NULL; 



    rc = sqlite3_step(query->stmt); 

    if ( rc == SQLITE_ROW ) {
        /* Return Matrix  */
        mat = (uf_matrix_t *) malloc(sizeof(uf_matrix_t) * (1)); 
        if ( mat == NULL ) return NULL; 
        memset(mat, 0, sizeof(uf_matrix_t)); 
        mat->col = query->col; 
        uf_matrix_from_stmt(mat, query->stmt); 
        return mat; 
    } else {
        /* Error or done   */
        query -> done = 1; 
        return NULL; 
    }
    return NULL;
}





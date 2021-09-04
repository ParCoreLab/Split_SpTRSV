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


uf_query_t * uf_query_spd(uf_collection_t *col)
{
    uf_query_t * query; 

    query = uf_query_sql(col, "nrows == ncols AND isReal==1 AND numerical_symmetry == 1.0 AND posdef==1"); 
    if (!query) {
        fprintf(stderr, "Failed to query for SPD matrices.\n");
        return NULL; 
    }
    return query;  
}

uf_query_t * uf_query_square(uf_collection_t *col, uf_field_t field)
{
    uf_query_t * query; 
    char fieldquery[100]; 
    char sql[200]; 

    if (field == UF_REAL) {
        strncpy(fieldquery, "isReal == 1", 100); 
    } else if (field == UF_COMPLEX) {
        strncpy(fieldquery, "isComplex == 1", 100); 
    } else if ( field == UF_PATTERN ) {
        strncpy(fieldquery, "isBinary == 1", 100); 
    } else {
        return NULL; 
    }
    snprintf(sql, 200, "nrows == ncols AND %s", fieldquery); 
    query = uf_query_sql(col, sql); 
    if (!query) {
        fprintf(stderr, "Failed to query for square matrices (%s) .\n", fieldquery);
        return NULL; 
    }
    return query;  

}

uf_query_t * uf_query_rectangular(uf_collection_t *col, uf_field_t field)
{
    uf_query_t * query; 
    char fieldquery[100]; 
    char sql[200]; 

    if (field == UF_REAL) {
        strncpy(fieldquery, "isReal == 1", 100); 
    } else if (field == UF_COMPLEX) {
        strncpy(fieldquery, "isComplex == 1", 100); 
    } else if ( field == UF_PATTERN ) {
        strncpy(fieldquery, "isBinary == 1", 100); 
    } else {
        return NULL; 
    }
    snprintf(sql, 200, "nrows != ncols AND %s", fieldquery); 
    query = uf_query_sql(col, sql); 
    if (!query) {
        fprintf(stderr, "Failed to query for rectangular matrices (%s) .\n", fieldquery);
        return NULL; 
    }
    return query;  

}


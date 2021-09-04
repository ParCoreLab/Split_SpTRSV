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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

#include "libufget.h"


#define MAX(a,b) ((a)<(b)?(b):(a))
#define MIN(a,b) ((a)<(b)?(a):(b))
int main(int argc, char **argv)
{
    uf_collection_t *col; 
    uf_matrix_t *mat = NULL; 
    int i; 

    printf("libufget - download all\n");
    printf("==================\n");



    /* Init Collection  */
    col = uf_collection_init(); 
    printf("Database contains %d matrices.\n", uf_collection_num_matrices(col));

    /* Download all matrices. Attention: IDs start with 1 do be consistent with the website and matlab  */
    for (i = 0; i < uf_collection_num_matrices(col); i++) {
        mat = uf_collection_get_by_id(col, i+1); 
        if (!mat) {
            printf("Failed to get information about %d.\n", i+1);
            continue; 
        }
        uf_matrix_get(mat); 
        uf_matrix_free(mat); 
        
    }

    /* Clean up  */
    uf_collection_finalize(col); 

    return 0;
}

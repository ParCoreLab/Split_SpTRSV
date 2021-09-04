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
#include <string.h>

#include <curl/curl.h>
#include <sqlite3.h>
#include <matio.h>
#include <archive.h>
#include <archive_entry.h>


#include "libufget.h"
#include "uf_internal.h"

void print_usage() {
    char *baseurl;
    char *cache_dir;

    if ( getenv("UF_COLLECTION_BASEURL") != NULL) {
        baseurl = strdup(getenv("UF_COLLECTION_BASEURL"));
    } else {
        baseurl = strdup(UF_INDEX_URL);
    }

    if ( getenv("UF_COLLECTION_CACHE_DIR") != NULL) {
        cache_dir = strdup(getenv("UF_COLLECTION_CACHE_DIR"));
    } else {
        size_t len = strlen(getenv("HOME"))+100;
        cache_dir = (char *) malloc(sizeof(char) * (len));
        snprintf(cache_dir, len, "%s/.ufget/", getenv("HOME"));
    }


    printf("libufget command line tool\n");
    printf("--------------------------\n");
    printf("Version %d.%d.%d\n", LIBUFGET_MAJOR, LIBUFGET_MINOR, LIBUFGET_PATCH);
    printf(UFGET_COPYRIGHT_STRING "\n");
    printf("This is free software under the terms of the GPLv2 (or any later version).\n");
    printf("This program is distributed in the hope that it will be useful,\n");
    printf("\n");
    printf("but WITHOUT ANY WARRANTY; without even the implied warranty of\n");
    printf("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n");
    printf("GNU General Public License for more details.\n");
    printf("\n");
    printf("ufget-update <command>\n");
    printf("commands:\n");
    printf(" update               force updateing the local collection database.\n");
    printf(" download-all         download all matrices from the collection into the local cache\n");
    printf(" download group name  download a specific matrix into the local cache.\n");
    printf("\n");
    printf("current cache dir: %s\n", cache_dir);
    printf("current base url:  %s\n", baseurl);
    printf("\n");
    free(baseurl);
    free(cache_dir);
    return ;
}

int main(int argc, char **argv)
{
    uf_collection_t *col = NULL ;
    uf_matrix_t *mat = NULL;
    int i;

    if ( argc <= 1 ) {
        print_usage();
        return 0;
    }

    /* Force Updating the data base.   */
    if (strcmp(argv[1], "update") == 0 ) {
        printf("Force update of the database. This might take a while. \n");
        col = uf_collection_init1(UF_COLLECTION_VERBOSE| UF_COLLECTION_FORCE_REFRESH);
        if ( col == NULL ) {
            printf("Refresh failed.\n");
            return -1;
        }
        printf("Collection now contains %d matrices.               \n", uf_collection_num_matrices(col));
        uf_collection_finalize(col);
        return 0;
    }
    /* Download the whole database.  */
    else if ( strcmp(argv[1], "download-all") == 0 ){
        printf("Downloading the whole database. This might take a while.\n");
        col = uf_collection_init1(UF_COLLECTION_VERBOSE);
        for (i = 0; i < uf_collection_num_matrices(col); i++) {
            mat = uf_collection_get_by_id(col, i+1);
            if (!mat) {
                printf("Failed to get information about matrix ID = %d.\n", i+1);
                continue;
            }
            if(uf_matrix_get(mat)) {
                printf("--> Download of %s/%s failed.\n", mat->group_name, mat->name);
            }
            uf_matrix_free(mat);
        }
        uf_collection_finalize(col);
    }
    /* Download a specific matrix   */
    else if ( strcmp(argv[1], "download") == 0 ) {
        if (argc != 4 ) {
            printf("Invalid number of arguments.\n");
            print_usage();
            return -1;
        }
        col = uf_collection_init1(UF_COLLECTION_VERBOSE);
        if ( col == NULL ) {
            printf("Initializing collection failed.\n");
            return -1;
        }

        mat = uf_collection_get_by_name(col, argv[2], argv[3]);
        if (!mat) {
            printf("Cannot find matrix %s/%s.\n", argv[2], argv[3]);
            uf_collection_finalize(col);
            return -1;
        }
        if(uf_matrix_get(mat)) {
            printf("--> Download of %s/%s failed.\n", mat->group_name, mat->name);
        }
        uf_matrix_free(mat);

    } else {
        printf("Invalid argument.\n\n");
        print_usage();
        return -1;
    }
    return 0;
}


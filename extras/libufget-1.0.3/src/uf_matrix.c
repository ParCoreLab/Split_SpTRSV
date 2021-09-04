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
 * Copyright (C) Martin Koehler, 2015, 2017
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <stdint.h>
#include <curl/curl.h>
#include <sqlite3.h>
#include <matio.h>
#include <archive.h>
#include <archive_entry.h>

#include "libufget.h"
#include "uf_internal.h"
#include "io/io.h"


void uf_matrix_free(uf_matrix_t * matrix)
{
    if ( matrix == NULL) return; 
    free(matrix->localpath); 
    free(matrix); 
}

void uf_matrix_print(uf_matrix_t * matrix, int brief)
{
    printf("ID:                     : %d\n", matrix->id);
    printf("Group                   : %s\n",matrix->group_name);
    printf("Name                    : %s\n",matrix->name );
    printf("Rows                    : %" PRId64 "\n",matrix->rows );
    printf("Columns                 : %" PRId64 "\n",matrix->cols );
    printf("Non Zeros               : %" PRId64 "\n",matrix->nnz);
    if (!brief){
    printf("Zeros stored            : %" PRId64 "\n",matrix->nzero );
    printf("Pattern Symmetry        : %lg\n",matrix->pattern_symmetry );
    printf("Numerical Symmetry      : %lg\n",matrix->numerical_symmetry );
    printf("Binary                  : %d\n",matrix->isBinary );
    printf("Real                    : %d\n",matrix->isReal );
    printf("Complex                 : %d\n",matrix->isComplex );
    printf("nnzdiag                 : %" PRId64 "\n",matrix->nnzdiag );
    printf("posdef                  : %d\n",matrix->posdef );
    printf("amd_lnz                 : %" PRId64 "\n",matrix->amd_lnz );
    printf("amd_flops               : %" PRId64 "\n",matrix->amd_flops );
    printf("amd_vnz                 : %" PRId64 "\n",matrix->amd_vnz );
    printf("amd_rnz                 : %" PRId64 "\n",matrix->amd_rnz );
    printf("nblocks                 : %" PRId64 "\n",matrix->nblocks );
    printf("sprank                  : %" PRId64 "\n",matrix->sprank );
    printf("RBtype                  : %c%c%c\n",matrix->RBtype[0], matrix->RBtype[1], matrix->RBtype[2] );
    printf("cholcand                : %d\n",matrix->cholcand );
    printf("ncc                     : %" PRId64 "\n",matrix->ncc );
    printf("isND                    : %d\n",matrix->isND );
    printf("isGraph                 : %d\n",matrix->isGraph);
    printf("lowerbandwidth          : %" PRId64 "\n",matrix->lowerbandwidth );
    printf("upperbandwidth          : %" PRId64 "\n",matrix->upperbandwidth );
    printf("rcm_lowerbandwidth      : %" PRId64 "\n",matrix->rcm_lowerbandwidth );
    printf("rcm_upperbandwidth      : %" PRId64 "\n",matrix->rcm_upperbandwidth );
    printf("xmin                    : %lg + %lg i\n",matrix->xmin_real, matrix->xmin_imag );
    printf("xmax                    : %lg + %lg i\n",matrix->xmax_real, matrix->xmax_imag );
    if (matrix->localpath) 
    printf("localpath:              : %s\n", matrix->localpath);
    }
}

struct curl_prg_data {
    char *name; 
    time_t start; 
}; 

/*-----------------------------------------------------------------------------
 *  Curl Helper
 *-----------------------------------------------------------------------------*/
static int xferinfo(void *p, curl_off_t dltotal, curl_off_t dlnow, curl_off_t ultotal, curl_off_t ulnow)
{
    static time_t lastcall= 0; 
    struct curl_prg_data * dat = (struct curl_prg_data * ) p; 
    time_t curcall = time(NULL); 
    if ( lastcall < curcall ){
        long speed; 
        if (curcall == dat->start) {
            speed = 0; 
        } else {
            speed = dlnow / (curcall - dat->start); 
        }
        printf("Downloading Matrix (%s) (%06ldKiB / %06ldKiB) @ %8ld KiB/s    \r",(char*) dat->name, (long) dlnow/1024, (long) dltotal/1024, (long) speed/1024);
        fflush(stdout); 
        lastcall = curcall; 
    }
    return 0;
}

static int older_progress(void *p,  double dltotal, double dlnow, double ultotal, double ulnow) 
{
  return xferinfo(p, (curl_off_t)dltotal, (curl_off_t)dlnow, (curl_off_t)ultotal,(curl_off_t)ulnow);
}



/*-----------------------------------------------------------------------------
 *  Download a matrix 
 *-----------------------------------------------------------------------------*/
int  uf_matrix_get(uf_matrix_t * mat)
{
    char *url = NULL;
    char *localpath = NULL; 
    size_t len_base, len_group, len_name, local_len; 
    CURL *curl; 
    char errbuf[CURL_ERROR_SIZE];
    FILE *output; 
    struct curl_prg_data dat; 


    if ( mat->localpath != NULL && uf_file_exist(mat->localpath) ) {
        return 0; 
    }


    len_base = strlen(mat->col->baseurl); 
    len_group = strlen(mat->group_name); 
    len_name = strlen(mat->name); 
    

    len_base += len_group + len_name +100; 
    local_len = strlen(mat->col->cache_dir) + len_name + len_group + 100; 

    url = (char *) malloc(sizeof(char) * (len_base)); 
    localpath = (char *) malloc(sizeof(char) * local_len); 

    snprintf(url, len_base, "%s/MM/%s/%s.tar.gz", mat->col->baseurl, mat->group_name, mat->name); 
    snprintf(localpath, local_len, "%s/MM/%s/", mat->col->cache_dir, mat->group_name); 
    if (uf_file_mkdir(localpath, S_IRWXU)) {
        fprintf(stderr, "mkdir (%s) failed.\n", localpath);
        goto failed; 
    }
    snprintf(localpath, local_len, "%s/MM/%s/%s.tar.gz", mat->col->cache_dir, mat->group_name, mat->name); 

   
    DPRINTF(1, "Download matrix from %s to %s\n", url, localpath); 

    output = fopen(localpath, "w"); 
    if ( !output) {
        fprintf(stderr, "Opening %s for writing failed.\n", localpath);
        goto failed; 
    }

    curl = curl_easy_init(); 
    if ( !curl ) {
        fprintf(stderr, "Error initializing CURL.\n");
        goto failed; 
    }


    dat.name = mat->name; 
    dat.start = time(NULL); 
    curl_easy_setopt(curl, CURLOPT_URL, url); 
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, output);
    curl_easy_setopt(curl, CURLOPT_ERRORBUFFER, errbuf);
    curl_easy_setopt(curl, CURLOPT_FAILONERROR, 1); 

    if ( uf_get_verbose_level > 0 ) {
        curl_easy_setopt(curl, CURLOPT_PROGRESSFUNCTION, older_progress);
        curl_easy_setopt(curl, CURLOPT_PROGRESSDATA, &dat);
#if LIBCURL_VERSION_NUM >= 0x072000
        curl_easy_setopt(curl, CURLOPT_XFERINFOFUNCTION, xferinfo);
        curl_easy_setopt(curl, CURLOPT_XFERINFODATA, &dat);
#endif
        curl_easy_setopt(curl, CURLOPT_NOPROGRESS, 0); 
    }
    
    CURLcode res = curl_easy_perform(curl); 
    if ( res ) {
        curl_easy_cleanup(curl); 
        curl_global_cleanup(); 
        fclose(output); 
        remove(localpath); 
        fflush(stdout); 
        fflush(stderr); 
        fprintf(stderr, "\nFailed to get %s - %s\n", url, errbuf);
        fflush(stdout); 
        fflush(stderr); 
        goto failed; 
    }

    /* New line after download  */
    if (uf_get_verbose_level > 0 ) printf("\n");
    curl_easy_cleanup(curl); 
    curl_global_cleanup(); 
    fclose(output); 

    mat->localpath = localpath; 
    uf_matrix_set_localpath(mat, localpath); 
    free(url); 

    return 0; 

failed: 
    free(localpath); 
    free(url); 
    return -1; 
}

/*-----------------------------------------------------------------------------
 *  Update the local path in a matrix
 *-----------------------------------------------------------------------------*/
int uf_matrix_set_localpath(uf_matrix_t *mat, const char * localpath)
{
    sqlite3 *db; 
    int ret; 

    if (mat == NULL) return -1; 
    if (mat->col == NULL) return -1; 

    db = mat->col->db; 
    ret = sql_execute(db, "UPDATE matrices SET localpath='%q' WHERE id=%d AND group_name='%q' AND name = '%q' ;", localpath,mat->id, mat->group_name, mat->name);
    if ( ret != SQLITE_OK )
        return -1; 
    return 0; 
}


/*-----------------------------------------------------------------------------
 *  Matrix from Statement 
 *-----------------------------------------------------------------------------*/
void uf_matrix_from_stmt(uf_matrix_t *matrix, sqlite3_stmt *res) 
{
/* Extra information from DB   */
    matrix -> id = sqlite3_column_int(res, COL_ID); 
    strncpy(matrix->group_name, (char *) sqlite3_column_text(res, COL_GROUP), 1023); 
    strncpy(matrix->name, (char *) sqlite3_column_text(res, COL_NAME), 1023); 
    matrix->rows = sqlite3_column_int64(res, COL_NROWS); 
    matrix->cols = sqlite3_column_int64(res, COL_NCOLS); 
    matrix->nnz   = sqlite3_column_int64(res, COL_NNZ); 
    matrix->nzero = sqlite3_column_int64(res, COL_NZERO); 
    matrix->pattern_symmetry = sqlite3_column_double(res, COL_PATTERN_SYMMETRY); 
    matrix->numerical_symmetry = sqlite3_column_double(res, COL_NUMERICAL_SYMMETRY); 
    matrix->isBinary  = sqlite3_column_int(res, COL_ISBINARY); 
    matrix->isReal    = sqlite3_column_int(res, COL_ISREAL); 
    matrix->isComplex   = sqlite3_column_int(res, COL_ISCOMPLEX);
    matrix->nnzdiag = sqlite3_column_int64(res, COL_NNZDIAG); 
    matrix->posdef  = sqlite3_column_int(res, COL_POSDEF); 
    matrix->amd_lnz = sqlite3_column_int64(res, COL_AMD_LNZ); 
    matrix->amd_flops = sqlite3_column_int64(res, COL_AMD_FLOPS); 
    matrix->amd_vnz = sqlite3_column_int64(res, COL_AMD_VNZ); 
    matrix->amd_rnz = sqlite3_column_int64(res, COL_AMD_RNZ); 
    matrix->nblocks = sqlite3_column_int64(res, COL_NBLOCKS); 
    matrix->sprank  = sqlite3_column_int64(res, COL_SPRANK); 
    strncpy(matrix->RBtype, (char *) sqlite3_column_text(res, COL_RBTYPE), 4); 
    matrix->cholcand  = sqlite3_column_int(res, COL_CHOLCAND); 
    matrix->ncc      = sqlite3_column_int64(res, COL_NCC); 
    matrix->isND = sqlite3_column_int(res, COL_ISND); 
    matrix->isGraph = sqlite3_column_int(res, COL_ISGRAPH); 
    matrix->lowerbandwidth = sqlite3_column_int64(res, COL_LOWERBANDWIDTH); 
    matrix->upperbandwidth = sqlite3_column_int64(res, COL_UPPERBANDWIDTH); 
    matrix->rcm_lowerbandwidth = sqlite3_column_int64(res, COL_RCM_LOWERBANDWIDTH); 
    matrix->rcm_upperbandwidth = sqlite3_column_int64(res, COL_RCM_UPPERBANDWIDTH); 
    matrix->xmin_real = sqlite3_column_double(res, COL_XMIN_REAL); 
    matrix->xmin_imag = sqlite3_column_double(res, COL_XMIN_IMAG); 
    matrix->xmax_real = sqlite3_column_double(res, COL_XMAX_REAL); 
    matrix->xmax_imag = sqlite3_column_double(res, COL_XMAX_IMAG); 

    const char *lp = (const char*)sqlite3_column_text(res, COL_LOCALPATH); 
    if (lp) 
        matrix->localpath = strdup(lp); 
}


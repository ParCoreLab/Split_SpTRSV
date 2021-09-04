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

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <complex.h>

#include <curl/curl.h>
#include <sqlite3.h>
#include <matio.h>
#include <archive.h>
#include <archive_entry.h>

#include "libufget.h"
#include "uf_internal.h"
#include "io/io.h"

static int updatedb(sqlite3 *db, const char *baseurl);
static long get_int_from_mat(matvar_t * mat, size_t idx);
static double get_double_from_mat(matvar_t * mat, size_t idx);
static double complex get_complex_from_mat(matvar_t * mat, size_t idx);
static void get_rbtype_from_mat(char *buf, matvar_t *mat, int idx);
static int get_matrix_from_sql(uf_matrix_t *matrix, sqlite3 *db, const char *sql);

/* Global Verbose Level  */
int uf_get_verbose_level = 0;

/* Support MATIO 1.3 */
#if MATIO_MAJOR_VERSION == 1 && MATIO_MINOR_VERSION == 3
static  matvar_t * Mat_VarGetStructFieldByName(matvar_t * mat, const char *name, size_t index)
{
	return Mat_VarGetStructField(mat, (void *) name, BY_NAME, index);
}

#endif

/* Set Flag from environment  */
static void set_from_env(unsigned long *flags, const char *envname, unsigned long flag)
{
    if ( getenv(envname) != NULL) {
        int v = atoi(getenv(envname));
        if (v > 0) {
            *flags |= flag;
        }
    }
}

uf_collection_t * uf_collection_init()
{
    unsigned long flags = UF_COLLECTION_DEFAULT;
    set_from_env(&flags, "UF_COLLECTION_CLEANUP", UF_COLLECTION_CLEANUP);
    set_from_env(&flags, "UF_COLLECTION_VERBOSE", UF_COLLECTION_VERBOSE);
    return uf_collection_init1(flags);
}

uf_collection_t * uf_collection_init1(unsigned long flags)
{
    char * baseurl;
    char *cache_dir;
    uf_collection_t * col;

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

    col = uf_collection_init2(baseurl, cache_dir,  flags);
    free(baseurl);
    free(cache_dir);
    return col;
}

uf_collection_t * uf_collection_init2(const char *baseurl, const char *cache_dir, unsigned long flags)
{
    uf_collection_t *col;
    int rc;
    long timestamp, current_time;
    int update = 0;
    int db_version = 0;
    char *dbpath = NULL;

    col = (uf_collection_t *) malloc(sizeof(uf_collection_t) * (1));
    col->baseurl = strdup(baseurl);
    col->cache_dir = strdup(cache_dir);
    col->flags = flags;

    /* Set verbose level  */
    if ( flags & UF_COLLECTION_VERBOSE ) {
        uf_get_verbose_level = 1;
    }

    /* Basic Dir  */
    uf_file_mkdir(col->cache_dir,S_IRWXU);
    /* Init database  */
    dbpath = (char *) malloc(sizeof(char) * (strlen(col->cache_dir)+100));
    sprintf(dbpath, "%s/ufmatrices.db", col->cache_dir);
    rc = sqlite3_open(dbpath, &(col->db));

    if (rc != SQLITE_OK) {
        fprintf(stderr, "Cannot open database: %s\n", sqlite3_errmsg(col->db));
        goto failed;
    }

    /* Check the database structure  */
    if ( !sql_table_exists(col->db, "meta")){
        rc = sql_execute(col->db, "CREATE TABLE meta (key TEXT UNIQUE, value TEXT);");
        if ( rc != SQLITE_OK ) {
            fprintf(stderr, "Failed to create meta table.\n");
            goto failed;
        }
    }
    /* Pick META information from DB  */
    else {
        db_version = sql_get_long(col->db, "SELECT value FROM meta WHERE key='db_version';");
        if (db_version != UF_SQLITE_TABLE_VERSION ) {
            fprintf(stderr, "Database version mismatch. Force re-download.\n");
            sql_execute(col->db, "DROP TABLE IF EXISTS matrices;");
            update = 1;
        }
    }

    /* Create New Matrix table  */
    if ( !sql_table_exists(col->db, "matrices")){
        rc =sql_execute(col->db, "CREATE TABLE IF NOT EXISTS matrices( id INTEGER PRIMARY KEY, group_name TEXT, name TEXT, nrows INTEGER, ncols INTEGER, nnz INTEGER, nzeros INTEGER, "
                            "pattern_symmetry REAL, numerical_symmetry REAL, isBinary INTEGER, isReal INTEGER, isComplex INTEGER, nnzdiag INTEGER, posdef INTEGER, amd_lnz INTEGER, amd_flops INTEGER,"
                            "amd_vnz INTEGER, amd_rnz INTEGER, nblocks INTEGER, sprank INTEGER, RBtype TEXT, cholcand INTEGER, ncc INTEGER, isND INTEGER, isGraph INTEGER, "
                            "lowerbandwidth INTEGER, upperbandwidth INTEGER, rcm_lowerbandwidth INTEGER, rcm_upperbandwidth INTEGER, xmin_real REAL, xmin_imag REAL, xmax_real REAL, xmax_imag REAL, localpath TEXT DEFAULT NULL);");
        if ( rc != SQLITE_OK ) {
            fprintf(stderr, "Failed to create matrices table.\n");
            goto failed;
        }

    }

    /* Get timestamp from meta   */
    current_time = (long) time(NULL);
    timestamp = sql_get_long(col->db, "SELECT value FROM meta WHERE key='timestamp';");

    if ( update || timestamp < current_time - 14*24*60*60 || flags & UF_COLLECTION_FORCE_REFRESH ) {
        update = 1;
    } else {
        update = 0;
    }

    // update = 1;
    /* Update Database */
    if ( update ) {
       rc = updatedb(col->db, baseurl);
       if ( rc != 0 ) {
            fprintf(stderr, "Failed to update the database.\n");
            goto failed;
       }
       rc = sql_execute(col->db, "DELETE FROM meta WHERE key='timestamp';");
       rc += sql_execute(col->db, "INSERT INTO meta VALUES('timestamp', '%ld');", current_time);
       if ( rc != SQLITE_OK ) {
           fprintf(stderr, "Failed to update timestamp.\n");
       }
       rc = sql_execute(col->db, "DELETE FROM meta WHERE key='db_version';");
       rc += sql_execute(col->db, "INSERT INTO meta VALUES('db_version', '%ld');", UF_SQLITE_TABLE_VERSION);
       if ( rc != SQLITE_OK ) {
           fprintf(stderr, "Failed to update timestamp.\n");
       }

    }

    free(dbpath);

    return col;

failed:
    sqlite3_close(col->db);
    free(col->baseurl);
    free(col->cache_dir);
    free(col);
    free(dbpath);
    return NULL;
}

void uf_collection_finalize(uf_collection_t * collection)
{
    if (collection == NULL) return;
    if (collection->flags & UF_COLLECTION_CLEANUP ) {
        uf_collection_cleanup(collection);
    }
    sqlite3_close(collection->db);
    free(collection->baseurl);
    free(collection->cache_dir);
    free(collection);
}

int uf_collection_num_matrices(uf_collection_t *col)
{
    int matrices;
    matrices = sql_get_long(col->db, "SELECT COUNT(*) FROM matrices;");
    return matrices;
}


uf_matrix_t * uf_collection_get_by_id(uf_collection_t * collection, int id)
{
    uf_matrix_t * mat;
    int rc;
    char * sql;

    mat = (uf_matrix_t *) malloc(sizeof(uf_matrix_t) * (1));
    if ( mat == NULL ) return NULL;
    memset(mat, 0, sizeof(uf_matrix_t));

    mat->col = collection;
    sql = sqlite3_mprintf("SELECT * FROM matrices WHERE id = %d;", id);
    rc = get_matrix_from_sql(mat, collection->db, sql);
    sqlite3_free(sql);
    if ( rc != 0 ) {
        fprintf(stderr, "Failed to get matrix with id = %d\n", id);
        uf_matrix_free(mat);
        return NULL;
    }
    return mat;
}


uf_matrix_t * uf_collection_get_by_name(uf_collection_t * collection, const char * group_name, const char *name)
{
    uf_matrix_t * mat;
    int rc;
    char * sql;

    mat = (uf_matrix_t *) malloc(sizeof(uf_matrix_t) * (1));
    if ( mat == NULL ) return NULL;
    memset(mat, 0, sizeof(uf_matrix_t));

    mat->col = collection;
    sql = sqlite3_mprintf("SELECT * FROM matrices WHERE UPPER(group_name)=UPPER('%q') AND UPPER(name)=UPPER('%q');", group_name, name);
    rc = get_matrix_from_sql(mat, collection->db, sql);
    sqlite3_free(sql);
    if ( rc != 0 ) {
        fprintf(stderr, "Failed to get matrix with group_name = %s and name = %s\n", group_name, name);
        uf_matrix_free(mat);
        return NULL;
    }
    return mat;

}


/*-----------------------------------------------------------------------------
 *  Cleanup
 *-----------------------------------------------------------------------------*/
void uf_collection_cleanup( uf_collection_t * col) {
    uf_query_t * query;
    uf_matrix_t *mat;
    if (col == NULL ) return;
    if (col->db == NULL ) return;

    query =  uf_query_sql(col, "localpath != \"\"");
    while ( (mat=uf_query_next(query)) != NULL ) {
        // printf("To Clean: %s\n", mat->name);
        remove(mat->localpath);
        sql_execute(col->db, "UPDATE matrices SET localpath=\"\" WHERE id=%d;", mat->id);
        uf_matrix_free(mat);
    }
    uf_query_free(query);
    return;
}


/*-----------------------------------------------------------------------------
 *  Curl Helper
 *-----------------------------------------------------------------------------*/
static int xferinfo(void *p, curl_off_t dltotal, curl_off_t dlnow, curl_off_t ultotal, curl_off_t ulnow)
{

    char *s = (char *)p;
    char spin = *s;
    switch(spin){
        case '|':
            spin = '/';
            break;
        case '/':
            spin = '-';
            break;
        case '-':
            spin = '\\';
            break;
        case '\\':
            spin = '|';
            break;
        default:
            spin = '|';
    }
    *s = spin;
    printf("Downloading Database: %c\r", spin);
    return 0;
}

static int older_progress(void *p,  double dltotal, double dlnow, double ultotal, double ulnow)
{
  return xferinfo(p, (curl_off_t)dltotal, (curl_off_t)dlnow, (curl_off_t)ultotal,(curl_off_t)ulnow);
}

/*-----------------------------------------------------------------------------
 *  Update DB
 *-----------------------------------------------------------------------------*/
static int updatedb(sqlite3 *db, const char *baseurl)
{
    CURL *curl;
    char *url; ;
    char *output_name = "/tmp/curl_matio_temp";
    FILE *output;
    char spin = '|';
    char errbuf[CURL_ERROR_SIZE];
    char *sErrMsg;

    DPRINTF(1,"Update Database\n");

    output = fopen(output_name, "w");
    if (!output) {
        fprintf(stderr, "Failed to open temporary file %s\n", output_name);
        return -1;
    }

    curl = curl_easy_init();
    if ( !curl ) {
        fprintf(stderr, "Error initializing CURL.\n");
        return -1;
    }

    url = (char *) malloc(sizeof(char) * (strlen(baseurl)+strlen(UF_INDEX_FILE)+10));
    sprintf(url,"%s/%s", baseurl, UF_INDEX_FILE);
    DPRINTF(1,"Database URL: %s\n", url);

    curl_easy_setopt(curl, CURLOPT_URL, url);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, output);
    curl_easy_setopt(curl, CURLOPT_ERRORBUFFER, errbuf);
    curl_easy_setopt(curl, CURLOPT_FAILONERROR, 1);
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION , 1);
    if (uf_get_verbose_level > 0 ) {
        curl_easy_setopt(curl, CURLOPT_PROGRESSFUNCTION, older_progress);
        curl_easy_setopt(curl, CURLOPT_PROGRESSDATA, &spin);
#if LIBCURL_VERSION_NUM >= 0x072000
        curl_easy_setopt(curl, CURLOPT_XFERINFOFUNCTION, xferinfo);
        curl_easy_setopt(curl, CURLOPT_XFERINFODATA, &spin);
#endif
        curl_easy_setopt(curl, CURLOPT_NOPROGRESS, 0);
    }

    CURLcode res = curl_easy_perform(curl);
    if ( res ) {
        fprintf(stderr, "Failed to get %s - %s\n", url, errbuf);
        curl_easy_cleanup(curl);
        curl_global_cleanup();
        fclose(output);
        remove (output_name);
        free(url);

        return -1;
    }
    curl_easy_cleanup(curl);
    curl_global_cleanup();
    fclose(output);
    free(url);


    /*-----------------------------------------------------------------------------
     *  Parse DB
     *-----------------------------------------------------------------------------*/
    int rc;

    mat_t * matfile = Mat_Open(output_name, MAT_ACC_RDONLY);
    if (! matfile) {
        fprintf(stderr, "Failed to open temporary database.\n");
        return -1;
    }
    matvar_t *ufindex = Mat_VarRead(matfile, UF_CONTAINER_STRUCT);
    if ( ! ufindex ) {
        fprintf(stderr, "Failed to pick the container structure \"" UF_CONTAINER_STRUCT "\" from the database.\n");
        Mat_Close(matfile);
        return -1;
    }
    matvar_t *revision = Mat_VarGetStructFieldByName(ufindex, "LastRevisionDate", 0);
    DPRINTF(1,"\rDownloaded Database with Last Revision Date: %s\n", (char *) revision->data);

    rc = sql_execute(db, "DELETE FROM meta WHERE key = 'lastrev'; ");
    rc |= sql_execute(db, "INSERT INTO meta VALUES ('lastrev', '%s');", (char *) revision->data);
    if ( rc != SQLITE_OK ) {
        fprintf(stderr, "Failed to update revivision data.\n");
        goto failed;
    }

    matvar_t *groups = Mat_VarGetStructFieldByName(ufindex, "Group", 0);
    matvar_t *names  = Mat_VarGetStructFieldByName(ufindex, "Name", 0);
    matvar_t *nrows  = Mat_VarGetStructFieldByName(ufindex, "nrows", 0);
    matvar_t *ncols  = Mat_VarGetStructFieldByName(ufindex, "ncols", 0);
    matvar_t *nnz  = Mat_VarGetStructFieldByName(ufindex, "nnz", 0);
    matvar_t *nzero  = Mat_VarGetStructFieldByName(ufindex, "nzero", 0);
    matvar_t *pattern_symmetry  = Mat_VarGetStructFieldByName(ufindex, "pattern_symmetry", 0);
    matvar_t *numerical_symmetry  = Mat_VarGetStructFieldByName(ufindex, "numerical_symmetry", 0);
    matvar_t *isBinary  = Mat_VarGetStructFieldByName(ufindex, "isBinary", 0);
    matvar_t *isReal  = Mat_VarGetStructFieldByName(ufindex, "isReal", 0);
    matvar_t *nnzdiag  = Mat_VarGetStructFieldByName(ufindex, "nnzdiag", 0);
    matvar_t *posdef  = Mat_VarGetStructFieldByName(ufindex, "posdef", 0);
    matvar_t *amd_lnz  = Mat_VarGetStructFieldByName(ufindex, "amd_lnz", 0);
    matvar_t *amd_flops  = Mat_VarGetStructFieldByName(ufindex, "amd_flops", 0);
    matvar_t *amd_vnz  = Mat_VarGetStructFieldByName(ufindex, "amd_vnz", 0);
    matvar_t *amd_rnz  = Mat_VarGetStructFieldByName(ufindex, "amd_rnz", 0);
    matvar_t *nblocks  = Mat_VarGetStructFieldByName(ufindex, "nblocks", 0);
    matvar_t *sprank  = Mat_VarGetStructFieldByName(ufindex, "sprank", 0);
    matvar_t *RBtype  = Mat_VarGetStructFieldByName(ufindex, "RBtype", 0);
    matvar_t *cholcand  = Mat_VarGetStructFieldByName(ufindex, "cholcand", 0);
    matvar_t *ncc  = Mat_VarGetStructFieldByName(ufindex, "ncc", 0);
    matvar_t *isND  = Mat_VarGetStructFieldByName(ufindex, "isND", 0);
    matvar_t *isGraph  = Mat_VarGetStructFieldByName(ufindex, "isGraph", 0);
    matvar_t *lowerbandwidth  = Mat_VarGetStructFieldByName(ufindex, "lowerbandwidth", 0);
    matvar_t *upperbandwidth  = Mat_VarGetStructFieldByName(ufindex, "upperbandwidth", 0);
    matvar_t *rcm_lowerbandwidth  = Mat_VarGetStructFieldByName(ufindex, "rcm_lowerbandwidth", 0);
    matvar_t *rcm_upperbandwidth  = Mat_VarGetStructFieldByName(ufindex, "rcm_upperbandwidth", 0);
    matvar_t *xmin  = Mat_VarGetStructFieldByName(ufindex, "xmin", 0);
    matvar_t *xmax  = Mat_VarGetStructFieldByName(ufindex, "xmax", 0);

    int nmats = 0, i ;
    if ( groups->rank == 1 ) {
        nmats = groups->dims[0];
    } else if (groups->rank == 2 ) {
        nmats = MAX(groups->dims[0], groups->dims[1]);
    }
    DPRINTF(1,"Found %d matrices.\n", nmats);

    /* rc = sql_execute(db, "DELETE FROM matrices; ");
    if ( rc != SQLITE_OK) {
        fprintf(stderr, "Failed to remove old entries.\n");
        goto failed;
    } */
    sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &sErrMsg);
    for (i = 0; i < nmats; i++) {
        char rbstr[4];
        matvar_t *group = Mat_VarGetCell(groups,i);
        matvar_t *name  = Mat_VarGetCell(names,i);
        get_rbtype_from_mat(rbstr, RBtype, i);
        long isComplex = get_int_from_mat(isReal,i) + get_int_from_mat(isBinary,i);
        long isReal;
        double dbl_numerical_symmetry, dbl_pattern_symmetry;
        // long mat_exist = 0;


        if ( isComplex == 0 ) {
            isComplex = 1;
        } else {
            isComplex = 0;
        }
        if ( isComplex == 0 && get_int_from_mat(isBinary, i) == 0) {
            isReal = 1;
        } else {
            isReal = 0;
        }

        dbl_pattern_symmetry = get_double_from_mat(pattern_symmetry,i);
        dbl_numerical_symmetry = get_double_from_mat(numerical_symmetry,i);

        if ( i % 25 == 0 ) { DPRINTF(1,"Updating local database: (%5d / %5d)\r", i+1, nmats); fflush(stdout); }
        double complex vxmin = get_complex_from_mat(xmin,i);
        double complex vxmax = get_complex_from_mat(xmax,i);
        rc = sql_execute(db,"INSERT OR REPLACE  INTO matrices VALUES (%d, '%s', '%s'"
                ",%ld, %ld, %ld, %ld"
                ",%lg, %lg, %ld, %ld, %ld, %ld, %ld, %ld, %ld"
                ",%ld, %ld, %ld, %ld, '%s', %ld, %ld, %ld, %ld"
                ",%ld, %ld, %ld, %ld, %lg, %lg, %lg, %lg, NULL);"
                , i+1, (char*) group->data, (char *) name->data
                , get_int_from_mat(nrows,i), get_int_from_mat(ncols, i), get_int_from_mat(nnz, i), get_int_from_mat(nzero,i)
                , dbl_pattern_symmetry, dbl_numerical_symmetry , get_int_from_mat(isBinary,i), isReal, isComplex, get_int_from_mat(nnzdiag,i), get_int_from_mat(posdef,i), get_int_from_mat(amd_lnz,i), get_int_from_mat(amd_flops,i)
                , get_int_from_mat(amd_vnz,i), get_int_from_mat(amd_rnz,i), get_int_from_mat(nblocks, i), get_int_from_mat(sprank,i), rbstr, get_int_from_mat(cholcand,i), get_int_from_mat(ncc,i), get_int_from_mat(isND, i), get_int_from_mat(isGraph, i)
                ,get_int_from_mat(lowerbandwidth, i), get_int_from_mat(upperbandwidth, i), get_int_from_mat(rcm_lowerbandwidth,i), get_int_from_mat(rcm_upperbandwidth,i), creal(vxmin), cimag(vxmin), creal(vxmax), cimag(vxmax));

        if ( rc != SQLITE_OK ){
            fprintf(stderr, "Error inserting matrix %s/%s\n", (char*) group->data, (char *) name->data);
            return -1;
        }


    }
    sqlite3_exec(db, "END TRANSACTION", NULL, NULL, &sErrMsg);

    Mat_Close(matfile);
    remove(output_name);
    return 0;

failed:
    Mat_Close(matfile);
    remove(output_name);
    return -1;
}




/*-----------------------------------------------------------------------------
 *  Get long from matlab
 *-----------------------------------------------------------------------------*/
static long get_int_from_mat(matvar_t * mat, size_t idx) {
    switch( mat -> data_type) {
        case MAT_T_INT8:
            {
                int8_t *ptr = mat->data;
                return ptr[idx];
            }
        case MAT_T_UINT8:
            {
                uint8_t *ptr = mat->data;
                return ptr[idx];
            }
        case MAT_T_INT16:
            {
                int16_t *ptr = mat->data;
                return ptr[idx];
            }
        case MAT_T_UINT16:
            {
                uint16_t *ptr = mat->data;
                return ptr[idx];
            }
        case MAT_T_INT32:
            {
                int32_t *ptr = mat->data;
                return ptr[idx];
            }
        case MAT_T_UINT32:
            {
                uint32_t *ptr = mat->data;
                return ptr[idx];
            }
        case MAT_T_INT64:
            {
                int64_t *ptr = mat->data;
                return ptr[idx];
            }
        case MAT_T_UINT64:
            {
                uint64_t *ptr = mat->data;
                return ptr[idx];
            }
        case MAT_T_SINGLE:
            {
                float * ptr = mat->data;
                return nearbyint(ptr[idx]);
            }
        case MAT_T_DOUBLE:
            {
                double * ptr = mat->data;
                return nearbyint(ptr[idx]);
            }
        default:
            return 0;

    }
    return 0;
}

/*-----------------------------------------------------------------------------
 *  Get double from matlab
 *-----------------------------------------------------------------------------*/
static double get_double_from_mat(matvar_t * mat, size_t idx) {
    switch( mat -> data_type) {
        case MAT_T_INT8:
            {
                int8_t *ptr = mat->data;
                return ptr[idx];
            }
        case MAT_T_UINT8:
            {
                uint8_t *ptr = mat->data;
                return ptr[idx];
            }
        case MAT_T_INT16:
            {
                int16_t *ptr = mat->data;
                return ptr[idx];
            }
        case MAT_T_UINT16:
            {
                uint16_t *ptr = mat->data;
                return ptr[idx];
            }
        case MAT_T_INT32:
            {
                int32_t *ptr = mat->data;
                return ptr[idx];
            }
        case MAT_T_UINT32:
            {
                uint32_t *ptr = mat->data;
                return ptr[idx];
            }
        case MAT_T_INT64:
            {
                int64_t *ptr = mat->data;
                return ptr[idx];
            }
        case MAT_T_UINT64:
            {
                uint64_t *ptr = mat->data;
                return ptr[idx];
            }
        case MAT_T_SINGLE:
            {
                float * ptr = mat->data;
                return ptr[idx];
            }
        case MAT_T_DOUBLE:
            {
                double * ptr = mat->data;
                return ptr[idx];

            }
        default:
            return 0;

    }
    return 0;
}

#if MATIO_MAJOR_VERSION == 1 && MATIO_MINOR_VERSION < 4
typedef struct ComplexSplit mat_complex_split_t;
#endif

static double complex get_complex_from_mat(matvar_t * mat, size_t idx)
{

    mat_complex_split_t *cpx_data;
    if ( ! mat->isComplex ) {
        return get_double_from_mat(mat, idx);
    }
    cpx_data = (mat_complex_split_t *) mat->data;
    switch(mat->data_type) {
        case MAT_T_SINGLE:
            {
                float *ptr_real = cpx_data->Re;
                float *ptr_imag = cpx_data->Im;

                return ptr_real[idx] + ptr_imag[idx]*I;
            }
        case MAT_T_DOUBLE:
            {
            double *ptr_real = cpx_data->Re;
            double *ptr_imag = cpx_data->Im;

            return ptr_real[idx] + ptr_imag[idx]*I;
            }
        default:
            return 0.0;
    }
    return 0.0;

}

/*-----------------------------------------------------------------------------
 *  Get the RB type from matlab
 *-----------------------------------------------------------------------------*/
static void get_rbtype_from_mat(char *buf, matvar_t *mat, int idx)
{
    int r = mat->rank;
    if ( r == 1) {
        strncpy(buf, mat->data, 3);
        buf[3] = 0;
        return;
    } else if (r==2){
        int ld = mat->dims[0];
        char *strarry = (char*) mat->data;
        buf[0] = strarry[idx];
        buf[1] = strarry[idx+ld];
        buf[2] = strarry[idx+2*ld];
        buf[3] = 0;
    } else {
        buf[0] = 0;
        return;
    }
}

/*-----------------------------------------------------------------------------
 *  Get Matrix infor from SQL
 *-----------------------------------------------------------------------------*/
static int get_matrix_from_sql(uf_matrix_t *matrix, sqlite3 *db, const char *sql)
{
    int rc;
    sqlite3_stmt *res;

    rc = sqlite3_prepare_v2(db, sql, -1, &res, NULL);
    if ( rc != SQLITE_OK ) {
         fprintf(stderr, "Failed to fetch data: %s\n", sqlite3_errmsg(db));
         return -1;
    }
    rc= sqlite3_step(res);
    if ( rc != SQLITE_ROW) {
        fprintf(stderr, "No entry matching '%s' found.\n", sql);
        sqlite3_finalize(res);
        return -1;
    }

    uf_matrix_from_stmt(matrix, res);

    sqlite3_finalize(res);
    return 0;
}


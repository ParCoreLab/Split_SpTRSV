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
#include <complex.h>
#include <ctype.h>

#include <curl/curl.h>
#include <sqlite3.h>
#include <matio.h>
#include <archive.h>
#include <archive_entry.h>

#include "libufget.h"
#include "uf_internal.h"
#include "io/io.h"

struct _libarchive_read_t {
    size_t size; 
    int64_t offset; 
    const void *buffer; 
    size_t pos;
    struct archive * a;
};

static int _libarchive_close(void **data)
{
    struct _libarchive_read_t * idata = (struct _libarchive_read_t *) *data;
   archive_read_close(idata->a);
#if ARCHIVE_VERSION_NUMBER < 3000000
   archive_read_finish(idata->a); 
#else
   archive_read_free(idata->a);
#endif

    free(*data);
    *data = NULL;
    return 0;
}


static size_t _libarchive_read(void *_data, void *buf, size_t len)
{
    struct _libarchive_read_t *data = (struct _libarchive_read_t*) _data;
    int r;
    size_t remain;
    const char *cbuffer;
    if ( data->pos == 0  && data->offset == 0) {
        r = archive_read_data_block(data->a, &(data->buffer), &(data->size), &(data->offset));
    } else if (data->pos == data->size) {
        r = archive_read_data_block(data->a, &(data->buffer), &(data->size), &(data->offset));
        data->pos = 0;
    } else {
        r = ARCHIVE_OK;
    }
    // printf("size: %ld \t pos: %ld \t offset: %ld\n", data->size, data->pos, data->offset);
    if (r!=ARCHIVE_OK) return 0;
    cbuffer = data->buffer;
    remain = data->size - data->pos;
    if (len >= remain) 
        len = remain;
    memcpy(buf, &cbuffer[data->pos], len);
    data->pos+=len;
    return len;
}

struct _compressed_io_handler_t libarchive_read_handler = {
    .extension = "",
    .type = UF_IO_FILE_HANDLE,
    .magic = "",
    .open = NULL,
    .close = _libarchive_close,
    .read  = _libarchive_read,
    .write = NULL
};

static uf_io_file_t* _archive_open(struct archive *a )
{
	uf_io_file_t  *f;
    struct _libarchive_read_t *data;

	f = (uf_io_file_t *) malloc ( sizeof(uf_io_file_t));
	if ( f == NULL ) {
		fprintf(stderr, "Allocation of the mess_file object failed.\n");
		return NULL;
	}
    data =  (struct _libarchive_read_t *) malloc(sizeof(struct _libarchive_read_t) * (1));
    if ( data == NULL) {
        fprintf(stderr, "Failed to allocate internal data\n");
        return NULL;
    }
    data->size = 0;
    data->offset =  0;
    data->pos = 0;
    data->a = a;

	/*-----------------------------------------------------------------------------
	 *  Open the file
	 *-----------------------------------------------------------------------------*/
	f->pos = 0;
	f->avail = 0;
	f->data  = data;
    f->mode = UF_IO_FILE_READ;
	f->eof = 0; 
    f->handler = &libarchive_read_handler;
	return f;
}

static uf_io_file_t * uf_matrix_open(uf_matrix_t *mat) 
{
    char *internal_filename; 
    size_t len; 
    struct archive *a;
    struct archive_entry         *entry;
    int r; 

    if ( mat == NULL) return NULL; 

    /* Get if not available  */
    if ( mat->localpath == NULL) {
        if ( uf_matrix_get(mat)) {
            return NULL; 
        }
    }

    len = 2*strlen(mat->name) + 15; 
    internal_filename = (char *) malloc(sizeof(char) * (len)); 
    snprintf(internal_filename, len, "%s/%s.mtx", mat->name, mat->name); 

    a = archive_read_new();
    if (!a) {
        fprintf(stderr, "Init libarchive failed.\n");
        free(internal_filename); 
        return NULL; 
    }

   archive_read_support_format_all(a);
#if ARCHIVE_VERSION_NUMBER < 3000000
   archive_read_support_compression_all(a);
#else 
   archive_read_support_filter_all(a); 
#endif
#define ARCHIVE_BUFFER_SIZE 1024*10240
   if ((r = archive_read_open_filename(a, mat->localpath, ARCHIVE_BUFFER_SIZE))){
       fprintf(stderr, "Failed to open - Error: %s\n", archive_error_string(a));
#if ARCHIVE_VERSION_NUMBER < 3000000
       archive_read_finish(a); 
#else
       archive_read_free(a);
#endif
       free(internal_filename);
       return NULL;
   }
   
   while (1 ) {
        r = archive_read_next_header(a, &entry); 
        if ( r == ARCHIVE_EOF) 
            goto failed;
        if ( r != ARCHIVE_OK) {
            fprintf(stderr, "ERROR: %s\n", archive_error_string(a));
            goto failed;
        }

        // s = archive_entry_size(entry); 
        if ( strcmp(archive_entry_pathname(entry), internal_filename) == 0) {
            // fprintf(stderr, "Found %s in file.\n", internal_filename);
            break;
        } else {
            archive_read_data_skip(a);
            continue;
        }
   }
   free(internal_filename);
   uf_io_file_t * file = _archive_open(a);
   return file;

failed:
   archive_read_close(a);
#if ARCHIVE_VERSION_NUMBER < 3000000
   archive_read_finish(a); 
#else
   archive_read_free(a);
#endif
   return NULL;
}



/*-----------------------------------------------------------------------------
 *  MTX definitions 
 *-----------------------------------------------------------------------------*/
#define LINE_LENGTH 1025
#define TOKEN_LENGTH 30

/* Sets the last entry of an array to value */
#define 	SET_LAST(array, len, value) 	(array)[(len)-1] = (value)
/* Definitions for the Matrix Market header */
#define MM_BANNER		"%%MatrixMarket"
#define MM_STR_MATRIX		"matrix"
#define MM_STR_COORD		"coordinate"
#define MM_STR_ARR		"array"
#define MM_STR_REAL		"real"
#define MM_STR_INT		"integer"
#define MM_STR_CPX		"complex"
#define MM_STR_PAT		"pattern"
#define MM_STR_GE		"general"
#define MM_STR_SYM		"symmetric"
#define MM_STR_SKSYM		"skew-symmetric"
#define MM_STR_HER		"hermitian"

static void _lowercase (char *str)
{
	char *p;
	if ( str == NULL) return;
	for (p=str; *p!='\0'; *p=tolower(*p),p++);
}


/*-----------------------------------------------------------------------------
 *  Read and parse the header of an mtx file.
 *-----------------------------------------------------------------------------*/
static mtx_info_t * uf_matrix_read_info(uf_io_file_t * file)
{
	char line[LINE_LENGTH];
	char mmbanner[TOKEN_LENGTH];
	char mtx[TOKEN_LENGTH];
	char type[TOKEN_LENGTH];
	char datatype[TOKEN_LENGTH];
	char shape[TOKEN_LENGTH];
	int num_items_read;
    char * cret;
    mtx_info_t *matrix_info;

	int64_t N,M,NZ;

    matrix_info = (mtx_info_t *) malloc(sizeof(mtx_info_t) * (1));
	if ( matrix_info == NULL ) {
		fprintf(stderr, "Failed to allocate matrix_info.\n");
		return NULL; 
	}

	cret = uf_io_gets(line, LINE_LENGTH,file);
	if ( !cret ) {
		fprintf(stderr, "Error reading the first line of the file.\n");
        free(matrix_info);
		return NULL;
	}

	// %%MatrixMarket matrix type datatype shape
	if (sscanf(line, "%s %s %s %s %s", mmbanner, mtx, type, datatype, shape) != 5){
        fprintf(stderr, "Failed to parse the header line: %s\n", line);
        free(matrix_info);
        return NULL;
	}

	// prevent memory overflow
	SET_LAST(mmbanner, TOKEN_LENGTH, '\0');
	SET_LAST(mtx, TOKEN_LENGTH, '\0');
	SET_LAST(type, TOKEN_LENGTH, '\0');
	SET_LAST(datatype, TOKEN_LENGTH, '\0');
	SET_LAST(shape, TOKEN_LENGTH, '\0');

	_lowercase(mtx);
	_lowercase(type);
	_lowercase(datatype);
	_lowercase(shape);

	/* check if the matrix market banner exists */
	if (strncmp(mmbanner, MM_BANNER, strlen(MM_BANNER)) != 0){
		fprintf(stderr, "wrong header information: %s\n", mmbanner);
        free(matrix_info);
		return NULL;
	}

	/* parse the rest of the 1st line*/
	if (strncmp(mtx, MM_STR_MATRIX, strlen(MM_STR_MATRIX)) != 0){
		fprintf(stderr, "wrong header information: %s\n", mtx);
        free(matrix_info);
		return NULL;
	}

	if (strncmp(type, MM_STR_ARR,strlen(MM_STR_ARR)) == 0){
		matrix_info->store_type = MTX_DENSE;
	} else if (strncmp(type, MM_STR_COORD, strlen(MM_STR_COORD)) == 0){
		matrix_info->store_type = MTX_COORD;
	} else {
	    fprintf(stderr, "wrong storage format: %s\n", type );
        free(matrix_info);
		return NULL;
 	}

	if (strncmp(datatype, MM_STR_REAL,strlen(MM_STR_REAL)) == 0) {
		matrix_info->data_type = UF_REAL;
	} else if (strncmp(datatype, MM_STR_INT, strlen(MM_STR_INT)) == 0){
		matrix_info->data_type = UF_REAL;
	} else if (strncmp(datatype, MM_STR_CPX, strlen(MM_STR_CPX)) == 0){
		matrix_info->data_type = UF_COMPLEX;
	}
	else if (strncmp(datatype, MM_STR_PAT, strlen(MM_STR_PAT)) == 0) {
        // printf("Read pattern\n");
		matrix_info->data_type = UF_PATTERN;
	} else {
        fprintf(stderr,"wrong datatype: %s\n", datatype);
        free(matrix_info);
		return NULL;
 	}

	if (strncmp(shape, MM_STR_GE,strlen(MM_STR_GE)) == 0){
		matrix_info->symmetry = MTX_GENERAL;
	} else if (strncmp(shape, MM_STR_SYM, strlen(MM_STR_SYM)) == 0){
		matrix_info->symmetry = MTX_SYMMETRIC;
	} else if (strncmp(shape, MM_STR_SKSYM, strlen(MM_STR_SKSYM)) == 0) {
		matrix_info->symmetry = MTX_SKEWSYMMETRIC;
	} else if (strncmp(shape, MM_STR_HER, strlen(MM_STR_HER)) == 0) {
		matrix_info->symmetry = MTX_HERMITIAN;
	}  else {
        fprintf(stderr, "Invaild symmetry information: %s \n", shape);
        free(matrix_info);
		return NULL;
 	}

	// read until matrix size is read.
	do{ 
		cret = uf_io_gets(line, LINE_LENGTH, file);
		if ( !cret  ) {
			fprintf(stderr, "error while reading a line from file");
            free(matrix_info);
			return NULL;
		}
	} while (line[0] == '%');


	if (MTX_IS_COORD(matrix_info)){
		do {
			long _n, _m, _nz; 
			num_items_read = sscanf(line, "%ld %ld %ld", &_n, &_m, &_nz);
			if (num_items_read != 3) {
				cret = uf_io_gets(line, LINE_LENGTH, file);
				if ( !cret ) {
					 fprintf(stderr, "error while reading a line from file\n");
                     free(matrix_info);
					 return  NULL;
				}
			}
			N = _n; 
			M = _m; 
			NZ = _nz; 
		} while (num_items_read != 3);
	}else{
		NZ = 0;
		do {
			long _n, _m; 
			num_items_read = sscanf(line,"%ld %ld", &_n, &_m);
			if (num_items_read != 2) {
				cret = uf_io_gets(line, LINE_LENGTH, file);
				if ( !cret   ) {
					fprintf(stderr, "error while reading a line from file\n");
                    free(matrix_info);
					return NULL;
				}
			}
			N = _n; 
			M = _m; 
		} while (num_items_read != 2);
		NZ = N * M;
	}

	matrix_info->rows = N;
	matrix_info->cols = M;
	matrix_info->nnz = NZ;
	return matrix_info;
}

typedef void (*set_int_func) (void *array, size_t pos, int64_t value);

static void set_int32(void *array, size_t pos, int64_t value)
{
    int32_t *iarray = (int32_t *) array;
    iarray[pos] = (int32_t) value;
}

static void set_int64(void *array, size_t pos, int64_t value)
{
    int64_t *iarray = (int64_t*) array;
    iarray[pos] = value;
}

/*-----------------------------------------------------------------------------
 *  Read matrix
 *-----------------------------------------------------------------------------*/
static int uf_matrix_coord_internal(uf_field_t *field, void *nrows, void *ncols, void *  nnz, void **rowptr, void **colptr, void **values, uf_matrix_t *mat, set_int_func si, size_t intsize)
{
    mtx_info_t *mat_info =  NULL;
    uf_io_file_t *file = NULL;
    char line[LINE_LENGTH];
    char *cret = NULL;

    /* Open */ 
    file = uf_matrix_open(mat);
    if ( file == NULL) {
        fprintf(stderr, "Open matrix failed.\n");
        return -1;
    }

    /* Parse Header  */
    mat_info = uf_matrix_read_info(file);
    if (mat_info == NULL ){
        fprintf(stderr, "Failed to read MTX header.\n");
        uf_io_close(file);
        return -1;
    }

    /* Copy header  */
    *field = mat_info->data_type;
    si(nrows, 0, mat_info->rows);
    si(ncols, 0, mat_info->cols);
    si(nnz, 0, mat_info->nnz);

    // printf("field: %d, nrows: %d, ncols: %d, nnz: %d\n", *field, *nrows, *ncols, *nnz);
    // printf("Intsize: %d\n", intsize*8);
    if (!MTX_IS_COORD(mat_info)) {
        fprintf(stderr, "Only coordinate matrices are supported at the moment.\n");
        goto failed;
    }

    /* Read Sparse Data */ 
    void *irowptr, *icolptr;
    int64_t innz;
    if (MTX_IS_GENERAL(mat_info)) 
        innz = mat_info->nnz;
    else 
        innz = 2  * mat_info->nnz;

    irowptr =  malloc(intsize  * (innz));
    icolptr =  malloc(intsize  * (innz));
    if ( irowptr == NULL || icolptr == NULL ){
        free(irowptr); free(icolptr);
        fprintf(stderr, "Alloate Row/Column pointer failed.\n");
        goto failed;
    }
    innz=0;

    /* Real data */
    if (*field == UF_REAL ) {
        double * ivalues;
        int64_t _r, _c, i;
        int ret; 

        ivalues = (double *) malloc(sizeof(double) * ((size_t)mat_info->nnz) * 2); 
        if (ivalues == NULL) {
            free(irowptr); free(icolptr);
            fprintf(stderr, "Allocate values failed.\n");
            goto failed;
        }

        for (i=0; i < mat_info->nnz; i++){
            cret = uf_io_gets(line, LINE_LENGTH,file);
            ret = sscanf(line, "%"PRId64" %"PRId64" %lg\n", &_r, &_c, &(ivalues[innz]));
            if ( ret != 3 || !cret){
                fprintf(stderr, "cannot read value for real format in line %ld  - \"%s\"\n", (long) i,line);
                free(irowptr); free(icolptr); free(ivalues);
                goto failed;
            }
            si(irowptr, innz, _r-1);
            si(icolptr, innz, _c-1);

            if ( MTX_IS_SYMMETRIC(mat_info) && _r != _c ){
                // printf("Symmetric\n");
                innz++; 
                si(irowptr, innz, _c-1);
                si(icolptr, innz, _r-1);
                ivalues[innz]=ivalues[innz-1]; 				
            }
            if ( MTX_IS_SKEWSYMMETRIC(mat_info) && _r != _c ){
                innz++; 
                si(irowptr, innz, _c-1);
                si(icolptr, innz, _r-1);
                ivalues[innz]=-ivalues[innz-1]; 				
            }
            innz ++; 
        }
        si(nnz, 0, innz);
        // printf("innz: %ld\n", innz);
        irowptr = realloc(irowptr, intsize *((size_t) innz)); 
        icolptr = realloc(icolptr, intsize *((size_t) innz)); 
        ivalues = (double * ) realloc(ivalues, sizeof(double) *((size_t) innz)); 
        *rowptr = irowptr;
        *colptr = icolptr;
        *values = ivalues;
        /* Complex data  */
    } else if (*field == UF_COMPLEX ) {
        double complex * ivalues;
        double re, im;
        int64_t _r, _c, i;
        int ret; 

        ivalues = (double complex *) malloc(sizeof(double complex) * ((size_t)mat_info->nnz) * 2); 
        if (ivalues == NULL) {
            free(irowptr); free(icolptr);
            fprintf(stderr, "Allocate values failed.\n");
            goto failed;
        }

        for (i=0; i < mat_info->nnz; i++){
            cret = uf_io_gets(line, LINE_LENGTH,file);
            ret = sscanf(line, "%"PRId64" %"PRId64" %lg %lg\n", &_r, &_c, &re, &im);
            if ( ret != 4 || !cret){
                fprintf(stderr, "cannot read value for complex format in line %ld  - \"%s\"\n", (long) i,line);
                free(irowptr); free(icolptr); free(ivalues);
                goto failed;
            }

            si(irowptr, innz, _r-1);
            si(icolptr, innz, _c-1);
            ivalues[innz] = re + im * I;

            if ( MTX_IS_SYMMETRIC(mat_info) && _r != _c ){
                innz++; 
                si(irowptr, innz, _c-1);
                si(icolptr, innz, _r-1);
                ivalues[innz]=ivalues[innz-1]; 
            }
            if ( MTX_IS_HERMITIAN(mat_info) && _r != _c ){
                innz++;
                si(irowptr, innz, _c-1);
                si(icolptr, innz, _r-1);
                ivalues[innz]=conj(ivalues[innz-1]); 
            }

            if ( MTX_IS_SKEWSYMMETRIC(mat_info) && _r != _c ){
                innz++; 
                si(irowptr, innz, _c-1);
                si(icolptr, innz, _r-1);
                ivalues[innz]=-ivalues[innz-1]; 				
            }
            innz ++; 
        }
        si(nnz, 0,innz);
        irowptr = realloc(irowptr, intsize *((size_t) innz)); 
        icolptr = realloc(icolptr, intsize *((size_t) innz)); 
        ivalues = (double complex* ) realloc(ivalues, sizeof(double complex) *((size_t) innz)); 
        *rowptr = irowptr;
        *colptr = icolptr;
        *values = ivalues;
    }
    /* Pattern */
    else if (*field == UF_PATTERN ) {
        double * ivalues;
        int64_t _r, _c, i;
        int ret; 

        ivalues = (double *) malloc(sizeof(double) * ((size_t)mat_info->nnz) * 2); 
        if (ivalues == NULL) {
            free(irowptr); free(icolptr);
            fprintf(stderr, "Allocate values failed.\n");
            goto failed;
        }

        for (i=0; i < mat_info->nnz; i++){
            cret = uf_io_gets(line, LINE_LENGTH,file);
            ret = sscanf(line, "%"PRId64" %"PRId64"\n", &_r, &_c);
            if ( ret != 2 || !cret){
                fprintf(stderr, "cannot read value for pattern format in line %ld  - \"%s\"\n", (long) i,line);
                free(irowptr); free(icolptr); free(ivalues);
                goto failed;
            }

            si(irowptr, innz, _r-1);
            si(icolptr, innz, _c-1);
            ivalues[innz] = 1.0;

            if ( (MTX_IS_SYMMETRIC(mat_info) || MTX_IS_SKEWSYMMETRIC(mat_info) || MTX_IS_HERMITIAN(mat_info)) && ( _r != _c ))
            {
                innz++;
                si(irowptr, innz, _c-1);
                si(icolptr, innz, _r-1);
                ivalues[innz]=ivalues[innz-1]; 
            }
            innz ++; 
        }
        si(nnz, 0, innz);
        irowptr =  realloc(irowptr, intsize *((size_t) innz)); 
        icolptr =  realloc(icolptr, intsize *((size_t) innz)); 
        ivalues = (double *) realloc(ivalues, sizeof(double ) *((size_t) innz)); 
        *rowptr = irowptr;
        *colptr = icolptr;
        *values = ivalues;
    }
    /* Unknown */
    else {
        fprintf(stderr, "Unknown data type.\n");
        free(icolptr); free(irowptr);
        goto failed;
    }
    free(mat_info);
    uf_io_close(file);
    return 0; 
failed:
    free(mat_info);
    uf_io_close(file);
    return -1;

}

/*-----------------------------------------------------------------------------
 *  Int32 and Int64 interface to the coordinate read function. 
 *-----------------------------------------------------------------------------*/
int uf_matrix_coord_int32(uf_field_t *field, int32_t *nrows, int32_t *ncols, int32_t *nnz, int32_t **rowptr, int32_t **colptr, void **values, uf_matrix_t *mat)
{
    return uf_matrix_coord_internal(field, (void *) nrows, (void *)ncols, (void *) nnz, 
               (void **) rowptr, (void **) colptr, (void **) values, mat, set_int32, sizeof(int32_t));
}

int uf_matrix_coord_int64(uf_field_t *field, int64_t *nrows, int64_t *ncols, int64_t *nnz, int64_t **rowptr, int64_t **colptr, void **values, uf_matrix_t *mat)
{
    return uf_matrix_coord_internal(field, (void *) nrows, (void *)ncols, (void *) nnz, 
               (void **) rowptr, (void **) colptr, (void **) values, mat, set_int64, sizeof(int64_t));

}


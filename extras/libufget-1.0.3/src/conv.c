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
#include <assert.h>
#include <complex.h>

#include <curl/curl.h>
#include <sqlite3.h>
#include <matio.h>
#include <archive.h>
#include <archive_entry.h>

#include "libufget.h"
#include "uf_internal.h"

static int __uf_matrix_conv_memsave = 0; 
int uf_matrix_conv_memsave(int save) 
{
    int old = __uf_matrix_conv_memsave;
    __uf_matrix_conv_memsave = save; 
    return old;
}

// #define LIBC_SORT
typedef void (*set_int_func) (void *array, size_t pos, int64_t value);
typedef int64_t (*get_int_func) (void *array, size_t pos); 

static void set_int32(void *array, size_t pos, int64_t value)
{
    int32_t *iarray = (int32_t *) array;
    iarray[pos] = (int32_t) value;
}

static int64_t get_int32(void *array, size_t pos) {
    int32_t *iarray = (int32_t *) array;
    int64_t ret = iarray[pos]; 
    return ret; 
}

static void set_int64(void *array, size_t pos, int64_t value)
{
    int64_t *iarray = (int64_t*) array;
    iarray[pos] = value;
}

static int64_t get_int64(void *array, size_t pos) {
    int64_t *iarray = (int64_t *) array;
    return iarray[pos]; 
}

#ifdef LIBC_SORT
#warning COORD -> CSR conversion does not work inplace 
typedef struct _mat_entry_t {
    int64_t row; 
    int64_t col; 
    double re; 
    double im; 
} mat_entry_t;

static int comp_csr(const void *_a, const void *_b)
{
    mat_entry_t* a = (mat_entry_t *) _a; 
    mat_entry_t* b = (mat_entry_t *) _b; 

    if ( a->row < b -> row ) {
        return -1; 
    } 
    if ( a->row > b->row ) {
        return 1; 
    }
    if ( a-> row == b->row ) {
        if ( a->col < b->col) {
            return -1; 
        } else if ( a->col > b->col ) {
            return 1; 
        }
        return 0; 
    }
    return 0; 
}
#else 
static int comp_csr2(int64_t arow, int64_t acol, int64_t brow, int64_t bcol)
{
    // printf("COMP: %ld %ld <-> %ld %ld\n", arow, acol, brow, bcol);
    if ( arow < brow ) {
        return -1; 
    } 
    if ( arow > brow ) {
        return 1; 
    }
    if ( arow == brow ) {
        if ( acol < bcol) {
            return -1; 
        } else if ( acol > bcol ) {
            return 1; 
        }
        return 0; 
    }
    return 0; 
}


#define EXCH(i,j) r = rowptr[i]; \
        c = colptr[i]; \
        v = values[i]; \
        rowptr[i] = rowptr[j]; \
        colptr[i] = colptr[j]; \
        values[i] = values[j]; \
        rowptr[j] = r; \
        colptr[j] = c; \
        values[j] = v; 

#define COMPEXCH(i,j) if ( comp_csr2(rowptr[i], colptr[i], rowptr[j],colptr[j]) == -1){\
        r = rowptr[i]; \
        c = colptr[i]; \
        v = values[i]; \
        rowptr[i] = rowptr[j]; \
        colptr[i] = colptr[j]; \
        values[i] = values[j]; \
        rowptr[j] = r; \
        colptr[j] = c; \
        values[j] = v; \
}

typedef struct {
    size_t len; 
    void *rptr; 
    void *cptr; 
    void *vptr; 
} stack_node; 

#define SORT_THRESH 32 
#define STACK_SIZE  64
#define PUSH(n, rowptr, colptr, values) ((void) ((top->len  = (n)), (top->rptr = (rowptr)), (top->cptr = (colptr)), (top->vptr = (values)), ++top))
#define POP(n, rowptr, colptr, values)  ((void) (--top, (n = top->len), (rowptr = top->rptr)), (colptr = top->cptr), (values = top->vptr))
#define STACK_NOT_EMPTY (stack < top)


/*-----------------------------------------------------------------------------
 *  Quick Sort 
 *-----------------------------------------------------------------------------*/
#define QUICKSORT(NAME,INTTYPE,FLOAT) static void qs ## NAME(size_t n, INTTYPE *rowptr, INTTYPE  *colptr, FLOAT *values) \
{\
    stack_node stack[STACK_SIZE];\
    stack_node *top = stack; \
    ssize_t i, j, p;\
    INTTYPE r,c; \
    FLOAT v;\
    if (n <= SORT_THRESH) return;\
    PUSH(0, NULL, NULL, NULL); \
    while ( STACK_NOT_EMPTY ) \
    {\
        if ( n < SORT_THRESH ){\
            POP(n, rowptr, colptr, values); \
            continue; \
        } \
        p = n-1;\
        if (n > 3) { \
            EXCH(n/2,p-1); \
            COMPEXCH(0,p-1);\
            COMPEXCH(0,p); \
            COMPEXCH(p-1,p);\
        } else {\
            EXCH(n/2,p);\
        }\
        for (i = 0, j = n - 1; ; i++, j--) { \
            while (comp_csr2(rowptr[i], colptr[i], rowptr[p],colptr[p]) == -1)\
                i++;\
            while (comp_csr2(rowptr[p], colptr[p], rowptr[j],colptr[j]) == -1){\
                j--;\
                if(j == i) break;\
            }\
            if (i >= j) break;\
            EXCH(i,j);\
        } \
        EXCH(i,p); \
        if ( i > n-i-1) { \
            PUSH(i, rowptr, colptr, values);\
            n = n-i-1; \
            rowptr += i+1; \
            colptr += i+1; \
            values += i+1; \
        } else {\
            PUSH(n-i-1, rowptr + i +1, colptr+i+1, values+i+1);\
            n = i; \
        }\
    }\
}

/*-----------------------------------------------------------------------------
 *  Insertion Sort 
 *-----------------------------------------------------------------------------*/
#define INSERTIONSORT(NAME,INTTYPE,FLOAT) static void is ## NAME(size_t n, INTTYPE *rowptr, INTTYPE  *colptr, FLOAT *values) \
{\
    ssize_t i, j; \
    INTTYPE r,c; \
    FLOAT v; \
    for(i = 1; i < n; i++) \
    {\
        r = rowptr[i]; \
        c = colptr[i]; \
        v = values[i]; \
        for (j = i - 1; j >= 0 && comp_csr2(r,c, rowptr[j], colptr[j]) == -1 ; j--)\
        {\
            rowptr[j+1] = rowptr[j]; \
            colptr[j+1] = colptr[j]; \
            values[j+1] = values[j]; \
        }\
        rowptr[j+1] = r; \
        colptr[j+1] = c; \
        values[j+1] = v; \
    }\
}


/*-----------------------------------------------------------------------------
 *  Hybrid Sort 
 *-----------------------------------------------------------------------------*/
#define HYBRIDSORT(NAME, INTTYPE, FLOAT) static void NAME(size_t n, INTTYPE *rowptr, INTTYPE  *colptr, FLOAT *values) \
{\
    qs ## NAME ( n, rowptr, colptr, values); \
    is ## NAME ( n, rowptr, colptr, values); \
}

QUICKSORT(sort_real_32, int32_t, double)
QUICKSORT(sort_real_64, int64_t, double)
QUICKSORT(sort_cpx_32, int32_t, double complex)
QUICKSORT(sort_cpx_64, int64_t, double complex)
/* Insert Sort   */
INSERTIONSORT(sort_real_32, int32_t, double)
INSERTIONSORT(sort_real_64, int64_t, double)
INSERTIONSORT(sort_cpx_32, int32_t, double complex)
INSERTIONSORT(sort_cpx_64, int64_t, double complex)
/* The sort function   */
HYBRIDSORT(sort_real_32, int32_t, double)
HYBRIDSORT(sort_real_64, int64_t, double)
HYBRIDSORT(sort_cpx_32, int32_t, double complex)
HYBRIDSORT(sort_cpx_64, int64_t, double complex)

#endif 

typedef struct _vmat_entry_t {
    int64_t col; 
    double re; 
    double im; 
} vmat_entry_t;

static int comp_row(const void *_a, const void *_b)
{
    vmat_entry_t* a = (vmat_entry_t *) _a; 
    vmat_entry_t* b = (vmat_entry_t *) _b; 

    if ( a->col < b -> col ) {
        return -1; 
    } 
    if ( a->col > b->col ) {
        return 1; 
    }
    return 0; 
}


/* #define PRINT(NAME, INTTYPE) void NAME ( size_t N, INTTYPE *rowptr, INTTYPE *colptr) {\
    size_t i; \
    for (i = 0; i < N; i++) {\
    printf("%6" PRId64 "\t %6"PRId64"\n", (int64_t) rowptr[i], (int64_t) colptr[i]);\
    }\
    }

    PRINT(pr32, int32_t)
    PRINT(pr64, int64_t) */

/*-----------------------------------------------------------------------------
 *  generic convert to csr function 
 *-----------------------------------------------------------------------------*/
static int uf_matrix_coord_to_csr_intX(uf_field_t field, int64_t nrows, int64_t ncols, int64_t nnz, void **rowptr, void **colptr, void **values, 
        set_int_func si, get_int_func gi,  size_t intsize )
{
    void *irowptr = *rowptr; 
    void *icolptr = *colptr; 

    if ( __uf_matrix_conv_memsave == 0 ) {
        int64_t *rowcount; 
        int64_t *orowptr, *ocolptr; 
        int64_t pos, i, j; 
        double *val, *oval; 
        vmat_entry_t *row; 


        rowcount = (int64_t *) malloc(sizeof(int64_t) * (nrows)); 
        orowptr  = (int64_t *) malloc(sizeof(int64_t) * (nrows+1)); 
        ocolptr  = (int64_t *) malloc(sizeof(int64_t) * (nnz)); 
        row      = (vmat_entry_t *) malloc(sizeof(vmat_entry_t) * (ncols)); 

        if ( field == UF_REAL || field == UF_PATTERN ) {
            oval = (double *) malloc(sizeof(double) * (nnz)); 
        } else {
            oval = (double *) malloc(sizeof(double) * (nnz*2)); 
        }
        val = (double *) *values; 

        // count entries in rows
        for (i = 0; i < nrows; i++)  rowcount[i] = 0;    
        for ( i = 0 ; i < nnz; i++) rowcount[gi(irowptr,i)]++;

        // fill rowptr
        orowptr[0] = 0;
        for ( i = 0; i < nrows; i++){
            orowptr[i+1] = orowptr[i] + rowcount[i];
            rowcount[i] = orowptr[i];
        }
        // fill values
        if (field == UF_REAL || field == UF_PATTERN) {
            for ( i = 0; i < nnz; i++){
                pos =  rowcount[gi(irowptr,i)];
                ocolptr[pos] = gi(icolptr,i);
                oval[pos] = val[i];
                rowcount[gi(irowptr,i)]++;
            }

            for (i = 0; i < nrows; i++) {
                pos = 0; 
                for (j = orowptr[i]; j < orowptr[i+1]; j++) {
                    row[pos].col = ocolptr[j]; 
                    row[pos].re  = oval[j]; 
                    pos++; 
                }       
                qsort(row, pos, sizeof(vmat_entry_t), comp_row); 
                pos = 0 ; 
                for (j = orowptr[i]; j < orowptr[i+1]; j++) {
                    si(icolptr,j, row[pos].col);  
                    val[j] = row[pos].re;  
                    pos++; 
                }       
            }
        } else {
            for ( i = 0; i < nnz; i++){
                pos =  rowcount[gi(irowptr,i)];
                ocolptr[pos] = gi(icolptr,i);
                oval[2*pos] = val[2*i];
                oval[2*pos+1] = val[2*i+1];
                rowcount[gi(irowptr,i)]++;
            }
            for (i = 0; i < nrows; i++) {
                pos = 0; 
                for (j = orowptr[i]; j < orowptr[i+1]; j++) {
                    row[pos].col = ocolptr[j]; 
                    row[pos].re  = oval[2*j]; 
                    row[pos].im  = oval[2*j+1]; 
                    pos++; 
                }       
                qsort(row, pos, sizeof(vmat_entry_t), comp_row); 
                pos = 0 ; 
                for (j = orowptr[i]; j < orowptr[i+1]; j++) {
                    si(icolptr,j, row[pos].col);  
                    val[2*j] = row[pos].re;  
                    val[2*j+1] = row[pos].re;  
                    pos++; 
                }       
            }
        }
        irowptr = realloc(irowptr, intsize * (nrows+1)); 
        for (i = 0; i < nrows+1; i++) {
            si(irowptr, i, orowptr[i]) ;         
        }

        *rowptr = irowptr; 
        free(rowcount); 
        free(orowptr); 
        free(ocolptr); 
        free(oval); 
        free(row); 
        return 0; 

    } else { /* Save memory but slower   */ 
        double *val; 
        double complex *valc; 
        int64_t i = 0, irow; 

#if defined(LIBC_SORT)
        mat_entry_t *mat = NULL; 
        mat = (mat_entry_t *) malloc(sizeof(mat_entry_t) * (nnz)); 
        if ( field == UF_REAL || field == UF_PATTERN ) {
            val = (double *) *values; 
            for (i = 0; i < nnz; i++) {
                mat[i].row = gi(irowptr, i); 
                mat[i].col = gi(icolptr, i); 
                mat[i].re  = val[i]; 
            }
            qsort(mat, nnz, sizeof(mat_entry_t), comp_csr); 
            for (i = 0; i < nnz; i++) {
                si(irowptr, i, mat[i].row);  
                si(icolptr, i, mat[i].col);  
                val[i] = mat[i].re; 
            }
            /* if (intsize == 4 ) {
               pr32(nnz, irowptr, icolptr); 
               } else {
               pr64(nnz, irowptr, icolptr); 
               } */
        } 
        else {
            valc= (double complex *) *values; 
            for (i = 0; i < nnz; i++) {
                mat[i].row = gi(irowptr, i); 
                mat[i].col = gi(icolptr, i); 
                mat[i].re  = creal(valc[i]); 
                mat[i].im  = cimag(valc[i]); 
            }
            qsort(mat, nnz, sizeof(mat_entry_t), comp_csr); 
            for (i = 0; i < nnz; i++) {
                si(irowptr, i, mat[i].row);  
                si(icolptr, i, mat[i].col);  
                valc[i] = mat[i].re + mat[i].im*I; 
            }

        }
        free(mat); 
#else
        if ( field == UF_REAL || field == UF_PATTERN ) {
            val = (double *) *values; 

            if (intsize == 4 ) {
                sort_real_32(nnz, irowptr, icolptr, val); 
            } else {
                sort_real_64(nnz, irowptr, icolptr, val); 
            }
        } else {
            valc= (double complex *) *values; 
            if (intsize == 4 ) {
                sort_cpx_32(nnz, irowptr, icolptr, valc); 
            } else {
                sort_cpx_64(nnz, irowptr, icolptr, valc); 
            }
        }
#endif

        /* Sum up  */
        void * irowptr2 = malloc(intsize * (nrows+1)); 
        if (nrows+1 >= nnz) {
            irowptr = realloc(irowptr, intsize * (nrows+1)); 
        }
        irow = 0; 
        i = 0;
        for (irow = 0; irow < nrows; irow++) {
            while ( i < nnz && gi(irowptr, i) == irow ) { 
                i++; 
            }
            si(irowptr2, irow+1, i); 
        }
        if ( i!= nnz ) {
            fprintf(stderr, "i!=nnz: %ld - %ld\n", i, nnz);
            free(irowptr2); 
            return -1; 
        }
        // assert(i==nnz); 
        si(irowptr2,0,0); 
        si(irowptr2,nrows, nnz);
        irowptr = realloc(irowptr, intsize * (nrows+1)); 
        memcpy(irowptr, irowptr2, intsize*(nrows+1)); 
        free(irowptr2); 
        *rowptr = irowptr; 

        return 0; 
    }
    return 0; 
}

/*-----------------------------------------------------------------------------
 *  Interface coord -> csr  with 4 byte integers
 *-----------------------------------------------------------------------------*/
int uf_matrix_coord_to_csr_int32(uf_field_t field, int32_t nrows, int32_t ncols, int32_t nnz, int32_t **rowptr, int32_t **colptr, void **values)
{
    return uf_matrix_coord_to_csr_intX(field, nrows, ncols, nnz, (void **) rowptr, (void **) colptr, (void **) values,  set_int32, get_int32, sizeof(int32_t)); 
}

/*-----------------------------------------------------------------------------
 *  Interface coord -> csr  with 8 byte integers
 *-----------------------------------------------------------------------------*/
int uf_matrix_coord_to_csr_int64(uf_field_t field, int64_t nrows, int64_t ncols, int64_t nnz, int64_t **rowptr, int64_t **colptr, void **values)
{
    return  uf_matrix_coord_to_csr_intX(field, nrows, ncols, nnz, (void **) rowptr, (void **) colptr, (void **) values,  set_int64, get_int64, sizeof(int64_t)); 

}

/*-----------------------------------------------------------------------------
 *  Interface coord -> csc  with 4 byte integers
 *-----------------------------------------------------------------------------*/
int uf_matrix_coord_to_csc_int32(uf_field_t field, int32_t nrows, int32_t ncols, int32_t nnz, int32_t **rowptr, int32_t **colptr, void **values)
{
    return uf_matrix_coord_to_csr_intX(field, ncols, nrows, nnz, (void **) colptr, (void **) rowptr, (void **) values,  set_int32, get_int32, sizeof(int32_t)); 
}

/*-----------------------------------------------------------------------------
 *  Interface coord -> csc  with 8 byte integers
 *-----------------------------------------------------------------------------*/
int uf_matrix_coord_to_csc_int64(uf_field_t field, int64_t nrows, int64_t ncols, int64_t nnz, int64_t **rowptr, int64_t **colptr, void **values)
{
    return  uf_matrix_coord_to_csr_intX(field, ncols, nrows, nnz, (void **) colptr, (void **) rowptr, (void **) values,  set_int64, get_int64, sizeof(int64_t)); 
}


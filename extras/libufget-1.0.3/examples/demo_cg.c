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

/* BLAS  */
#define BlasInt int
void dcopy_(BlasInt * N, double *b, BlasInt *INCX, double *p, BlasInt *INCY);
void daxpy_(BlasInt * N, double *da, double *x, BlasInt *INCX, double *y, BlasInt *INCY);
double dnrm2_(BlasInt *N, double *X, BlasInt *INCX);
double ddot_(BlasInt *N, double *X, BlasInt *incx, double *y, BlasInt *incy);
void dscal_(BlasInt *n , double *alpha, double *v, BlasInt * INCX);

/* Sparse MVP  */
void sparse_mvp(int64_t rows, int64_t *rowptr, int64_t *colptr, double *values, double * x, double *y ) {
    int64_t i , j;
    for ( i = 0 ; i < rows; i++) {
	    y[i]=0;
        for ( j = rowptr[i]; j < rowptr[i+1]; j++) {
            y[i] += values[j] * x[colptr[j]];
        }
    }
}
/* CG   */
int cg(int64_t rows, int64_t*rowptr, int64_t*colptr, double *values, double *x, double *b, int maxit, double tol)
{
    int n = rows, INCX=1, INCY=1,m;
    double *p, *r, *v, alpha, alphaOld, lambda,lambdamin, k, da=-1.0, one=1;
    double nrmB;
    double r0;
    p = (double *)malloc((n)* sizeof(double));
    r = (double *)malloc((n)* sizeof(double));
    v = (double *)malloc((n)* sizeof(double));


    nrmB = dnrm2_(&n, b, &INCX);
    sparse_mvp(rows, rowptr, colptr, values, x, v);
    dcopy_(&n, b, &INCX, p, &INCY);
    daxpy_(&n, &da, v, &INCX, p, &INCY);
    dcopy_(&n, p, &INCX, r, &INCY);
    r0 = dnrm2_(&n,r,&INCX);
    alpha = r0*r0;
    for(m = 0; m < maxit; m++)
    {
        if ( alpha == 0.0 ) break;

        sparse_mvp(rows, rowptr, colptr, values, p, v);
        lambda = alpha / ddot_(&n, v,&INCX, p, &INCY);
        lambdamin=-lambda;
        daxpy_(&n, &lambda, p, &INCX, x, &INCY);
        daxpy_(&n, &lambdamin, v, &INCX, r, &INCY);
        alphaOld = alpha;
        r0 = dnrm2_(&n,r,&INCX);
        alpha = r0 * r0;
        k=alpha / alphaOld;
        dscal_(&n, &k, p, &INCX);
        daxpy_(&n, &one, r, &INCX, p, &INCY);
        // printf("it = %5d \t res = %lg\n",m, r0/nrmB);
        if ( r0 < nrmB *tol ) break;
    }
    free(p);
    free(r);
    free(v);
    return m;
}


#define MAX(a,b) ((a)<(b)?(b):(a))
#define MIN(a,b) ((a)<(b)?(a):(b))
int main(int argc, char **argv)
{
    uf_collection_t *col;
    uf_matrix_t *mat = NULL;
    uf_query_t *query;
    int i,j, it;

    printf("libufget - CG demo\n");
    printf("==================\n");



    /* Init Collection  */

    /* Default Init  */
    col = uf_collection_init();

    /* Init with flags */
    // col = uf_collection_init1(UF_COLLECTION_VERBOSE|UF_COLLECTION_CLEANUP);

    printf("Database contains %d matrices.\n", uf_collection_num_matrices(col));

    /*  Pre download all matrices for testing */
    for (i = 0; i < 30; i++) {
        mat = uf_collection_get_by_id(col, i+1);
        if (mat) uf_matrix_get(mat);
        uf_matrix_free(mat);
    }

    /* SQL Query */
    query = uf_query_sql(col, "posdef = 1.0 AND numerical_symmetry = 1.0 AND isReal=1 AND nnz < 20000000 AND id < 30 ");
    /*   */

    while ( (mat = uf_query_next(query)) != NULL ) {

        uf_field_t field;
        int64_t rows, cols, nnz;
        int64_t *rowptr, *colptr;
        double *values;
        double * x , *y;
        double err, nrm;

        /* Load the matrix  */
        printf("Perform CG on...\n");
        uf_matrix_print(mat, 1);

        /* Download marix and skip if not possible.  */
        if ( uf_matrix_get(mat) != 0 ){
            printf("Cannot download Matrix (%s/%s) skipping.\n", mat->group_name, mat->name);
            uf_matrix_free(mat);
            continue;
        }


        if ( uf_matrix_coord_int64(&field, &rows, &cols, &nnz, &rowptr, &colptr, (void **) &values, mat) != 0) {
            printf("Cannot open Matrix (%s/%s) skipping.\n", mat->group_name, mat->name);
            uf_matrix_free(mat);
            continue;

        }
        if (field != UF_REAL) {
            printf("--> matrix is not real. next.\n");
            uf_matrix_free(mat);
            continue;
        }

        uf_matrix_coord_to_csr_int64(field, rows, cols, nnz, &rowptr, &colptr, (void **) &values);

        /* Setup RHS ...  */
        x = malloc(sizeof(double)*rows);
        y = malloc(sizeof(double)*rows);
        for (j = 0 ; j < rows ; j++) {
            x[j]=j+1; y[j]=0;
        }
        sparse_mvp(rows, rowptr, colptr, values, x, y);

        for ( i = 0 ; i < rows ; i ++) {
            x[i] = 1;
        }


        /* Solve   */
        it = cg (rows, rowptr, colptr, values, x, y, MIN(1000,rows*10),1e-7);


        /* Check */
        err = 0.0;
        nrm = 0.0;
        for (j = 0; j < rows; j++) {
            err += pow(x[j]-(j+1),2);
            nrm += pow(j+1,2);
        }
        nrm = sqrt(nrm);
        err = sqrt(err);
        printf("IT = %6d \t\tForward Error: %lg\n\n", it, err/nrm);


        free(x);
        free(y);
        free(rowptr);
        free(colptr);
        free(values);
        uf_matrix_free(mat);
        uf_collection_cleanup(col);

    }
    // uf_collection_cleanup(col);
    /* Clean up  */
    uf_query_free(query);
    uf_collection_finalize(col);

    return 0;
}

#include <stdio.h>
#include <stdlib.h>

#include "libufget.h"

#include <unistd.h>
#include <sys/time.h>

#ifndef __USE_POSIX199309
#define __USE_POSIX199309
#include <time.h>
#undef __USE_POSIX199309
#else
#include <time.h>
#endif


double wtime ()
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  return tv.tv_sec + tv.tv_usec / 1e6;
}



int main(int argc, char **argv)
{
    uf_collection_t *col; 
    uf_matrix_t *mat = NULL; 
    int i; 
    double ts, te, te_read; 

    col = uf_collection_init1(UF_COLLECTION_VERBOSE);

    printf("Database contains %d matrices.\n", uf_collection_num_matrices(col));

    for (i = 0; i < uf_collection_num_matrices(col) ; i++) {
        mat = uf_collection_get_by_id(col, i+1); 

        /*  Load as Int32  */
        { 
            uf_field_t field;
            int32_t rows, cols, nnz;
            int32_t *rowptr, *colptr;
            double *values;

            printf("Matrix[%4d] %s/%s\n", mat->id, mat->group_name, mat->name);
            ts = wtime(); 
            uf_matrix_coord_int32(&field, &rows, &cols, &nnz, &rowptr, &colptr, (void **) &values, mat);
            te_read = wtime() - ts; 
            ts = wtime();         
            uf_matrix_coord_to_csr_int32(field, rows, cols, nnz, &rowptr, &colptr, (void **) &values); 
            te = wtime()-ts; 

            printf(" - field: %2d \t rows: %10d\t cols: %10d\t nnz: %12d  \t time_read: %10lgs \ttime_conv: %10lg rates( %8lg / %8lg )\n", field, (int) rows, (int) cols, (int)nnz, te_read, te, nnz/te_read, nnz/te);
            free(rowptr);
            free(colptr);
            free(values);
        }
        /* Load as Int64 */
        { 
            uf_field_t field;
            int64_t rows, cols, nnz;
            int64_t *rowptr, *colptr;
            double *values;

            ts = wtime(); 
            uf_matrix_coord_int64(&field, &rows, &cols, &nnz, &rowptr, &colptr, (void **) &values, mat);
            te_read = wtime() - ts; 
            ts = wtime(); 
            uf_matrix_coord_to_csr_int64(field, rows, cols, nnz, &rowptr, &colptr, (void **) &values); 
            te = wtime()-ts; 
            printf(" - field: %2d \t rows: %10d\t cols: %10d\t nnz: %12d  \t time_read: %10lgs \ttime_conv: %10lg rates( %8lg / %8lg )\n", field, (int) rows, (int) cols, (int)nnz, te_read, te, nnz/te_read, nnz/te);

            free(rowptr);
            free(colptr);
            free(values);

        }


        if (mat) uf_matrix_get(mat); 
        uf_matrix_free(mat);        
    }

    uf_collection_finalize(col); 

    return 0;
}

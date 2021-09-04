// dummyFactorize.h
// Dummy ILU0 Factorization
// Created: 05-12-2019
// Author: Najeeb

#ifndef __HEADER_DUMMYFACTORIZE__
#define __HEADER_DUMMYFACTORIZE__

#include "matrix.h"
#include "common.h"

int dummyFactorize(csr_matrix *A, csr_matrix *L, csr_matrix *U, int substitution)
{
	ind_type *csrRowPtr_tmp_L = (ind_type *)malloc((A->m+1) * sizeof(ind_type));
    ind_type *csrColIdx_tmp_L = (ind_type *)malloc((A->m+A->nnz) * sizeof(ind_type));
    val_type *csrVal_tmp_L    = (val_type *)malloc((A->m+A->nnz) * sizeof(val_type));

    ind_type *csrRowPtr_tmp_U = (ind_type *)malloc((A->m+1) * sizeof(ind_type));
    ind_type *csrColIdx_tmp_U = (ind_type *)malloc((A->m+A->nnz) * sizeof(ind_type));
    val_type *csrVal_tmp_U    = (val_type *)malloc((A->m+A->nnz) * sizeof(val_type));

    int nnz_pointer_L = 0;
    int nnz_pointer_U = 0;
    csrRowPtr_tmp_L[0] = 0;
    csrRowPtr_tmp_U[0] = 0;
    val_type valTemp = 0.0;
	printf("A->nnz: %d\n",A->nnz);
    
    for (int i = 0; i < A->m; i++)
    {   
        for (int j = A->csrRowPtr[i]; j < A->csrRowPtr[i+1]; j++)
        {   
            valTemp = rand() % 10 + 1;
            A->csrVal[j] = valTemp;

            if (substitution == SUBSTITUTION_FORWARD)
            {
                if (A->csrColIdx[j] < i)
                {
                    csrColIdx_tmp_L[nnz_pointer_L] = A->csrColIdx[j];                 
                    //csrColIdx_tmp_L[nnz_pointer_L];                    
                    //csrVal_tmp_L[nnz_pointer_L] = A->csrVal[j]; //csrValA[j]; 
                    csrVal_tmp_L[nnz_pointer_L] = valTemp;
                    nnz_pointer_L++;                    
                }
                
            }            
            else if (substitution == SUBSTITUTION_BACKWARD)
            {
                if (A->csrColIdx[j] > i)
                {
                    csrColIdx_tmp_U[nnz_pointer_U] = A->csrColIdx[j];
                    //csrVal_tmp_U[nnz_pointer_U] = A->csrVal[j]; //csrValA[j]; 
                    csrVal_tmp_U[nnz_pointer_U] = valTemp;
                    nnz_pointer_U++;
                }                
            }
            else if (substitution == SUBSTITUTION_BOTH)
            {
                if (A->csrColIdx[j] < i)
                {
                    csrColIdx_tmp_L[nnz_pointer_L] = A->csrColIdx[j];
                    //csrVal_tmp_L[nnz_pointer_L] = A->csrVal[j]; //csrValA[j]; 
                    csrVal_tmp_L[nnz_pointer_L] = valTemp;
                    nnz_pointer_L++;
                }	
                if (A->csrColIdx[j] > i)
                {
                    csrColIdx_tmp_U[nnz_pointer_U] = A->csrColIdx[j];
                    //csrVal_tmp_U[nnz_pointer_U] = A->csrVal[j]; //csrValA[j]; 
                    csrVal_tmp_U[nnz_pointer_U] = valTemp;
                    nnz_pointer_U++;
                }                
            }

        }
        // add dia nonzero
        if(substitution == SUBSTITUTION_FORWARD || substitution == SUBSTITUTION_BOTH)
		{
			csrColIdx_tmp_L[nnz_pointer_L] = i;
        	csrVal_tmp_L[nnz_pointer_L] = 1.0;            
        	nnz_pointer_L++;
        	csrRowPtr_tmp_L[i+1] = nnz_pointer_L;	
		}
        if(substitution == SUBSTITUTION_BACKWARD || substitution == SUBSTITUTION_BOTH)
		{
			csrColIdx_tmp_U[nnz_pointer_U] = i;
        	csrVal_tmp_U[nnz_pointer_U] = 1.0;
        	nnz_pointer_U++;

        	csrRowPtr_tmp_U[i+1] = nnz_pointer_U;	
		}        
    }
    //printf("A->m:%d\n", nnz_pointer_L);
    if(substitution == SUBSTITUTION_FORWARD || substitution == SUBSTITUTION_BOTH)
    {
    	int nnz_tmp_L = csrRowPtr_tmp_L[A->m];
    	int nnzTR = nnz_tmp_L;
    	csrColIdx_tmp_L = (ind_type *)realloc(csrColIdx_tmp_L, sizeof(val_type) * nnzTR);
        csrVal_tmp_L = (val_type *)realloc(csrVal_tmp_L, sizeof(val_type) * nnzTR);

        L->csrRowPtr = csrRowPtr_tmp_L;
        L->csrColIdx = csrColIdx_tmp_L;
        L->csrVal = csrVal_tmp_L;
        L->m = A->m;
        L->n = A->n;
        L->nnz = nnzTR;	
    }
    
    if(substitution == SUBSTITUTION_BACKWARD || substitution == SUBSTITUTION_BOTH)
    {
    	int nnz_tmp_U = csrRowPtr_tmp_U[A->m];
    	int nnzTR = nnz_tmp_U;	
    	csrColIdx_tmp_U = (ind_type *)realloc(csrColIdx_tmp_U, sizeof(val_type) * nnzTR);
        csrVal_tmp_U = (val_type *)realloc(csrVal_tmp_U, sizeof(val_type) * nnzTR);	
        
        U->csrRowPtr = csrRowPtr_tmp_U;
        U->csrColIdx = csrColIdx_tmp_U;
        U->csrVal = csrVal_tmp_U;
        U->m = A->m;
        U->n = A->n;
        U->nnz = nnzTR;	
    }
    
    return EXIT_SUCCESS;
}

#endif
// SILU_GPU.cpp
// Contains SILU GPU CUDA kernels and calls
// Author: Najeeb Ahmad
// Created: 24-12-2019

#include "SILU_GPU.h"

int SILU::SyncFreeAnalyzer(device_data_csc *device_csc_matrix, int *inDegree)
{
	int num_threads = 128;
    int num_blocks = ceil ((double)device_csc_matrix->nnz / (double)num_threads);
    cudaMemset(inDegree, 0, device_csc_matrix->m * sizeof(int));
       
    if(num_blocks >= 1)
    sptrsv_syncfree_cuda_analyser<<<num_blocks, num_threads>>>
                                (device_csc_matrix->d_cscRowIdx, 
                                device_csc_matrix->m, 
                                device_csc_matrix->nnz, 
                                inDegree);
    cudaError_t err = cudaGetLastError();
    if(err != cudaSuccess)
    {
        printf("CUDA SyncFree/SLFC Analyzer Error: %s\n", cudaGetErrorString(err));
        return EXIT_FAILURE;
    }
    else
        return EXIT_SUCCESS;
}

int SILU::SyncFreeExecutor(device_data_csc *device_csc_matrix, val_type *d_x, val_type *d_b, int *inDegree, int direction)
{
	int num_threads = WARP_PER_BLOCK * WARP_SIZE;
    int num_blocks = ceil ((double)device_csc_matrix->m / (double)(num_threads/WARP_SIZE));
    if(num_blocks >= 1)
    sptrsv_syncfree_cuda_executor<<<num_blocks, num_threads >>>
                                (device_csc_matrix->d_cscColPtr, device_csc_matrix->d_cscRowIdx,
                                 device_csc_matrix->d_cscVal, inDegree, 
                                 device_csc_matrix->d_left_sum,
                                 device_csc_matrix->m, direction, d_b, d_x);
    cudaError_t err = cudaGetLastError();
    if(err != cudaSuccess)
    {
        printf("CUDA SyncFree Executor Error: %s\n", cudaGetErrorString(err));
        return EXIT_FAILURE;
    }
    else
        return EXIT_SUCCESS;
}

int SILU::SLFCExecutor(device_data_csc *device_csc_matrix, val_type *d_x,  
                       val_type *diag, ind_type *jlev, int *inDegree)
{
    int num_threads = WARP_PER_BLOCK * WARP_SIZE;
    int num_blocks = ceil ((double)device_csc_matrix->m / (double)(num_threads/WARP_SIZE));

    if(num_blocks >= 1)
    {
        SLFCKernel<<<num_blocks, num_threads>>>
                                (device_csc_matrix->m, d_x, device_csc_matrix->d_cscVal, 
                                 device_csc_matrix->d_cscColPtr, device_csc_matrix->d_cscRowIdx,
                                 diag, inDegree, device_csc_matrix->jlev);
        cudaDeviceSynchronize();
    }
    
    cudaError_t err = cudaGetLastError();
    if(err != cudaSuccess)
    {
        printf("CUDA SLFC Executor Error: %s\n", cudaGetErrorString(err));
        return EXIT_FAILURE;
    }
    else
        return EXIT_SUCCESS;
}

int SILU::ELMRExecutor(device_data_csr *device_csr_matrix, val_type *d_x, 
                       val_type *d_b, char *ready)
{

    int num_threads = WARP_PER_BLOCK * WARP_SIZE;
    int num_blocks = ceil ((double)device_csr_matrix->m / (double)(num_threads/WARP_SIZE));


    if(num_blocks >= 1)
    {
        ELMRKernel<<<num_blocks, num_threads>>>
                            (device_csr_matrix->n, d_b, d_x, device_csr_matrix->diag, device_csr_matrix->d_csrVal, 
                             device_csr_matrix->d_csrColIdx, device_csr_matrix->d_csrRowPtr,
                             device_csr_matrix->jlev, ready);
    }    
    cudaError_t err = cudaGetLastError();
    if(err != cudaSuccess)
    {
        printf("CUDA ELMR Executor Error: %s\n", cudaGetErrorString(err));
        return EXIT_FAILURE;
    }
    else
    {
        //printf("ELMR Kernel success\n");
        return EXIT_SUCCESS;
    }
    
    //return EXIT_SUCCESS;
}

int SILU::ELMCExecutor(device_data_csc *device_csc_matrix, val_type *d_x, 
                       val_type *d_b, ind_type *count)
{
    int num_threads = WARP_PER_BLOCK * WARP_SIZE;
    int num_blocks = ceil ((double)device_csc_matrix->m / (double)(num_threads/WARP_SIZE));

    if(num_blocks >= 1)
    {
        ELMCKernel<<<num_blocks, num_threads>>>
                    (device_csc_matrix->n, d_x, device_csc_matrix->diag, device_csc_matrix->d_cscVal,
                     device_csc_matrix->d_cscRowIdx, device_csc_matrix->d_cscColPtr,
                     device_csc_matrix->jlev, count);
    }
    cudaError_t err = cudaGetLastError();
    if(err != cudaSuccess)
    {
        printf("CUDA ELMC Executor Error: %s\n", cudaGetErrorString(err));
        return EXIT_FAILURE;
    }
    else
    {
        //printf("ELMC Kernel success\n");
        return EXIT_SUCCESS;
    }
}
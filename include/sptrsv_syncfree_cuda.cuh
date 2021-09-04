// sptrsv_syncfree_cuda.h
// Created: 12-12-2019

#include "common.h"
#include <cuda_runtime.h>

__global__
void SLFCKernel(sz_type n, val_type *x, val_type *val, ind_type *jb, ind_type *ib,
          val_type *d, ind_type *dp, ind_type *jlev)
{
   int wid = (blockIdx.x * blockDim.x + threadIdx.x)/WARP_SIZE;
   int wlane = threadIdx.x /WARP_SIZE;
   int lane = threadIdx.x & (WARP_SIZE - 1);
   volatile __shared__ double s_xi[BLOCKDIM/WARP_SIZE];
   
   if(wid >= n) return;
   int i = jlev[wid];
   int p1 = jb[i], q1 = jb[i+1];
   double ti, xi;

   if(lane == 0)
   {
       d[i] = 1.0;
       ti = 1.0/d[i];
       //printf("dp[%d]=%d, jlev[%d]=%d\n", i, dp[i], i, jlev[i]);
       while(((volatile int *)dp)[i] != 1);
       xi = x[i] * ti; s_xi[wlane] = xi;
       
       //printf("warp %d is here, jlev[%d]=%d, ti=%f\n", i, i, jlev[i], ti);
   }
   __syncwarp();   
   xi = s_xi[wlane];
   for(int j = p1 + lane; j < q1; j+=WARP_SIZE)
   {
       atomicAdd(&x[ib[j]], (double)(-xi*val[j]));
       __threadfence();       
       atomicSub(&dp[ib[j]], 1);
       //printf("Subtracting from dp[%d]=%d, %d, %d\n", ib[j], dp[ib[j]], p1, q1);
   }
   if(lane == 0) 
   {
      x[i] = xi;      
   }
}

__global__ 
void ELMRKernel(sz_type n, val_type *f, volatile val_type *x, val_type *d, val_type *val, 
          ind_type *colind, ind_type *rowptr,ind_type *jlev, volatile char *ready)
{
    int wid = (blockIdx.x * blockDim.x + threadIdx.x) >> 5;
    int lane = threadIdx.x & (WARP_SIZE - 1);
    
    if(wid >= n) return;
    int p = 0, q = 0, i = -1;
    val_type sum = 0.0, diag;
   
    if(lane < 2)
    {
        i = jlev[wid];
        p = rowptr[i + lane];        
    }    
    q = __shfl_sync(-1, p, 1)-1;
    p = __shfl_sync(-1, p, 0);

    if(lane == 0)
    {
        sum = f[i]; diag = d[i];
    }

    //printf("wid:%d, p:%d, q:%d\n", wid, p, q);
    for(p += lane; p < q; p += WARP_SIZE)
    {
        while(ready[colind[p]] == 0);
        sum -= val[p] * x[colind[p]];
    }

    //printf("After\n");
    #pragma unroll
    for(int d = WARP_SIZE/2; d > 0; d >>= 1)
        sum += __shfl_down_sync(-1, sum, d);

    if(lane == 0)
    {
        x[i] = sum / diag;
        //printf("wid: %d, x[%d]=%f, diag:%f, ready[%d]=%d\n", wid, wid, x[i], diag, wid, ready[wid]);
        __threadfence();
        ready[i] = 1;
    }
}

__global__
void ELMCKernel(int n, volatile val_type *x, val_type *d, val_type *val,
          ind_type *rowind, ind_type *colptr, ind_type *jlev, volatile ind_type *count)
{
    int wid = (blockIdx.x * blockDim.x + threadIdx.x) >> 5;
    int lane = threadIdx.x & (WARP_SIZE - 1);

    if(wid >= n) return;
    int p = 0, q = 0;
    val_type t = 0.0;
    int i = 0;

    if(lane < 2)
    {
        i = jlev[wid];
        p = colptr[i+lane];
    }

    q=__shfl_sync(-1,p,1);
    p=__shfl_sync(-1,p,0)+1;
    
    if(lane == 0)
    {
        t = 1.0 / d[i]; 
        //printf("waiting[%d]:%d\n", i, count[i]);       
        while(count[i] != 1);
        //printf("Done %d, p:%d, q:%d, count[%d]:%d\n", i, p, q, i, count[i]);
        x[i] = t = x[i] * t;
        //printf("x[%d]=%f\n", i, x[i]);
    }

    t = __shfl_sync(-1,t,0);
    
    for(p+=lane; p < q; p += WARP_SIZE)
    {
        atomicAdd((val_type *)&x[rowind[p]], -t*val[p]);
        //printf("x[%d]=%f, %d, p:%d\n", rowind[p], x[rowind[p]], threadIdx.x, p);
        __threadfence();
        atomicSub((ind_type *)&count[rowind[p]], 1);        
    }    
}

__global__
void sptrsv_syncfree_cuda_analyser(const ind_type   *d_cscRowIdx,
                                   const sz_type    m,
                                   const sz_type    nnz,
                                         int   *d_graphInDegree)
{
    const int global_id = blockIdx.x * blockDim.x + threadIdx.x; //get_global_id(0);
    if (global_id < nnz)
    {
        atomicAdd(&d_graphInDegree[d_cscRowIdx[global_id]], 1);
    }
}

__global__
void sptrsv_syncfree_cuda_executor(const ind_type* __restrict__        d_cscColPtr,
                                   const ind_type* __restrict__        d_cscRowIdx,
                                   const val_type* __restrict__ d_cscVal,
                                         int*                     d_graphInDegree,
                                         val_type*              d_left_sum,
                                   const sz_type                  m,
                                   const int                      substitution,
                                   const val_type* __restrict__ d_b,
                                         val_type*              d_x
                                         )
{
    const int global_id = blockIdx.x * blockDim.x + threadIdx.x;
    int global_x_id = global_id / WARP_SIZE;
    if (global_x_id >= m) return;

    // substitution is forward or backward
    global_x_id = substitution == SUBSTITUTION_FORWARD ? 
                  global_x_id : m - 1 - global_x_id;

    volatile __shared__ int s_graphInDegree[WARP_PER_BLOCK];
    volatile __shared__ val_type s_left_sum[WARP_PER_BLOCK];

    // Initialize
    const int local_warp_id = threadIdx.x / WARP_SIZE;
    const int lane_id = (WARP_SIZE - 1) & threadIdx.x;
    int starting_x = (global_id / (WARP_PER_BLOCK * WARP_SIZE)) * WARP_PER_BLOCK;
    starting_x = substitution == SUBSTITUTION_FORWARD ? 
                  starting_x : m - 1 - starting_x;
    
    // Prefetch
    const int pos = substitution == SUBSTITUTION_FORWARD ?
                    d_cscColPtr[global_x_id] : d_cscColPtr[global_x_id+1]-1;
    const val_type coef = (val_type)1 / d_cscVal[pos];
    //asm("prefetch.global.L2 [%0];"::"d"(d_cscVal[d_cscColPtr[global_x_id] + 1 + lane_id]));
    //asm("prefetch.global.L2 [%0];"::"r"(d_cscRowIdx[d_cscColPtr[global_x_id] + 1 + lane_id]));

    if (threadIdx.x < WARP_PER_BLOCK) { s_graphInDegree[threadIdx.x] = 1; s_left_sum[threadIdx.x] = 0; }
    __syncthreads();
    //printf("Thread %d with dependencies %d going into wait, s=%d\n", global_x_id, d_graphInDegree[global_x_id],s_graphInDegree[local_warp_id]);
    clock_t start;
    // Consumer
    do {
        start = clock();
        __syncwarp();
    }
    while ((volatile int *)s_graphInDegree[local_warp_id] != (volatile int *)d_graphInDegree[global_x_id]);
    //printf("Thread %d reached here\n", global_x_id);
    //// Consumer
    //int graphInDegree;
    //do {
    //    //bypass Tex cache and avoid other mem optimization by nvcc/ptxas
    //    asm("ld.global.u32 %0, [%1];" : "=r"(graphInDegree),"=r"(d_graphInDegree[global_x_id]) :: "memory"); 
    //}
    //while (s_graphInDegree[local_warp_id] != graphInDegree );

    val_type xi = d_left_sum[global_x_id] + s_left_sum[local_warp_id];
    xi = (d_b[global_x_id] - xi) * coef;

    // Producer
    const int start_ptr = substitution == SUBSTITUTION_FORWARD ? 
                          d_cscColPtr[global_x_id]+1 : d_cscColPtr[global_x_id];
    const int stop_ptr  = substitution == SUBSTITUTION_FORWARD ? 
                          d_cscColPtr[global_x_id+1] : d_cscColPtr[global_x_id+1]-1;
    for (int jj = start_ptr + lane_id; jj < stop_ptr; jj += WARP_SIZE)
    {
        const int j = substitution == SUBSTITUTION_FORWARD ? jj : stop_ptr - 1 - (jj - start_ptr);
        const int rowIdx = d_cscRowIdx[j];
        const bool cond = substitution == SUBSTITUTION_FORWARD ? 
                    (rowIdx < starting_x + WARP_PER_BLOCK) : (rowIdx > starting_x - WARP_PER_BLOCK);
        if (cond) {
            const int pos = substitution == SUBSTITUTION_FORWARD ? 
                            rowIdx - starting_x : starting_x - rowIdx;
            atomicAdd((val_type *)&s_left_sum[pos], xi * d_cscVal[j]);
            //__threadfence_block();
            atomicAdd((int *)&s_graphInDegree[pos], 1);
        }
        else {
            atomicAdd(&d_left_sum[rowIdx], xi * d_cscVal[j]);
            //__threadfence();
            atomicSub(&d_graphInDegree[rowIdx], 1);
        }
    }

    //finish
    if (!lane_id) d_x[global_x_id] = xi;
}


// policy.h
// Level set calculations and data structures
// Created: 09-12-2019
// Author: Najeeb Ahmad

#ifndef __HEADER_POLICY__
#define __HEADER_POLICY__

#include "common.h"

#define DAG_SPLIT      1
#define MATRIX_SPLIT   2
#define SYNC_FREE      3
#define PARK_PARALLEL  4
#define MKL_SEQ_CSR    5
#define SPMV_CSR       6
#define SPMV_CSC       7
#define SPMV_CPU       8
#define SPMV_GPU       9
#define CUSPARSE_V2    10
#define MKL            11
#define SPTRSV_GPU     12
#define SPTRSV_CPU     13
#define SLFC           14
#define ELMR           15
#define ELMC           16
#define PARK           17   // For future extension
#define HTS            18 

#define MIN_ROWS                        9999
#define GPU_FRIENDLY_ROWS_MIN_PERC      10
#define CPU_FRIENDLY_LEVELS_MIN_PERC    50
#define CPU_FRIENDLY_MAX_ROWS_THRESHOLD 200
#define MIN_LEVELS                      40
#define MAX_LEVELS                      500
#define MIN_LVL_PERC_FOR_CPU_EXEC       90
#define MAX_GPU_FRIENDLY_PERC_FOR_CPU_EXEC 50
#define MAX_COL_LENGTH_THRESHOLD  30.0
#define MAX_ROW_LENGTH_THRESHOLD  10.0
#define MAX_CPU_FRIENDLY_LVLS_IN_SPLIT_PERC 10.0
#define MIN_GPU_FRIENDLY_ROWS_IN_SPLIT_PERC 85
#define WARP_COLWISE_ITERATIONS_THRESHOLD 100
#define WARP_ROWISE_ITERATIONS_THRESHOLD 75  

struct execution_policy
{
        int splitting_policy;       // Splitting DAG or matrix
	int sptrsv_method1;         // sptrsv method for top triangle 
	int sptrsv_method2;         // sprtsv method for bottom triangle
	int spmv_method;            // SpMV method ()
	int spmv_platform;          // SpMV on CPU or GPU
	int sptrsv1_platform;
	int sptrsv2_platform; 
	int split_at_level;
	int split_at_row;
        short description;   // 0=CPU, 1=GPU, 2=CPU-GPU, 3=GPU-CPU, 4=GPU-GPU, 5=CPU-CPU 
	int cpu_friendly_lvls;
	int gpu_friendly_lvls;
	int gpu_friendly_rows;
	float gpu_friendly_row_perc;
	float cpu_friendly_lvl_perc;
	int rhs;
	int nosplit;                // TRUE or FALSE. Whether split or not
};











#endif

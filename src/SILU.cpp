// SILU.cpp
// Implementation of Split Execution Framework Class
// Created: 03-12-2019
// Author: Najeeb Ahmad

#include "SILU.h"
#include "dummyFactorize.h"
#include "sptrsv_syncfree_serialref.h"
#include "quicksort.h"
#include <cuda_runtime.h>
#include <cusparse_v2.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
using namespace std::chrono;

#define CUDA_RT_CALL(call)                                                                  \
    {                                                                                       \
        cudaError_t cudaStatus = call;                                                      \
        if (cudaSuccess != cudaStatus)                                                      \
            fprintf(stderr,                                                                 \
                    "ERROR: CUDA RT call \"%s\" in line %d of file %s failed "              \
                    "with "                                                                 \
                    "%s (%d).\n",                                                           \
                    #call, __LINE__, __FILE__, cudaGetErrorString(cudaStatus), cudaStatus); \
    }

#define CUSPARSE_RT_CALL(call)                                                                  \
    {                                                                                       \
        cusparseStatus_t cusparseStatus = call;                                                      \
        if (cudaSuccess != cusparseStatus)                                                      \
            fprintf(stderr,                                                                 \
                    "ERROR: CUSPARSE RT call \"%s\" in line %d of file %s failed "              \
                    "with "                                                                 \
                    "%s (%d).\n",                                                           \
                    #call, __LINE__, __FILE__, cusparseGetErrorString(cusparseStatus), cusparseStatus); \
    }


SILU::SILU()
{
    GPU_device_id = 0;
    cusparseHandle = 0;
    descr = 0;
    descrL = 0;
    descrL1 = 0;
    ELMC_only=0;
    MKL_only=0;
    CUSPARSE_RT_CALL(cusparseCreate(&cusparseHandle));

    //cusparseMatDescr_t descr = 0;
    CUSPARSE_RT_CALL(cusparseCreateMatDescr(&descr));
    CUSPARSE_RT_CALL(cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL));
    cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);

    CUSPARSE_RT_CALL(cusparseCreateMatDescr(&descrL));
    CUSPARSE_RT_CALL(cusparseSetMatType(descrL,CUSPARSE_MATRIX_TYPE_GENERAL));
    CUSPARSE_RT_CALL(cusparseSetMatIndexBase(descrL,CUSPARSE_INDEX_BASE_ZERO));
    CUSPARSE_RT_CALL(cusparseSetMatFillMode(descrL, CUSPARSE_FILL_MODE_LOWER));
    CUSPARSE_RT_CALL(cusparseSetMatDiagType(descrL, CUSPARSE_DIAG_TYPE_UNIT));

    CUSPARSE_RT_CALL(cusparseCreateMatDescr(&descrL1));
    CUSPARSE_RT_CALL(cusparseSetMatType(descrL1,CUSPARSE_MATRIX_TYPE_GENERAL));
    CUSPARSE_RT_CALL(cusparseSetMatIndexBase(descrL1,CUSPARSE_INDEX_BASE_ZERO));
    CUSPARSE_RT_CALL(cusparseSetMatFillMode(descrL1, CUSPARSE_FILL_MODE_LOWER));
    CUSPARSE_RT_CALL(cusparseSetMatDiagType(descrL1, CUSPARSE_DIAG_TYPE_UNIT));


    pBuffer = NULL;
    pBuffer1 = NULL;
    infoA = 0;
    infoA1 = 0;
    CUSPARSE_RT_CALL(cusparseCreateCsrsv2Info(&infoA));
    CUSPARSE_RT_CALL(cusparseCreateCsrsv2Info(&infoA1));

    doubleone = 1.0;
    minusdoubleone = -1.0;
    doublezero = 0.0;
    m = 0;
    n = 0;

    // MKL
    ipar = new MKL_INT[128];
    dpar = new double[128];
    MKLdescrA.type = SPARSE_MATRIX_TYPE_GENERAL;
    MKLdescrA.mode = SPARSE_FILL_MODE_UPPER; 
    MKLdescrA.diag = SPARSE_DIAG_NON_UNIT;
    MKLdescrL.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
    MKLdescrL.mode = SPARSE_FILL_MODE_LOWER;
    MKLdescrL.diag = SPARSE_DIAG_UNIT;
    MKLdescrL1.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
    MKLdescrL1.mode = SPARSE_FILL_MODE_LOWER;
    MKLdescrL1.diag = SPARSE_DIAG_UNIT;

    ipar[30] = 1;
    dpar[30] = 1.E-32;
    dpar[31] = 1.E-32;

    time_sptrsv1 = 0.0;
    time_sptrsv2 = 0.0;
    time_spmv = 0.0;
    time_d2h = 0.0;
    time_h2d = 0.0;
    time_d2h_final = 0.0;    
    time_levels = 0.0;
    time_dag_analysis = 0.0;
    time_split = 0.0;
    time_data_structures = 0.0;
    time_algo_analysis = 0.0;
    time_factorize = 0.0;

    time_d2h_sptrsv1 = 0.0;
    time_d2h_spmv = 0.0;
    time_d2h_final = 0.0;
    time_h2d = 0.0;
    time_h2d_sptrsv1 = 0.0;
    time_h2d_spmv = 0.0;

    timing = 1;
    policy.splitting_policy = DAG_SPLIT; //MATRIX_SPLIT,DAG_SPLIT
    if(timing==1)
    {
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
    }
    init_ptrs();    
}

SILU::~SILU()
{	
    free_memory();
}

int SILU::SaveMatrixFeatures(csr_matrix *tri_mat, std::string filename)
{
  std::string file="matrix_features/"+filename+".bin";
  //char file[100];
  //sprintf(file, "%s/%s%s", "matrix_features",filename,".txt");
  //printf("%s",filename);
  std::cout << "FILE:" << file << std::endl;
#ifdef STORE_BINARY
    std::ofstream outfile(file, std::ios::out | std::ios::binary);
#else
    std::ofstream outfile(file, std::ios::out );
#endif

    if(outfile.is_open())
    {
    }
    else
    {        
        return EXIT_FAILURE;
    }
    //printf("Reached hererrr:%d\n", tri_mat->m);
#ifdef STORE_BINARY
#else
    outfile << "row_no,tot_levels,row_level,level_delta,parallelism,row_nnzs,row_nnzs,col_nnzs,avg_row_nnzs,avg_col_nnzs,avg_level_delta,col_center,col_nnzs\n";
#endif
    for(long i = 0; i < tri_mat->m; i++)
    {
#ifdef STORE_BINARY
    	long row_no = (i+1);
    	outfile.write((char *)&row_no, sizeof(long));
    	outfile.write((char *)&levels.nlevel, sizeof(int));
    	outfile.write((char *)&levels.level_of_row[i], sizeof(ind_type));
    	outfile.write((char *)&levels.level_delta[i], sizeof(ind_type));
        outfile.write((char *)&levels.parallelism[i], sizeof(float));
    	outfile.write((char *)&levels.nnz_per_row[i], sizeof(ind_type));
    	outfile.write((char *)&levels.avg_row_nnzs[i], sizeof(float));
    	outfile.write((char *)&levels.avg_col_nnzs[i], sizeof(float));
    	outfile.write((char *)&levels.avg_level_delta[i], sizeof(float));
    	outfile.write((char *)&levels.col_center[i], sizeof(float));
    	outfile.write((char *)&levels.nnz_per_col[i], sizeof(ind_type));
#else
        outfile << i+1 << ",";
        outfile << levels.nlevel << ","; 
        outfile << levels.level_of_row[i] << ",";
        outfile << levels.level_delta[i] << ",";
        outfile << std::fixed << std::setprecision(1) << levels.parallelism[i] << ",";
        outfile << levels.nnz_per_row[i] << ",";
        outfile << levels.avg_row_nnzs[i] << ",";
        outfile << levels.avg_col_nnzs[i] << ",";
        outfile << levels.avg_level_delta[i] << ",";
        outfile << levels.col_center[i] << ",";
        outfile << levels.nnz_per_col[i] << "\n";        
#endif
    }

    outfile.close();    
}


int SILU::SavePolicy(std::string filename, int mat_id)
{
  char file[50];
  sprintf(file,"%s/%s%s","policy",filename, ".txt");
  std::ofstream outfile(file, std::ios::out | std::ios::app);
  if(outfile.is_open())
    {
    }
  else
    {
      printf("Failure to open file to save\n");
      return EXIT_FAILURE;
    }
  outfile << mat_id << ",";
  outfile << policy.description << ",";
  outfile << policy.sptrsv_method1 << ",";
  outfile << policy.sptrsv_method2 << ",";
  outfile << policy.spmv_method << ",";
  outfile << policy.split_at_level << ",";
  outfile << levels.nlevel-1 << "\n";
  return 0;
}
int SILU::SaveLevelFeatures(int mat_id)
{
  char file[50];
  sprintf(file, "%s/%d%s", "level_info",mat_id,".txt");

  printf("%s\n", file);
#ifdef STORE_BINARY
    std::ofstream outfile(file, std::ios::out | std::ios::binary);
#else
    std::ofstream outfile(file, std::ios::out );
#endif

    if(outfile.is_open())
    {
    }
    else
    {
      printf("Failure to open file to save\n");
        return EXIT_FAILURE;
    }

    for(long i = 0; i < levels.nlevel; i++)
    {
 #ifdef STORE_BINARY
        int row_no = (i+1);
        outfile.write((char *)&row_no, sizeof(ind_type));
        outfile.write((char *)&levels.num_levelRows[i], sizeof(ind_type));
        outfile.write((char *)&levels.num_levelNnzs[i], sizeof(ind_type));
        outfile.write((char *)&levels.num_indegree[i], sizeof(ind_type));
        outfile.write((char *)&levels.num_outdegree[i], sizeof(ind_type));
	outfile.write((char *)&levels,max_levelRowLength[i], sizeof(ind_type));
	outfile.write((char *)&levels,max_levelColLength[i], sizeof(ind_type));

 #else       
        outfile << i+1 << ","; 
        outfile << levels.num_levelRows[i] << ",";
        outfile << levels.num_levelNnzs[i] << ",";
        //outfile << levels.num_indegree[i] << ",";
	outfile << levels.max_levelRowLength[i] << ",";
	outfile << levels.max_levelColLength[i] << "\n";
        //outfile << levels.num_outdegree[i] << "\n";
#endif
    }

    outfile.close();
}

void SILU::init_ptrs()
{
    levels.levelPtr=NULL; 
    levels.levelItem=NULL;
    levels.levelItemNewRowIdx=NULL;
    levels.num_levelRows=NULL;
    levels.num_levelNnzs=NULL;
    levels.num_indegree=NULL;
    levels.num_outdegree=NULL;
    levels.max_levelRowLength=NULL;
    levels.max_levelColLength=NULL;
    levels.cum_Rows=NULL;
    levels.cum_Nnzs=NULL;           
    TopTri_csc.cscRowIdx=NULL; 
    TopTri_csc.cscColPtr=NULL; 
    TopTri_csc.cscVal=NULL;
    TopTri_csr.csrRowPtr=NULL; 
    TopTri_csr.csrColIdx=NULL; 
    TopTri_csr.csrVal=NULL;    
    BotTri_csc.cscRowIdx=NULL; 
    BotTri_csc.cscColPtr=NULL; 
    BotTri_csc.cscVal=NULL; 
    BotTri_csr.csrRowPtr=NULL; 
    BotTri_csr.csrColIdx=NULL; 
    BotTri_csr.csrVal=NULL;
    Mvp_csr.csrRowPtr=NULL; 
    Mvp_csr.csrColIdx=NULL; 
    Mvp_csr.csrVal=NULL;
    d_graphInDegree_1=NULL;
    device_csc_matrix_top.d_cscColPtr=NULL;
    device_csc_matrix_top.d_cscRowIdx=NULL;
    device_csc_matrix_top.d_cscVal=NULL;          
    device_csc_matrix_top.d_left_sum=NULL;
    device_csc_matrix_top.diag=NULL;
    device_csc_matrix_top.jlev=NULL;        
    device_csr_matrix_top.d_csrRowPtr=NULL;
    device_csr_matrix_top.d_csrColIdx=NULL;
    device_csr_matrix_top.d_csrVal=NULL;
    d_b=NULL;      
    device_csr_matrix_top.diag=NULL;
    device_csr_matrix_top.jlev=NULL;
    ready_1=NULL;
    d_x=NULL;
    device_csr_matrix_mvp.d_csrRowPtr=NULL;
    device_csr_matrix_mvp.d_csrColIdx=NULL;
    device_csr_matrix_mvp.d_csrVal=NULL;      
    device_csc_matrix_bottom.d_cscColPtr=NULL;
    device_csc_matrix_bottom.d_cscRowIdx=NULL;
    device_csc_matrix_bottom.d_cscVal=NULL;
    d_x_1=NULL;
    device_csc_matrix_bottom.d_left_sum=NULL;
    d_graphInDegree_2=NULL;   
    device_csc_matrix_bottom.diag=NULL;
    device_csc_matrix_bottom.jlev=NULL;        
    device_csr_matrix_bottom.d_csrVal=NULL;
    device_csr_matrix_bottom.d_csrColIdx=NULL;
    device_csr_matrix_bottom.d_csrRowPtr=NULL;
    d_b_1=NULL;    
    device_csr_matrix_bottom.diag=NULL;
    device_csr_matrix_bottom.jlev=NULL;
    ready_2=NULL;           
    pBuffer1=NULL;
    pBuffer=NULL;   
}

int SILU::free_memory()
{
        free(levels.levelPtr); 
    free(levels.levelItem);
    free(levels.levelItemNewRowIdx);
    free(levels.num_levelRows);
    free(levels.num_levelNnzs);
    free(levels.num_indegree);
    free(levels.num_outdegree);
    free(levels.max_levelRowLength);
    free(levels.max_levelColLength);
    free(levels.cum_Rows);
    free(levels.cum_Nnzs);

    if(policy.sptrsv_method1 == SYNC_FREE || policy.sptrsv_method1 == SLFC || 
       policy.sptrsv_method1 == ELMC || policy.sptrsv_method1 == ELMR)
    {       
       free(TopTri_csc.cscRowIdx); 
       free(TopTri_csc.cscColPtr); 
       free(TopTri_csc.cscVal);
    }
    else
    {
       free(TopTri_csr.csrRowPtr); 
       free(TopTri_csr.csrColIdx); 
       free(TopTri_csr.csrVal); 
    }

    if(policy.sptrsv_method2 == SYNC_FREE  || policy.sptrsv_method2 == SLFC ||
       policy.sptrsv_method2 == ELMC || policy.sptrsv_method2 == ELMR)
    {
       free(BotTri_csc.cscRowIdx); 
       free(BotTri_csc.cscColPtr); 
       free(BotTri_csc.cscVal); 
    }
    else
    {
       free(BotTri_csr.csrRowPtr); 
       free(BotTri_csr.csrColIdx); 
       free(BotTri_csr.csrVal); 
    }

    free(Mvp_csr.csrRowPtr); 
    free(Mvp_csr.csrColIdx); 
    free(Mvp_csr.csrVal);

    
    if(policy.sptrsv_method1 == SYNC_FREE || policy.sptrsv_method1 == SLFC || policy.sptrsv_method1 == ELMC)
    {
        cudaFree(d_graphInDegree_1);
        cudaFree(device_csc_matrix_top.d_cscColPtr);
        cudaFree(device_csc_matrix_top.d_cscRowIdx);
        cudaFree(device_csc_matrix_top.d_cscVal);
           
        if(policy.sptrsv_method1 == SYNC_FREE)              
        {
            cudaFree(device_csc_matrix_top.d_left_sum);
            cudaFree(d_b);        
        }  
        else
        {        
            cudaFree(device_csc_matrix_top.diag);
            cudaFree(device_csc_matrix_top.jlev);        
        }
    }
    if(policy.sptrsv_method1 == CUSPARSE_V2 || policy.sptrsv_method1 == ELMR)
    {
        cudaFree(device_csr_matrix_top.d_csrRowPtr);
        cudaFree(device_csr_matrix_top.d_csrColIdx);
        cudaFree(device_csr_matrix_top.d_csrVal);
        cudaFree(d_b);
      
        if(policy.sptrsv_method1 == ELMR)
        {
            cudaFree(device_csr_matrix_top.diag);
            cudaFree(device_csr_matrix_top.jlev);
            cudaFree(ready_1);        
        }
    }

    if(policy.spmv_method = SPMV_CSR)
    {
        cudaFree(d_x);
        cudaFree(d_b_1);                                    
    }

    if(policy.spmv_method == SPMV_CSR && policy.spmv_platform == SPMV_GPU)
    {
        cudaFree(device_csr_matrix_mvp.d_csrRowPtr);
        cudaFree(device_csr_matrix_mvp.d_csrColIdx);
        cudaFree(device_csr_matrix_mvp.d_csrVal);      
    }
    
    if(policy.sptrsv_method2 == SYNC_FREE || policy.sptrsv_method2 == SLFC || policy.sptrsv_method2 == ELMC)
    {
        
        cudaFree(device_csc_matrix_bottom.d_cscColPtr);
        cudaFree(device_csc_matrix_bottom.d_cscRowIdx);
        cudaFree(device_csc_matrix_bottom.d_cscVal);
        cudaFree(d_x_1);
        cudaFree(device_csc_matrix_bottom.d_left_sum);
        cudaFree(d_graphInDegree_2);   

        if(policy.sptrsv_method2 == SLFC || policy.sptrsv_method2 == ELMC)
        {
            cudaFree(device_csc_matrix_bottom.diag);
            cudaFree(device_csc_matrix_bottom.jlev);        
        }        
    }

    if(policy.sptrsv_method2 == CUSPARSE_V2 || policy.sptrsv_method2 == ELMR)
      {
        cudaFree(device_csr_matrix_bottom.d_csrVal);
        cudaFree(device_csr_matrix_bottom.d_csrColIdx);
        cudaFree(device_csr_matrix_bottom.d_csrRowPtr);
        cudaFree(d_b_1);
        cudaFree(d_x_1);
        if(policy.sptrsv_method2 == ELMR)
        {
           cudaFree(device_csr_matrix_bottom.diag);
           cudaFree(device_csr_matrix_bottom.jlev);
           cudaFree(ready_2);           
        }
    }
        
    if(policy.sptrsv_method1 == CUSPARSE_V2)
    {
        cudaFree(pBuffer1);
    }
    if(policy.sptrsv_method2 == CUSPARSE_V2)
    {
        cudaFree(pBuffer);   
    }
}

int SILU::free_memory_1()
{
    
    if(levels.levelPtr!=NULL)
    {
        free(levels.levelPtr); levels.levelPtr=NULL;     
    }
    if(levels.levelItem!=NULL)
    {
        free(levels.levelItem); levels.levelItem=NULL;
    }    
    if(levels.levelItemNewRowIdx!=NULL)
    {
        free(levels.levelItemNewRowIdx); levels.levelItemNewRowIdx=NULL;
    }
    
    if(levels.num_levelRows!=NULL)
    {
        free(levels.num_levelRows); levels.num_levelRows=NULL;
    }
    
    if(levels.num_levelNnzs!=NULL)
    {
        free(levels.num_levelNnzs); levels.num_levelNnzs=NULL;
    }
    
    if(levels.num_indegree!=NULL)
    {
        free(levels.num_indegree); levels.num_indegree=NULL;
    }
    
    if(levels.num_outdegree!=NULL)
    {
        free(levels.num_outdegree); levels.num_outdegree=NULL;
    }
    
    if(levels.max_levelRowLength!=NULL)
    {
        free(levels.max_levelRowLength); levels.max_levelRowLength=NULL;
    }

    if(levels.max_levelColLength!=NULL)
    {
        free(levels.max_levelColLength); levels.max_levelColLength=NULL;
    }
    
    if(levels.cum_Rows!=NULL)
    {
        free(levels.cum_Rows); levels.cum_Rows=NULL;
    }
    
    if(levels.cum_Nnzs!=NULL)
    {
        free(levels.cum_Nnzs); levels.cum_Nnzs=NULL;
    }
    
    if(TopTri_csc.cscRowIdx!=NULL)
    {
        free(TopTri_csc.cscRowIdx); TopTri_csc.cscRowIdx=NULL;
    }       
    
    if(TopTri_csc.cscColPtr!=NULL)
    {
        free(TopTri_csc.cscColPtr); TopTri_csc.cscColPtr=NULL;
    }
    
    if(TopTri_csc.cscVal!=NULL)
    {
        free(TopTri_csc.cscVal); TopTri_csc.cscVal=NULL;    
    }
    
    if(TopTri_csr.csrRowPtr!=NULL)
    {
        free(TopTri_csr.csrRowPtr); TopTri_csr.csrRowPtr=NULL;     
    }
    if(TopTri_csr.csrColIdx!=NULL)
    {
        free(TopTri_csr.csrColIdx); TopTri_csr.csrColIdx=NULL;
    }
    
    if(TopTri_csr.csrVal!=NULL)
    {
        free(TopTri_csr.csrVal); TopTri_csr.csrVal=NULL;     
    }
    
    if(BotTri_csc.cscRowIdx!=NULL)
    {
        free(BotTri_csc.cscRowIdx); BotTri_csc.cscRowIdx=NULL;     
    }
    
    if(BotTri_csc.cscColPtr!=NULL)
    {
        free(BotTri_csc.cscColPtr); BotTri_csc.cscColPtr=NULL;    
    }

    if(BotTri_csc.cscVal!=NULL)
    {
        free(BotTri_csc.cscVal); BotTri_csc.cscVal=NULL;     
    }
    
    if(BotTri_csr.csrRowPtr!=NULL)
    {
        free(BotTri_csr.csrRowPtr); BotTri_csr.csrRowPtr=NULL;    
    }
    
    if(BotTri_csr.csrColIdx!=NULL)
    {
        free(BotTri_csr.csrColIdx); BotTri_csr.csrColIdx=NULL;    
    }
    
    if(BotTri_csr.csrVal!=NULL)
    {
        free(BotTri_csr.csrVal); BotTri_csr.csrVal=NULL;    
    }
    if(Mvp_csr.csrColIdx)
    {
        free(Mvp_csr.csrColIdx); Mvp_csr.csrColIdx=NULL;
    }
    if(Mvp_csr.csrVal!=NULL)
    {
        free(Mvp_csr.csrVal); Mvp_csr.csrVal=NULL;   
    }
    
    if(d_graphInDegree_1!=NULL)
    {
        cudaFree(d_graphInDegree_1); d_graphInDegree_1=NULL;   
    }
    
    if(device_csc_matrix_top.d_cscColPtr!=NULL)
    {
        cudaFree(device_csc_matrix_top.d_cscColPtr);  device_csc_matrix_top.d_cscColPtr=NULL;  
    }
    
    if(device_csc_matrix_top.d_cscRowIdx!=NULL)
    {
        cudaFree(device_csc_matrix_top.d_cscRowIdx); device_csc_matrix_top.d_cscRowIdx=NULL;   
    }
    
    if(device_csc_matrix_top.d_cscVal!=NULL)
    {
        cudaFree(device_csc_matrix_top.d_cscVal); device_csc_matrix_top.d_cscVal=NULL;   
    }
    
    if(device_csc_matrix_top.d_left_sum!=NULL)
    {
        cudaFree(device_csc_matrix_top.d_left_sum); device_csc_matrix_top.d_left_sum=NULL;   
    }       
    
    if(d_b!=NULL)
    {
        cudaFree(d_b); d_b=NULL;           
    }
    
    if(device_csc_matrix_top.diag!=NULL)
    {
        cudaFree(device_csc_matrix_top.diag); device_csc_matrix_top.diag=NULL;   
    }
    
    if(device_csc_matrix_top.jlev!=NULL)
    {
        cudaFree(device_csc_matrix_top.jlev); device_csc_matrix_top.jlev=NULL;           
    }
    
    if(device_csr_matrix_top.d_csrRowPtr!=NULL)
    {
        cudaFree(device_csr_matrix_top.d_csrRowPtr); device_csr_matrix_top.d_csrRowPtr=NULL;   
    }
    
    if(device_csr_matrix_top.d_csrColIdx!=NULL)
    {
        cudaFree(device_csr_matrix_top.d_csrColIdx); device_csr_matrix_top.d_csrColIdx=NULL;   
    }
    
    if(device_csr_matrix_top.d_csrVal!=NULL)
    {
        cudaFree(device_csr_matrix_top.d_csrVal); device_csr_matrix_top.d_csrVal=NULL;   
    }
    
    if(device_csr_matrix_top.diag!=NULL)
    {
        cudaFree(device_csr_matrix_top.diag); device_csr_matrix_top.diag=NULL;   
    }  
    
    if(device_csr_matrix_top.jlev!=NULL)
    {
        cudaFree(device_csr_matrix_top.jlev); device_csr_matrix_top.jlev=NULL;   
    }
    
    if(ready_1!=NULL)
    {
        cudaFree(ready_1); ready_1=NULL;   
    }
    if(d_x!=NULL)
    {
        cudaFree(d_x); d_x=NULL;   
    }
    
    if(d_b_1!=NULL)
    {
        cudaFree(d_b_1); d_b_1==NULL;    
    }
    if(device_csr_matrix_mvp.d_csrRowPtr!=NULL)
    {
        cudaFree(device_csr_matrix_mvp.d_csrRowPtr); device_csr_matrix_mvp.d_csrRowPtr=NULL;
    }    
    if(device_csr_matrix_mvp.d_csrColIdx!=NULL)
    {
        cudaFree(device_csr_matrix_mvp.d_csrColIdx); device_csr_matrix_mvp.d_csrColIdx=NULL;
    }    
    if(device_csr_matrix_mvp.d_csrVal!=NULL)
    {
        cudaFree(device_csr_matrix_mvp.d_csrVal); device_csr_matrix_mvp.d_csrVal=NULL;     
    }    
    
    if(device_csc_matrix_bottom.d_cscColPtr!=NULL)
    {
        cudaFree(device_csc_matrix_bottom.d_cscColPtr); device_csc_matrix_bottom.d_cscColPtr=NULL;
    }    
    if(device_csc_matrix_bottom.d_cscRowIdx!=NULL)
    {
        cudaFree(device_csc_matrix_bottom.d_cscRowIdx); device_csc_matrix_bottom.d_cscRowIdx=NULL;  
    }
    if(device_csc_matrix_bottom.d_cscVal!=NULL)
    {
        cudaFree(device_csc_matrix_bottom.d_cscVal); device_csc_matrix_bottom.d_cscVal=NULL;
    }   
    if(d_x_1!=NULL)
    {
        cudaFree(d_x_1); d_x_1=NULL;
    }   
    if(device_csc_matrix_bottom.d_left_sum!=NULL)
    {
        cudaFree(device_csc_matrix_bottom.d_left_sum); device_csc_matrix_bottom.d_left_sum=NULL;
    }   
    if(d_graphInDegree_2!=NULL)
    {
        cudaFree(d_graphInDegree_2); d_graphInDegree_2=NULL;
    }      
    if(device_csc_matrix_bottom.diag!=NULL)
    {
        cudaFree(device_csc_matrix_bottom.diag); device_csc_matrix_bottom.diag=NULL;
    }      
    if(device_csc_matrix_bottom.jlev!=NULL)
    {
        cudaFree(device_csc_matrix_bottom.jlev); device_csc_matrix_bottom.jlev=NULL;       
    }   
    
   if(device_csr_matrix_bottom.d_csrVal!=NULL)
    {
        cudaFree(device_csr_matrix_bottom.d_csrVal); device_csr_matrix_bottom.d_csrVal=NULL;
    }   
    if(device_csr_matrix_bottom.d_csrColIdx!=NULL)
    {
        cudaFree(device_csr_matrix_bottom.d_csrColIdx); device_csr_matrix_bottom.d_csrColIdx=NULL;
    }   
    if(device_csr_matrix_bottom.d_csrRowPtr!=NULL)
    {
        cudaFree(device_csr_matrix_bottom.d_csrRowPtr); device_csr_matrix_bottom.d_csrRowPtr=NULL;
    }  
    if(d_b_1!=NULL)
    {
        cudaFree(d_b_1); d_b_1=NULL;
    }        
    if(device_csr_matrix_bottom.diag!=NULL)
    {
        cudaFree(device_csr_matrix_bottom.diag); device_csr_matrix_bottom.diag=NULL;
    }   
    if(device_csr_matrix_bottom.jlev!=NULL)
    {
        cudaFree(device_csr_matrix_bottom.jlev); device_csr_matrix_bottom.jlev=NULL; 
    } 
    if(ready_2!=NULL)
    {
        cudaFree(ready_2); ready_2=NULL;          
    }   
    if(pBuffer1!=NULL)
    {
        cudaFree(pBuffer1); pBuffer1=NULL;
    }  
    
    if(pBuffer!=NULL)
    {
        cudaFree(pBuffer); pBuffer=NULL;  
    }    
    //printf("here\n");    
}

// Smart ILU factorization routine
// @return Factorization status
// @param  A  input csr matrix
// @param  L  output lower triangular hybrid matrix
// @param  U  output upper triangular hybrid matrix
// @param  method ILU0 factorization method
//         Options: dummy: Upper and lower parts of A used as ILU0 factors
//                  MKL: MKL ILU0 factorization is used
//                  cuSPARSE: cuSPARSE ILU0 factorization is used
int SILU::Factorize(csr_matrix *A, csr_matrix *L, csr_matrix *U, std::string method)
{
	start_profile_CPU();
    if(method == "dummy")
	{
		dummyFactorize(A, L, U, SUBSTITUTION_FORWARD);		
	}
	else if (method == "MKL")
	{

	}
	else if (method == "cuSPARSE")
	{

	}
	else
	{
		printf("Unknown factorization method specified\n");
		return EXIT_FAILURE;
	}
    end_profile_CPU(&time_factorize);

	return EXIT_SUCCESS;
}

// Smart ILU Analysis routine
// @return Analysis status
// @param  tri_mat is the triangular matrix. 
// @param  method  Whether to use DAG partitioning or
///                matrix partitioning
//         Options: level-set, matrix
// It is stored in TopTri in hybrid format
int SILU::Analyze(csr_matrix *mat, csr_matrix *tri_mat, val_type *b)
{
	if(policy.splitting_policy == DAG_SPLIT)
	{
		//printf("DAG split\n");
        start_profile_CPU();        
        findlevels_csr(tri_mat);
        end_profile_CPU(&time_levels);
		
	//printf("levels found\n");
        start_profile_CPU();
	//analyse_dag(tri_mat);
	analyse_dag_revised(tri_mat);
        end_profile_CPU(&time_dag_analysis);
	

	//return 0;
	//printf("dag analyzed\n");
        start_profile_CPU();
	split_dag(tri_mat);
        end_profile_CPU(&time_split);

	//printf("dag split\n");
        setup_execution_devices();	
        
        start_profile_CPU();		
	setup_device_datastructures(b);
        end_profile_CPU(&time_data_structures);

	//printf("setup data structures done\n");
	//printf("MY VALUE IS: %d\n", device_csr_matrix_top.nnz);
        start_profile_CPU();
        algorithm_analysis();
        end_profile_CPU(&time_algo_analysis);
	//printf("algorithm analysis done\n");
	}
	else
	{
        printf("Matrix split\n");
        start_profile_CPU();
        findlevels_csr(tri_mat);
        end_profile_CPU(&time_levels);
        
        start_profile_CPU();
        analyse_dag(tri_mat);
        end_profile_CPU(&time_dag_analysis);
        
        start_profile_CPU();
        split_matrix(mat, tri_mat);
        end_profile_CPU(&time_split);
        
        setup_execution_devices();
        
        start_profile_CPU();
        setup_device_datastructures(b);
        end_profile_CPU(&time_data_structures);

        start_profile_CPU();
        algorithm_analysis();
        end_profile_CPU(&time_algo_analysis);       

	}
	return EXIT_SUCCESS;
}

int SILU::setup_execution_devices(void)
{
	if(policy.sptrsv_method1 == SYNC_FREE || policy.sptrsv_method2 == SYNC_FREE ||
       policy.sptrsv_method1 == SLFC || policy.sptrsv_method2 == SLFC || 
       policy.sptrsv_method1 == ELMC || policy.sptrsv_method2 == ELMC ||
       policy.sptrsv_method1 == ELMR || policy.sptrsv_method2 == ELMR ||
       policy.sptrsv_method1 == CUSPARSE_V2 || policy.sptrsv_method2 == CUSPARSE_V2)
	{		
		printf("Querying GPU ...\n");
		cudaSetDevice(0);
    	cudaDeviceProp deviceProp;
    	cudaGetDeviceProperties(&deviceProp, 0);    	
    	printf("GPU [%i] %s @ %4.2f MHz\n", 0, deviceProp.name, deviceProp.clockRate * 1e-3f);		
	}
	return EXIT_SUCCESS;
}

int SILU::algorithm_analysis()
{
	if(policy.sptrsv_method1 == SYNC_FREE || policy.sptrsv_method1 == SLFC || policy.sptrsv_method1 == ELMC)
	{
		CUDA_RT_CALL(cudaMalloc((void **)&d_graphInDegree_1, (device_csc_matrix_top.m) * sizeof(int)));
        ind_type *temp = (ind_type *)malloc(sizeof(ind_type) * device_csc_matrix_top.m);
        for(int i = 0; i < device_csc_matrix_top.m; i++)
            temp[i] = 1.0;
        CUDA_RT_CALL(cudaMemcpy(d_graphInDegree_1, temp, device_csc_matrix_top.m * sizeof(ind_type), cudaMemcpyHostToDevice));
        free(temp);
        SyncFreeAnalyzer(&device_csc_matrix_top, d_graphInDegree_1);
        cudaDeviceSynchronize();                          
		
	}
    
    if(policy.sptrsv_method1 == CUSPARSE_V2)
    {        
        CUSPARSE_RT_CALL(cusparseDcsrsv2_bufferSize(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, device_csr_matrix_top.m,
                                                    device_csr_matrix_top.nnz, descrL1, device_csr_matrix_top.d_csrVal, 
                                                    device_csr_matrix_top.d_csrRowPtr,
                                                    device_csr_matrix_top.d_csrColIdx, infoA1, &pBufferSizeInBytesL1));
        

        //pBufferSize1 = pBufferSizeInBytesL1;
        CUDA_RT_CALL(cudaMalloc((void **)&pBuffer1, pBufferSizeInBytesL1));
        //printf("%d\n", pBufferSizeInBytesL1);
        //printf("%X:%X:%X||||%d:%d:%d\n",device_csr_matrix_top.d_csrVal, 
        //    device_csr_matrix_top.d_csrRowPtr,
        //    device_csr_matrix_top.d_csrColIdx, device_csr_matrix_top.m, device_csr_matrix_top.n, device_csr_matrix_top.nnz);
        
        
        CUSPARSE_RT_CALL(cusparseDcsrsv2_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, device_csr_matrix_top.m,
                                                  device_csr_matrix_top.nnz, descrL1, device_csr_matrix_top.d_csrVal,
                                                  device_csr_matrix_top.d_csrRowPtr, device_csr_matrix_top.d_csrColIdx, infoA1, 
                          CUSPARSE_SOLVE_POLICY_USE_LEVEL, pBuffer1));

    }

    if(policy.sptrsv_method1 == MKL)
    {        
        if(TopTri_csr.m > 0)
        {
            mkl_sparse_set_sv_hint(csrLU1, SPARSE_OPERATION_NON_TRANSPOSE, MKLdescrL1, n_iter_LU1);
            mkl_sparse_optimize(csrLU1);    
        }
    }

    if(policy.sptrsv_method2 == SYNC_FREE  || policy.sptrsv_method2 == SLFC || policy.sptrsv_method2 == ELMC)
    {
        cudaError_t cudaStatus;
        cudaStatus=cudaMalloc((void **)&d_graphInDegree_2, device_csc_matrix_bottom.m * sizeof(int));
        ind_type *temp = (ind_type *)malloc(sizeof(ind_type) * device_csc_matrix_bottom.m);
        for(int i = 0; i < device_csc_matrix_bottom.m; i++)
            temp[i] = 1.0;
        cudaMemcpy(d_graphInDegree_2, temp, device_csc_matrix_bottom.m * sizeof(ind_type), cudaMemcpyHostToDevice);
        free(temp);
        if(cudaStatus == 0)
        {
            SyncFreeAnalyzer(&device_csc_matrix_bottom, d_graphInDegree_2);            
        }        
        
    }
    
    if(policy.spmv_method == SPMV_CSR && policy.spmv_platform == SPMV_CPU)
    {        
        if(Mvp_csr.m > 0)
        {
            mkl_sparse_set_sv_hint(csrA, SPARSE_OPERATION_NON_TRANSPOSE, MKLdescrA, n_iter);
            mkl_sparse_optimize(csrA);
        }
    }

    if(policy.sptrsv_method2 == CUSPARSE_V2)
    {
        CUSPARSE_RT_CALL(cusparseDcsrsv2_bufferSize(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, device_csr_matrix_bottom.m,
                                                    device_csr_matrix_bottom.nnz, descrL, device_csr_matrix_bottom.d_csrVal, 
                                                    device_csr_matrix_bottom.d_csrRowPtr,
                                                    device_csr_matrix_bottom.d_csrColIdx, infoA, &pBufferSizeInBytesL));
        CUDA_RT_CALL(cudaMalloc((void **)&pBuffer, pBufferSizeInBytesL));

        CUSPARSE_RT_CALL(cusparseDcsrsv2_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, device_csr_matrix_bottom.m,
                                                  device_csr_matrix_bottom.nnz, descrL, device_csr_matrix_bottom.d_csrVal,
                                                  device_csr_matrix_bottom.d_csrRowPtr, device_csr_matrix_bottom.d_csrColIdx, infoA, 
						  CUSPARSE_SOLVE_POLICY_USE_LEVEL, pBuffer));        	
    }

    if(policy.sptrsv_method2 == MKL)
    {        
        if(BotTri_csr.m > 0)
        {
            mkl_sparse_set_sv_hint(csrLU, SPARSE_OPERATION_NON_TRANSPOSE, MKLdescrL, n_iter_LU);
            mkl_sparse_optimize(csrLU);                
        }
    }    
            
	return EXIT_SUCCESS;
}

int SILU::start_profile_CPU()
{
    if(timing == 1)
    {
        hrt1 = high_resolution_clock::now();
    }

    return EXIT_SUCCESS;
}
// int SILU::start_profile_CPU()
// {
//     if(timing == 1)
//         gettimeofday(&t1, NULL);
//     return EXIT_SUCCESS;
// }

int SILU::start_profile_GPU()
{
    if(timing == 1)
    {
        cudaEventRecord(start);
    }
    return EXIT_SUCCESS;
}

int SILU::end_profile_CPU(float *elapsed_time)
{
    if(timing == 1)
    {
        hrt2 = high_resolution_clock::now();
        duration<float> t_diff = hrt2 - hrt1;
        nanoseconds t_diff_ns = duration_cast<nanoseconds>(t_diff);
        *elapsed_time = t_diff_ns.count()/1e6;
    }
    return EXIT_SUCCESS;
}

// int SILU::end_profile_CPU(float *elapsed_time)
// {
//     if(timing == 1)
//     {       
//        gettimeofday(&t2, NULL);
//        *elapsed_time += (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;
//     }
//     return EXIT_SUCCESS;
// }


int SILU::end_profile_GPU(float *elapsed_time)
{
    if(timing == 1)
    {       
       float temp_time = 0.0;
       cudaEventRecord(stop);
       cudaEventSynchronize(stop);
       cudaEventElapsedTime(&temp_time, start, stop);

       *elapsed_time = temp_time;
    }
    return EXIT_SUCCESS;
}

 int SILU::HTSExecutor(device_data_csr *device_csr_matrix)
 {
   
 }

int SILU::trsv(val_type *rhs, val_type *x, int direction)
{
  if(direction == SUBSTITUTION_FORWARD)
	{
		//printf("Forward solving\n");
	
		if (policy.rhs == 1)
        {
            if(policy.sptrsv_method1 == SYNC_FREE)
            {                
                start_profile_GPU();
                SyncFreeExecutor(&device_csc_matrix_top, d_x, d_b, d_graphInDegree_1, direction);                
                end_profile_GPU(&time_sptrsv1);                    
                
                start_profile_GPU();
                cudaMemcpy(x, d_x, device_csc_matrix_top.m * policy.rhs * sizeof(val_type), cudaMemcpyDeviceToHost);
                end_profile_GPU(&time_d2h_sptrsv1);                

            }   

            if(policy.sptrsv_method1 == SLFC)
            {   
                start_profile_GPU();
                SLFCExecutor(&device_csc_matrix_top, d_x, device_csc_matrix_top.diag,
                device_csc_matrix_top.jlev, d_graphInDegree_1);
                end_profile_GPU(&time_sptrsv1);
                
                start_profile_GPU();
                cudaMemcpy(x, d_x, device_csc_matrix_top.m * policy.rhs * sizeof(val_type), cudaMemcpyDeviceToHost);             
                end_profile_GPU(&time_d2h_sptrsv1);
                
            }

            if(policy.sptrsv_method1 == ELMC)
            {
                start_profile_GPU();
                ELMCExecutor(&device_csc_matrix_top, d_x, d_b, d_graphInDegree_1);
                end_profile_GPU(&time_sptrsv1);

                start_profile_GPU();
                cudaMemcpy(x, d_x, device_csc_matrix_top.m * policy.rhs * sizeof(val_type), cudaMemcpyDeviceToHost);                
                end_profile_GPU(&time_d2h_sptrsv1);                
            }

            if(policy.sptrsv_method1 == CUSPARSE_V2)
            {
                start_profile_GPU();              
                cusparseStatus = cusparseDcsrsv2_solve(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, device_csr_matrix_top.m,
                               device_csr_matrix_top.nnz, &doubleone, descrL1,device_csr_matrix_top.d_csrVal,
                               device_csr_matrix_top.d_csrRowPtr, device_csr_matrix_top.d_csrColIdx, infoA1,
                               d_b, d_x, CUSPARSE_SOLVE_POLICY_USE_LEVEL, pBuffer1);
                
                end_profile_GPU(&time_sptrsv1);

                if(cusparseStatus != cudaSuccess)
                {
                    printf("Error performing cuSPARSE v2 triangular solve! Error code: %d\n", cusparseStatus);
                }

                start_profile_GPU();
                cudaMemcpy(x, d_x, device_csr_matrix_top.m * policy.rhs * sizeof(val_type), cudaMemcpyDeviceToHost);
                if(policy.spmv_platform == SPMV_CPU)
                    cudaDeviceSynchronize();                                
                end_profile_GPU(&time_d2h_sptrsv1);
            }
            if(policy.sptrsv_method1 == ELMR)
            {
                start_profile_GPU();
                ELMRExecutor(&device_csr_matrix_top, d_x, d_b, ready_1);
                end_profile_GPU(&time_sptrsv1);

                start_profile_GPU();
                cudaMemcpy(x, d_x, device_csr_matrix_top.m * policy.rhs * sizeof(val_type), cudaMemcpyDeviceToHost);
                if(policy.spmv_platform == SPMV_CPU)
                    cudaDeviceSynchronize();
                end_profile_GPU(&time_d2h_sptrsv1);                
            }

            if(policy.sptrsv_method1 == MKL)
            {
                cudaError_t cudaStatus;
                
                start_profile_CPU();                
                int stat = mkl_sparse_d_trsv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, 
                                              csrLU1, MKLdescrL1, rhs, x);
                if(stat != SPARSE_STATUS_SUCCESS)
                    printf("Error MKL TRSV\n");
                end_profile_CPU(&time_sptrsv1);
                
                if(policy.spmv_platform == SPMV_GPU)
                {
                    start_profile_GPU();
                    cudaStatus = cudaMemcpy(d_x, x, TopTri_csr.m * policy.rhs * sizeof(val_type), cudaMemcpyHostToDevice);
                    if(cudaStatus != cudaSuccess)
                    {
                        printf("Error copying\n");
                    }
                    end_profile_GPU(&time_h2d_sptrsv1);    
                }
            }
	    
	    if(policy.sptrsv_method1 == HTS)
	      {
		int stat = HTSExecutor(&device_csr_matrix_top);
	      }
            
            if(policy.spmv_method == SPMV_CSR && policy.spmv_platform == SPMV_GPU)
            {
                int status = 0;                
                start_profile_GPU();
                if(policy.sptrsv_method2 == SLFC || policy.sptrsv_method2 == ELMC)
                {
                    cusparseStatus = cusparseDcsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE, 
                                            device_csr_matrix_mvp.m, device_csr_matrix_mvp.n, 
                                            device_csr_matrix_mvp.nnz, &minusdoubleone, descrL, 
                                            device_csr_matrix_mvp.d_csrVal, device_csr_matrix_mvp.d_csrRowPtr, 
                                            device_csr_matrix_mvp.d_csrColIdx, d_x, 
                                            &doubleone, d_x_1);                   
                }
                else
                {
                    cusparseStatus = cusparseDcsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE, 
                                            device_csr_matrix_mvp.m, device_csr_matrix_mvp.n, 
                                            device_csr_matrix_mvp.nnz, &minusdoubleone, descrL, 
                                            device_csr_matrix_mvp.d_csrVal, device_csr_matrix_mvp.d_csrRowPtr, 
                                            device_csr_matrix_mvp.d_csrColIdx, d_x, 
                                            &doubleone, d_b_1);
                }
                end_profile_GPU(&time_spmv);

                if(cusparseStatus != cudaSuccess)
                {
                    printf("GPU SpMV failed\n");
                }                
                
                CUDA_RT_CALL(cudaDeviceSynchronize());
                start_profile_GPU();
                if(policy.sptrsv_method1 == SYNC_FREE || policy.sptrsv_method1 == SLFC || policy.sptrsv_method1 == ELMC)
                {
                    status = cudaMemcpy(rhs+device_csc_matrix_top.m, d_b_1, device_csr_matrix_mvp.m * policy.rhs * sizeof(val_type), cudaMemcpyDeviceToHost);
                }
                else
                {                
                    status = cudaMemcpy(rhs+device_csr_matrix_top.m, d_b_1, device_csr_matrix_mvp.m * policy.rhs * sizeof(val_type), cudaMemcpyDeviceToHost);                
                }
                end_profile_GPU(&time_d2h_spmv);
                int bytes = device_csr_matrix_mvp.m * 8;
                printf("Bytes:%d, Expected Time:%f, Actual Time:%f, BW: %f GB/sec\n", device_csr_matrix_mvp.m*8, (device_csr_matrix_mvp.m*8)/(32*1e6), time_d2h_spmv, bytes * 1e-6/time_d2h_spmv );
                if(status != cudaSuccess)
                {
                    printf("Error performing memory copy\n");
                }
                cudaDeviceSynchronize();                
            }           
            
            if(policy.spmv_method == SPMV_CSR && policy.spmv_platform == SPMV_CPU)
            {
                start_profile_CPU();
                if(policy.sptrsv_method1 == SYNC_FREE || policy.sptrsv_method1 == SLFC || policy.sptrsv_method1 == ELMC)
                {   
                    if(Mvp_csr.m > 0)
                    {
                        sp_status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, -1.0, csrA, MKLdescrA, x, 1.0, rhs+device_csc_matrix_top.m);
                        if(sp_status != SPARSE_STATUS_SUCCESS)
                        {
                            printf("Error doing MKL spMV, Error:%d\n", sp_status);
                        }
                    }
                }
                else
                {                
                    int status;
                    if(Mvp_csr.m > 0)
                    {
                         status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, -1.0, csrA, MKLdescrA, x, 1.0, rhs+device_csr_matrix_top.m);                    
                        if(status != SPARSE_STATUS_SUCCESS)
                        {
                            printf("Error doing MKL spMV 1, Error:%d\n", sp_status);
                        }
                    }   
                }
                end_profile_CPU(&time_spmv);

                start_profile_GPU();
                if(policy.sptrsv_method1 == SYNC_FREE || policy.sptrsv_method1 == SLFC || policy.sptrsv_method1 == ELMC)
                {
                    int status;
                    if(policy.sptrsv_method2 != MKL)
                    {
                        if(policy.sptrsv_method2 == SLFC || policy.sptrsv_method2 == ELMC)
                        {
                            status = cudaMemcpy(d_x_1, rhs+device_csc_matrix_top.m, Mvp_csr.m * policy.rhs * sizeof(val_type), cudaMemcpyHostToDevice);
                        }   
                       else
                           status = cudaMemcpy(d_b_1, rhs+device_csc_matrix_top.m, Mvp_csr.m * policy.rhs * sizeof(val_type), cudaMemcpyHostToDevice);                   

                        if(status != SPARSE_STATUS_SUCCESS)
                        {
                            printf("Memory copy error\n");
                        }
                    }
                }
                else
                {
                    int status;                  
                    if(policy.sptrsv_method2 != MKL)
                    {
                        if(policy.sptrsv_method2 == SLFC || policy.sptrsv_method2 == ELMC)
                            status = cudaMemcpy(d_x_1, rhs+device_csr_matrix_top.m, Mvp_csr.m * policy.rhs * sizeof(val_type), cudaMemcpyHostToDevice);
                        else
                            status = cudaMemcpy(d_b_1, rhs+device_csr_matrix_top.m, Mvp_csr.m * policy.rhs * sizeof(val_type), cudaMemcpyHostToDevice);                    

                        if(status != SPARSE_STATUS_SUCCESS)
                        {
                            printf("Memory copy error\n");
                        }    
                    }
                    
                }
                end_profile_GPU(&time_h2d_spmv);                
            }
            
          
          if(policy.sptrsv_method2 == SYNC_FREE)
	      {		        
                start_profile_GPU();
                SyncFreeExecutor(&device_csc_matrix_bottom, d_x_1, d_b_1, d_graphInDegree_2, direction);
                end_profile_GPU(&time_sptrsv2);

                start_profile_GPU();
                if(policy.sptrsv_method1 == SYNC_FREE || policy.sptrsv_method1 == SLFC || policy.sptrsv_method1 == ELMC)
                {
                    cudaMemcpy(x+device_csc_matrix_top.m, d_x_1, device_csc_matrix_bottom.m * policy.rhs * sizeof(val_type), cudaMemcpyDeviceToHost);    
                }
                else
                {  

                    cudaMemcpy(x+device_csr_matrix_top.m, d_x_1, device_csc_matrix_bottom.m * policy.rhs * sizeof(val_type), cudaMemcpyDeviceToHost);
                }
                end_profile_GPU(&time_d2h_final);                
	      }
	    
	    if(policy.sptrsv_method2 == CUSPARSE_V2)
	      {
		        start_profile_GPU();
                cusparseStatus = cusparseDcsrsv2_solve(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, device_csr_matrix_bottom.m,
						       device_csr_matrix_bottom.nnz, &doubleone, descrL,device_csr_matrix_bottom.d_csrVal,
						       device_csr_matrix_bottom.d_csrRowPtr, device_csr_matrix_bottom.d_csrColIdx, infoA,
						       d_b_1, d_x_1, CUSPARSE_SOLVE_POLICY_USE_LEVEL,
						       pBuffer);
                end_profile_GPU(&time_sptrsv2);

                if(cusparseStatus != cudaSuccess)
                {
                    printf("Error performing cuSPARSE v2 triangular solve! Error code: %d\n", cusparseStatus);
                }

                start_profile_GPU();
        		if(policy.sptrsv_method1 == SYNC_FREE || policy.sptrsv_method1 == SLFC || policy.sptrsv_method1 == ELMC)
        		{
                    cudaMemcpy(x+device_csc_matrix_top.m, d_x_1, 
                               device_csr_matrix_bottom.m * policy.rhs * sizeof(val_type), 
                               cudaMemcpyDeviceToHost);                    
                }  
        		else
                {                    
        		    cudaMemcpy(x+device_csr_matrix_top.m, d_x_1, 
                               device_csr_matrix_bottom.m * policy.rhs * sizeof(val_type), 
                               cudaMemcpyDeviceToHost);
                }
                end_profile_GPU(&time_d2h_final);
          }
        
        if(policy.sptrsv_method2 == MKL)
        {
            start_profile_CPU(); 
            if(policy.sptrsv_method1 == SYNC_FREE || policy.sptrsv_method1 == SLFC || policy.sptrsv_method1 == ELMC)
            {                
                if(device_csr_matrix_bottom.m > 0)               
                sp_status = mkl_sparse_d_trsv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, 
                                              csrLU, MKLdescrL, rhs+device_csc_matrix_top.m, x+device_csc_matrix_top.m);                    
            }
            else
            {
                if(device_csr_matrix_bottom.m > 0)                
                sp_status = mkl_sparse_d_trsv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, 
                                              csrLU, MKLdescrL, rhs+device_csr_matrix_top.m, x+device_csr_matrix_top.m);             
            }            
            end_profile_CPU(&time_sptrsv2); 
        }
        
        if(policy.sptrsv_method2 == SLFC)
        {            
            start_profile_GPU();            
            SLFCExecutor(&device_csc_matrix_bottom, d_x_1, device_csc_matrix_bottom.diag,
                        device_csc_matrix_bottom.jlev, d_graphInDegree_2);
            end_profile_GPU(&time_sptrsv2);
            
            start_profile_GPU();
            if(policy.sptrsv_method1 == SYNC_FREE || policy.sptrsv_method1 == SLFC || policy.sptrsv_method1 == ELMC)
            {                
                cudaMemcpy(x+device_csc_matrix_top.m, d_x_1, 
                           device_csc_matrix_bottom.m * policy.rhs * sizeof(val_type), 
                           cudaMemcpyDeviceToHost);                        
            }  
            else
            {   
                cudaMemcpy(x+device_csr_matrix_top.m, d_x_1, 
                           device_csc_matrix_bottom.m * policy.rhs * sizeof(val_type), 
                           cudaMemcpyDeviceToHost);                
            }
            end_profile_GPU(&time_d2h_final);                         
        }

        if(policy.sptrsv_method2 == ELMC)
        {
            start_profile_GPU();
            ELMCExecutor(&device_csc_matrix_bottom, d_x_1, d_b_1, d_graphInDegree_2);
            end_profile_GPU(&time_sptrsv2);
            
            start_profile_GPU();
            if(policy.sptrsv_method1 == SYNC_FREE || policy.sptrsv_method1 == SLFC || policy.sptrsv_method1 == ELMC)
            {                
                cudaMemcpy(x+device_csc_matrix_top.m, d_x_1, 
                           device_csc_matrix_bottom.m * policy.rhs * sizeof(val_type), 
                           cudaMemcpyDeviceToHost);                        
            }  
            else
            {   
                cudaMemcpy(x+device_csr_matrix_top.m, d_x_1, 
                           device_csc_matrix_bottom.m * policy.rhs * sizeof(val_type), 
                           cudaMemcpyDeviceToHost);                
            }
            end_profile_GPU(&time_d2h_final);   
        }

        if(policy.sptrsv_method2 == ELMR)
        {
            start_profile_GPU();
            ELMRExecutor(&device_csr_matrix_bottom, d_x_1, d_b_1, ready_2);
            end_profile_GPU(&time_sptrsv2);
            
            start_profile_GPU();
            if(policy.sptrsv_method1 == SYNC_FREE || policy.sptrsv_method1 == SLFC || policy.sptrsv_method1 == ELMC)
            {
                cudaMemcpy(x+device_csc_matrix_top.m, d_x_1, 
                           device_csr_matrix_bottom.m * policy.rhs * sizeof(val_type), 
                           cudaMemcpyDeviceToHost);    
            }  
            else
            {                    
                cudaMemcpy(x+device_csr_matrix_top.m, d_x_1, 
                           device_csr_matrix_bottom.m * policy.rhs * sizeof(val_type), 
                           cudaMemcpyDeviceToHost);
            }
            end_profile_GPU(&time_d2h_final);
        }

        }
        else
        {

        }
	}
 
	return EXIT_SUCCESS;
}

void SILU::PrintCumNnzs()
{
    for(int j = 0; j < levels.nlevel; j++)
            printf("Level %d: %d\n", j, levels.cum_Nnzs[j]);
}

int SILU::SpTRSV_Time_Breakup()
{
    printf("Factorization time: %0.6f\n", time_factorize);
    printf("Level analysis time: %0.6f\n", time_levels);
    printf("DAG analysis time: %0.6f\n", time_dag_analysis);
    printf("Split time: %0.6f\n", time_split);
    printf("Data structure setup time: %0.6f\n", time_data_structures);
    printf("Algo analysis time: %0.6f\n", time_algo_analysis);
    printf("SpTRSV 1 time: %0.6f msec\n", time_sptrsv1);
    printf("SpTRSV 2 time: %0.6f msec\n", time_sptrsv2);
    printf("SpMV Time: %0.6f msec\n", time_spmv);
    printf("Device to host memory transfer time - intermediate: %0.6f msec\n", time_d2h);
    printf("Device to host memory transfer time - end: %0.6f msec\n", time_d2h_final);
    printf("Host to device memory transfer time: %0.6f msec\n", time_h2d);    
    printf("Total pre-processing: %0.6f msec\n", time_levels + time_dag_analysis +
            time_split + time_data_structures + time_algo_analysis);
    printf("Total SpTRSV Time: %0.3f msec\n", time_sptrsv1 + time_sptrsv2 + time_spmv + time_d2h + time_h2d);
    return EXIT_SUCCESS;
}

int SILU::Print_SubMatrix_stats()
{
    long totalrows = 0;
    long totalnnz = 0;
    if(policy.sptrsv_method1 == MKL || policy.sptrsv_method1 == CUSPARSE_V2)
    {
        printf("Top: rows=%d\n", TopTri_csr.m);
        printf("Top: nnzs=%d\n", TopTri_csr.nnz);        
        totalrows += TopTri_csr.m;
        totalnnz += TopTri_csr.nnz;
    }
    else
    {
        printf("Top: rows=%d\n", TopTri_csc.m);
        printf("Top: nnzs=%d\n", TopTri_csc.nnz);
        totalrows += TopTri_csc.m;
        totalnnz += TopTri_csc.nnz;
    }
    printf("SpMV: rows=%d\n", Mvp_csr.m);
    printf("SpMV: nnzs=%d\n", Mvp_csr.nnz);
    totalrows += Mvp_csr.m;
    totalnnz += Mvp_csr.nnz;

    if(policy.sptrsv_method2 == MKL || policy.sptrsv_method2 == CUSPARSE_V2)
    {
        printf("Bottom: rows=%d\n", BotTri_csr.m);
        printf("Bottom: nnzs=%d\n", BotTri_csr.nnz);
        
        totalnnz += BotTri_csr.nnz;        
    }
    else
    {
        printf("Bottom: rows=%d\n", BotTri_csc.m);
        printf("Bottom: nnzs=%d\n", BotTri_csc.nnz);
        //totalrows += BotTri_csr.m;
        totalnnz += BotTri_csr.nnz;
    }

    printf("Total rows: %d\n", totalrows);
    printf("Total nnzs: %d\n", totalnnz);
}

int SILU::setup_device_datastructures(val_type *b)
{
  printf("Submatrix Rows\n");	
    if(policy.sptrsv_method1 == SYNC_FREE || policy.sptrsv_method1 == SLFC || policy.sptrsv_method1 == ELMC)
    {
      device_csc_matrix_top.m = TopTri_csc.m;
      device_csc_matrix_top.n = TopTri_csc.n;
      device_csc_matrix_top.nnz = TopTri_csc.nnz;
      m += device_csc_matrix_top.m;
      n += device_csc_matrix_top.n;
      printf("%d:", device_csc_matrix_top.m); 
    }
    else
    {
      //printf("%d:::%d:::%d\n", TopTri_csr.m, TopTri_csr.n, TopTri_csr.nnz);
      device_csr_matrix_top.m = TopTri_csr.m;
      device_csr_matrix_top.n = TopTri_csr.n;
      device_csr_matrix_top.nnz = TopTri_csr.nnz;
      m += device_csr_matrix_top.m;
      n += device_csr_matrix_top.n;
      printf("%d:", device_csr_matrix_top.m);    
    }

    if(policy.sptrsv_method2 == SYNC_FREE  || policy.sptrsv_method2 == SLFC || policy.sptrsv_method2 == ELMC)
    {
      device_csc_matrix_bottom.m = BotTri_csc.m;
      device_csc_matrix_bottom.n = BotTri_csc.n;
      device_csc_matrix_bottom.nnz = BotTri_csc.nnz;
      m += device_csc_matrix_bottom.m;
      n += device_csc_matrix_bottom.n;
      printf(":%d\n", device_csc_matrix_bottom.m); 
    }
    else
    {

      device_csr_matrix_bottom.m = BotTri_csr.m;
      device_csr_matrix_bottom.n = BotTri_csr.n;
      device_csr_matrix_bottom.nnz = BotTri_csr.nnz;
      m += device_csr_matrix_bottom.m;
      n += device_csr_matrix_bottom.n;
      printf(":%d\n", device_csr_matrix_bottom.m);
    }


    if(policy.splitting_policy == DAG_SPLIT)
    {
      //val_type *b_ = (val_type *)malloc(sizeof(val_type) * (m));
      val_type *b_;
      cudaError_t status = cudaMallocHost((void **)&b_, sizeof(val_type)*m);
      printf("Pinned: %d\n", m);
      if(status !=cudaSuccess)
	printf("Error allocating pinned host memory\n");
      //for(int i = 0; i < m; i++)
      //	{
      //  b_[i]=b[i];
      //	}
      memcpy(b_, b, sizeof(val_type)*(m));   
        for(int i = 0; i < m; i++)
        {
           b[i] = b_[levels.levelItem[i]];                
        }
	cudaFree(b_);
    }

    if(policy.sptrsv_method1 == SYNC_FREE || policy.sptrsv_method1 == SLFC || policy.sptrsv_method1 == ELMC)
    {
      CUDA_RT_CALL(cudaMalloc((void **)&device_csc_matrix_top.d_cscColPtr, (TopTri_csc.n+1) * sizeof(ind_type)));
      CUDA_RT_CALL(cudaMalloc((void **)&device_csc_matrix_top.d_cscRowIdx, TopTri_csc.nnz  * sizeof(ind_type)));
      CUDA_RT_CALL(cudaMalloc((void **)&device_csc_matrix_top.d_cscVal,    TopTri_csc.nnz  * sizeof(val_type)));
      CUDA_RT_CALL(cudaMemcpy(device_csc_matrix_top.d_cscColPtr, TopTri_csc.cscColPtr, (TopTri_csc.n+1) * sizeof(ind_type),   cudaMemcpyHostToDevice));
      CUDA_RT_CALL(cudaMemcpy(device_csc_matrix_top.d_cscRowIdx, TopTri_csc.cscRowIdx, (TopTri_csc.nnz)  * sizeof(ind_type),   cudaMemcpyHostToDevice));
      CUDA_RT_CALL(cudaMemcpy(device_csc_matrix_top.d_cscVal, TopTri_csc.cscVal,    (TopTri_csc.nnz)  * sizeof(val_type),   cudaMemcpyHostToDevice));
       
      if(policy.sptrsv_method1 == SYNC_FREE)              
      {
	CUDA_RT_CALL(cudaMalloc((void **)&device_csc_matrix_top.d_left_sum, sizeof(val_type) * TopTri_csc.m * policy.rhs));
	CUDA_RT_CALL(cudaMalloc((void **)&d_b, TopTri_csc.m * policy.rhs * sizeof(val_type)));
	CUDA_RT_CALL(cudaMemcpy(d_b, b, TopTri_csc.m * policy.rhs * sizeof(val_type), cudaMemcpyHostToDevice));
      }  
      else
      {        
        CUDA_RT_CALL(cudaMalloc((void **)&device_csc_matrix_top.diag, sizeof(val_type) * TopTri_csc.m * policy.rhs));
        val_type *diag = (val_type *)malloc(sizeof(val_type) * TopTri_csc.m);
        for(int i = 0; i < TopTri_csc.m; i++)
            diag[i] = 1.0;
        CUDA_RT_CALL(cudaMemcpy(device_csc_matrix_top.diag, diag, TopTri_csc.m * sizeof(val_type), cudaMemcpyHostToDevice));
        free(diag);
        CUDA_RT_CALL(cudaMalloc((void **)&device_csc_matrix_top.jlev, sizeof(val_type) * TopTri_csc.m * policy.rhs));
        CUDA_RT_CALL(cudaMemcpy(device_csc_matrix_top.jlev, levelsTop.levelItem, (TopTri_csc.m)  * sizeof(ind_type),   cudaMemcpyHostToDevice));
      }
    }    

    if(policy.sptrsv_method1 == CUSPARSE_V2 || policy.sptrsv_method1 == ELMR)
	{
	  //PrintCsrMatrix(&TopTri_csr);
      //printf("%d:::%d:::%d\n", TopTri_csr.m, TopTri_csr.n, TopTri_csr.nnz);
      CUDA_RT_CALL(cudaMalloc((void **)&device_csr_matrix_top.d_csrRowPtr, (TopTri_csr.m+1) * sizeof(ind_type)));
	  CUDA_RT_CALL(cudaMalloc((void **)&device_csr_matrix_top.d_csrColIdx, TopTri_csr.nnz * sizeof(ind_type)));
	  CUDA_RT_CALL(cudaMalloc((void **)&device_csr_matrix_top.d_csrVal,    TopTri_csr.nnz  * sizeof(val_type)));
	  CUDA_RT_CALL(cudaMemcpy(device_csr_matrix_top.d_csrRowPtr, TopTri_csr.csrRowPtr, (TopTri_csr.m+1) * sizeof(ind_type),   cudaMemcpyHostToDevice));
	  CUDA_RT_CALL(cudaMemcpy(device_csr_matrix_top.d_csrColIdx, TopTri_csr.csrColIdx, (TopTri_csr.nnz)  * sizeof(ind_type),   cudaMemcpyHostToDevice));
	  CUDA_RT_CALL(cudaMemcpy(device_csr_matrix_top.d_csrVal, TopTri_csr.csrVal, (TopTri_csr.nnz)  * sizeof(val_type),   cudaMemcpyHostToDevice));
	  CUDA_RT_CALL(cudaMalloc((void **)&d_b, TopTri_csr.m * policy.rhs * sizeof(val_type)));
	  CUDA_RT_CALL(cudaMemcpy(d_b, b, TopTri_csr.m * policy.rhs * sizeof(val_type), cudaMemcpyHostToDevice));
      
      if(policy.sptrsv_method1 == ELMR)
      {
        cudaMalloc((void **)&device_csr_matrix_top.diag, sizeof(val_type) * TopTri_csr.m * policy.rhs);
        val_type *diag = (val_type *)malloc(sizeof(val_type) * TopTri_csr.m);
        for(int i = 0; i < TopTri_csr.m; i++)
            diag[i] = 1.0;
        cudaMemcpy(device_csr_matrix_top.diag, diag, TopTri_csr.m * sizeof(val_type), cudaMemcpyHostToDevice);
               
        cudaMalloc((void **)&device_csr_matrix_top.jlev, sizeof(val_type) * TopTri_csr.m * policy.rhs);
        CUDA_RT_CALL(cudaMemcpy(device_csr_matrix_top.jlev, levelsTop.levelItem, (TopTri_csr.m)  * sizeof(ind_type),   cudaMemcpyHostToDevice));

        cudaMalloc((void **)&ready_1, sizeof(char) * TopTri_csr.m * policy.rhs);
        for(int i = 0; i < TopTri_csr.m; i++)
             diag[i] = 0.0;
        free(diag);
        cudaMemcpy(ready_1, diag, TopTri_csr.m * sizeof(val_type), cudaMemcpyHostToDevice);
        
      }
	          
	}

    if(policy.spmv_method = SPMV_CSR)
    {
        if(policy.sptrsv_method1 == SYNC_FREE || policy.sptrsv_method1 == SLFC || policy.sptrsv_method1 == ELMC)
        {
	  CUDA_RT_CALL(cudaMalloc((void **)&d_x, TopTri_csc.m * policy.rhs * sizeof(val_type)));
            if(policy.sptrsv_method1 == SYNC_FREE)
            {
	      CUDA_RT_CALL(cudaMemset(d_x, 0, TopTri_csc.m * policy.rhs * sizeof(val_type)));    
            }
            else
            {
	      CUDA_RT_CALL(cudaMemcpy(d_x, b, TopTri_csc.m * policy.rhs * sizeof(val_type), cudaMemcpyHostToDevice));
            }
            
            if(policy.sptrsv_method2 == SYNC_FREE || policy.sptrsv_method2 == SLFC || policy.sptrsv_method2 == ELMC)
            {    
	      CUDA_RT_CALL(cudaMalloc((void **)&d_b_1, BotTri_csc.m * policy.rhs * sizeof(val_type)));
	      CUDA_RT_CALL(cudaMemcpy(d_b_1, b+TopTri_csc.m, BotTri_csc.m * policy.rhs * sizeof(val_type), cudaMemcpyHostToDevice));
            }
            else
            {
	      CUDA_RT_CALL(cudaMalloc((void **)&d_b_1, BotTri_csr.m * policy.rhs * sizeof(val_type)));
	      CUDA_RT_CALL(cudaMemcpy(d_b_1, b+TopTri_csc.m, BotTri_csr.m * policy.rhs * sizeof(val_type), cudaMemcpyHostToDevice));
            }            
        }
        else
        {
            CUDA_RT_CALL(cudaMalloc((void **)&d_x, TopTri_csr.m * policy.rhs * sizeof(val_type)));
            CUDA_RT_CALL(cudaMemset(d_x, 0, TopTri_csr.m * policy.rhs * sizeof(val_type)));
            if(policy.sptrsv_method2 == SYNC_FREE || policy.sptrsv_method2 == SLFC || policy.sptrsv_method2 == ELMC)
            {    
	      CUDA_RT_CALL(cudaMalloc((void **)&d_b_1, BotTri_csc.m * policy.rhs * sizeof(val_type)));
	      CUDA_RT_CALL(cudaMemcpy(d_b_1, b+TopTri_csr.m, BotTri_csc.m * policy.rhs * sizeof(val_type), cudaMemcpyHostToDevice));                
            }
            else
            {
                if(BotTri_csr.m > 0)
                {
                    //printf("But not here\n");   
                    CUDA_RT_CALL(cudaMalloc((void **)&d_b_1, BotTri_csr.m * policy.rhs * sizeof(val_type)));
                    CUDA_RT_CALL(cudaMemcpy(d_b_1, b+TopTri_csr.m, BotTri_csr.m * policy.rhs * sizeof(val_type), cudaMemcpyHostToDevice));
                }                
            }
        }

    }

    if(policy.sptrsv_method1 == MKL)
    {
        sp_status = mkl_sparse_d_create_csr(&csrLU1, SPARSE_INDEX_BASE_ZERO, TopTri_csr.m,
                                            TopTri_csr.n, TopTri_csr.csrRowPtr, TopTri_csr.csrRowPtr + 1,
                                            TopTri_csr.csrColIdx, TopTri_csr.csrVal);
        if(sp_status == SPARSE_STATUS_SUCCESS)
        {
            //printf("MKL LU create succeed\n");
        }
    }    
    
    if(policy.spmv_method == SPMV_CSR && policy.spmv_platform == SPMV_CPU)
    {
        if(Mvp_csr.m > 0)
        sp_status = mkl_sparse_d_create_csr(&csrA, SPARSE_INDEX_BASE_ZERO, Mvp_csr.m, Mvp_csr.n,
                                            Mvp_csr.csrRowPtr, Mvp_csr.csrRowPtr + 1, 
                                            Mvp_csr.csrColIdx, Mvp_csr.csrVal);        
    }
    
    if(policy.spmv_method == SPMV_CSR && policy.spmv_platform == SPMV_GPU)
      {
        cudaError_t cudaStatus;
        device_csr_matrix_mvp.m = Mvp_csr.m;
        device_csr_matrix_mvp.n = Mvp_csr.n;
        device_csr_matrix_mvp.nnz = Mvp_csr.nnz;
	    //printf("%d::%d::%d\n", Mvp_csr.m, Mvp_csr.n, Mvp_csr.nnz);
        if(Mvp_csr.m > 0)
        {
            CUDA_RT_CALL(cudaMalloc((void **)&device_csr_matrix_mvp.d_csrRowPtr, (Mvp_csr.m+1) * sizeof(ind_type)));
            CUDA_RT_CALL(cudaMemcpy(device_csr_matrix_mvp.d_csrRowPtr, Mvp_csr.csrRowPtr, (Mvp_csr.m+1) * sizeof(ind_type),   cudaMemcpyHostToDevice));    
        }
        if(Mvp_csr.nnz > 0)
        {
            CUDA_RT_CALL(cudaMalloc((void **)&device_csr_matrix_mvp.d_csrColIdx, Mvp_csr.nnz  * sizeof(ind_type)));
            CUDA_RT_CALL(cudaMalloc((void **)&device_csr_matrix_mvp.d_csrVal, Mvp_csr.nnz  * sizeof(val_type)));
            CUDA_RT_CALL(cudaMemcpy(device_csr_matrix_mvp.d_csrColIdx, Mvp_csr.csrColIdx, (Mvp_csr.nnz)  * sizeof(ind_type),   cudaMemcpyHostToDevice));
            CUDA_RT_CALL(cudaMemcpy(device_csr_matrix_mvp.d_csrVal, Mvp_csr.csrVal,    (Mvp_csr.nnz)  * sizeof(val_type),   cudaMemcpyHostToDevice));     
        }
                
      }
    
    if(policy.sptrsv_method2 == SYNC_FREE || policy.sptrsv_method2 == SLFC || policy.sptrsv_method2 == ELMC)
      {
        if(BotTri_csc.m > 0)
        {CUDA_RT_CALL(cudaMalloc((void **)&device_csc_matrix_bottom.d_cscColPtr, (BotTri_csc.m+1) * sizeof(ind_type)));
        CUDA_RT_CALL(cudaMalloc((void **)&device_csc_matrix_bottom.d_cscRowIdx, BotTri_csc.nnz  * sizeof(ind_type)));
        CUDA_RT_CALL(cudaMalloc((void **)&device_csc_matrix_bottom.d_cscVal, BotTri_csc.nnz  * sizeof(val_type)));
        CUDA_RT_CALL(cudaMemcpy(device_csc_matrix_bottom.d_cscColPtr, BotTri_csc.cscColPtr, (BotTri_csc.n+1) * sizeof(ind_type),   cudaMemcpyHostToDevice));
       CUDA_RT_CALL(cudaMemcpy(device_csc_matrix_bottom.d_cscRowIdx, BotTri_csc.cscRowIdx, (BotTri_csc.nnz)  * sizeof(ind_type),   cudaMemcpyHostToDevice));
       CUDA_RT_CALL(cudaMemcpy(device_csc_matrix_bottom.d_cscVal, BotTri_csc.cscVal,    (BotTri_csc.nnz)  * sizeof(val_type),   cudaMemcpyHostToDevice));
        
        CUDA_RT_CALL(cudaMalloc((void **)&d_x_1, BotTri_csc.n * policy.rhs * sizeof(val_type)));
        CUDA_RT_CALL(cudaMalloc((void **)&device_csc_matrix_bottom.d_left_sum, sizeof(val_type) * BotTri_csc.m * policy.rhs));
        
        if(policy.sptrsv_method2 == SYNC_FREE || policy.sptrsv_method2 == SLFC || policy.sptrsv_method2 == ELMC)              
        {
	  CUDA_RT_CALL(cudaMemset(d_x_1, 0, BotTri_csc.n * policy.rhs * sizeof(val_type)));
        }    
        else
        {
	  CUDA_RT_CALL(cudaMemset(d_x_1, 0, BotTri_csr.n * policy.rhs * sizeof(val_type)));                            
        }
        if(policy.sptrsv_method2 == SLFC || policy.sptrsv_method2 == ELMC)
        {
	  CUDA_RT_CALL(cudaMalloc((void **)&device_csc_matrix_bottom.diag, sizeof(val_type) * BotTri_csc.m * policy.rhs));
            val_type *diag = (val_type *)malloc(sizeof(val_type) * BotTri_csc.m);
            for(int i = 0; i < BotTri_csc.m; i++)
                diag[i] = 1.0;
            CUDA_RT_CALL(cudaMemcpy(device_csc_matrix_bottom.diag, diag, BotTri_csc.m * sizeof(val_type), cudaMemcpyHostToDevice));
            
            CUDA_RT_CALL(cudaMalloc((void **)&device_csc_matrix_bottom.jlev, sizeof(val_type) * BotTri_csc.m * policy.rhs));
            CUDA_RT_CALL(cudaMemcpy(device_csc_matrix_bottom.jlev, levelsBot.levelItem, (BotTri_csc.m)  * sizeof(ind_type),   
                                 cudaMemcpyHostToDevice));
            if(policy.sptrsv_method1 == SYNC_FREE || policy.sptrsv_method1 == SLFC || policy.sptrsv_method1 == ELMC)
	      {
		CUDA_RT_CALL(cudaMemcpy(d_x_1, b+TopTri_csc.m, BotTri_csc.m * policy.rhs * sizeof(val_type), cudaMemcpyHostToDevice));
	      }
            else
	      {
		CUDA_RT_CALL(cudaMemcpy(d_x_1, b+TopTri_csr.m, BotTri_csc.m * policy.rhs * sizeof(val_type), cudaMemcpyHostToDevice));
	      }
        }
        }        
      }    
    
    if(policy.sptrsv_method2 == CUSPARSE_V2 || policy.sptrsv_method2 == ELMR)
      {
	CUDA_RT_CALL(cudaMalloc((void **)&device_csr_matrix_bottom.d_csrVal, BotTri_csr.nnz  * sizeof(val_type)));
	CUDA_RT_CALL(cudaMalloc((void **)&device_csr_matrix_bottom.d_csrColIdx, BotTri_csr.nnz  * sizeof(ind_type)));
        CUDA_RT_CALL(cudaMalloc((void **)&device_csr_matrix_bottom.d_csrRowPtr, (BotTri_csr.m+1) * sizeof(ind_type)));
        CUDA_RT_CALL(cudaMemcpy(device_csr_matrix_bottom.d_csrColIdx, BotTri_csr.csrColIdx, (BotTri_csr.nnz) * sizeof(ind_type),   cudaMemcpyHostToDevice));
        CUDA_RT_CALL(cudaMemcpy(device_csr_matrix_bottom.d_csrRowPtr, BotTri_csr.csrRowPtr, (BotTri_csr.n+1)  * sizeof(ind_type),   cudaMemcpyHostToDevice));
        CUDA_RT_CALL(cudaMemcpy(device_csr_matrix_bottom.d_csrVal, BotTri_csr.csrVal,    (BotTri_csr.nnz)  * sizeof(val_type),   cudaMemcpyHostToDevice));
        CUDA_RT_CALL(cudaMalloc((void **)&d_b_1, BotTri_csr.m * policy.rhs * sizeof(val_type)));
        if(policy.sptrsv_method1 == SYNC_FREE || policy.sptrsv_method1 == SLFC || policy.sptrsv_method1 == ELMC)
	  {
	    CUDA_RT_CALL(cudaMemcpy(d_b_1, b+TopTri_csc.m, BotTri_csr.m * policy.rhs * sizeof(val_type), cudaMemcpyHostToDevice));     
        }
        else
	  {
           CUDA_RT_CALL(cudaMemcpy(d_b_1, b+TopTri_csr.m, BotTri_csr.m * policy.rhs * sizeof(val_type), cudaMemcpyHostToDevice));                 
	  }
        CUDA_RT_CALL(cudaMalloc((void **)&d_x_1, BotTri_csr.n * policy.rhs * sizeof(val_type)));
        CUDA_RT_CALL(cudaMemset(d_x_1, 0, BotTri_csr.n * policy.rhs * sizeof(val_type)));
	
        if(policy.sptrsv_method2 == ELMR)
	  {
	    CUDA_RT_CALL(cudaMalloc((void **)&device_csr_matrix_bottom.diag, sizeof(val_type) * BotTri_csr.m * policy.rhs));
	    val_type *diag = (val_type *)malloc(sizeof(val_type) * BotTri_csr.m);
           for(int i = 0; i < BotTri_csr.m; i++)
             diag[i] = 1.0;
           CUDA_RT_CALL(cudaMemcpy(device_csr_matrix_bottom.diag, diag, BotTri_csr.m * sizeof(val_type), cudaMemcpyHostToDevice));
           CUDA_RT_CALL(cudaMalloc((void **)&device_csr_matrix_bottom.jlev, sizeof(val_type) * BotTri_csr.m * policy.rhs));
           CUDA_RT_CALL(cudaMemcpy(device_csr_matrix_bottom.jlev, levelsBot.levelItem, (BotTri_csr.m)  * sizeof(ind_type),   cudaMemcpyHostToDevice));
           CUDA_RT_CALL(cudaMalloc((void **)&ready_2, sizeof(char) * BotTri_csr.m * policy.rhs));
           for(int i = 0; i < BotTri_csr.m; i++)
             diag[i] = 0.0;
           CUDA_RT_CALL(cudaMemcpy(ready_2, diag, BotTri_csr.m * sizeof(val_type), cudaMemcpyHostToDevice));

           if(policy.sptrsv_method1 == SYNC_FREE || policy.sptrsv_method1 == SLFC || policy.sptrsv_method1 == ELMC)
           {
                CUDA_RT_CALL(cudaMemcpy(d_x_1, b+TopTri_csc.m, BotTri_csr.m * policy.rhs * sizeof(val_type), cudaMemcpyHostToDevice)); 
           }
            else
	      {
                CUDA_RT_CALL(cudaMemcpy(d_x_1, b+TopTri_csr.m, BotTri_csr.m * policy.rhs * sizeof(val_type), cudaMemcpyHostToDevice));
	      }
           
        }
      }    
    
    if(policy.sptrsv_method2 == MKL)
      {        
        if(BotTri_csr.m > 0)
        sp_status = mkl_sparse_d_create_csr(&csrLU, SPARSE_INDEX_BASE_ZERO, BotTri_csr.m,
                                            BotTri_csr.n, BotTri_csr.csrRowPtr, BotTri_csr.csrRowPtr + 1,
                                            BotTri_csr.csrColIdx, BotTri_csr.csrVal);
        if(sp_status == SPARSE_STATUS_SUCCESS)
        {
            //printf("MKL LU create succeed\n");
        }        
    }      

    cudaError_t err = cudaGetLastError();
    if(err != cudaSuccess)
    {
        printf("CUDA Error: (%s)\n", cudaGetErrorString(err));
        return EXIT_FAILURE;
    }

	return EXIT_SUCCESS;
}

void SILU::PrintCsrMatrix(csr_matrix *mat, char rows, int norows)
{
	if(rows == 'a')
    {
        for(int i = 0; i < mat->m; i++)
        {            
            for(int j = mat->csrRowPtr[i]; j < mat->csrRowPtr[i+1]; j++)
            {
                printf("[%d, %ld] = %f ", i, mat->csrColIdx[j], mat->csrVal[j]);
            }
            printf("\n");      
        }    
    }
    else
    {
        for(int i = 0; i < norows; i++)
        {
            printf("entries: %d:", mat->csrRowPtr[i+1] -mat->csrRowPtr[i]);
            for(int j = mat->csrRowPtr[i]; j < mat->csrRowPtr[i+1]; j++)
            {
                printf("[%d, %ld] = %f ", i, mat->csrColIdx[j], mat->csrVal[j]);
            }
            printf("\n.");      
        }   
    }
    
}

void SILU::PrintCscMatrix(csc_matrix *mat)
{
	for(int i = 0; i < mat->m; i++)
	{
		for(int j = mat->cscColPtr[i]; j < mat->cscColPtr[i+1]; j++)
		{
			printf("[%ld, %d] = %f ", mat->cscRowIdx[j], i, mat->cscVal[j]);
		}
		printf("\n");
	}	
}

void SILU::PrintLevels()
{
	for(int i = 0; i < levels.nlevel; i++)
	{
		for(int j = levels.levelPtr[i]; j <levels.levelPtr[i+1]; j++)
		{
			printf("%ld, ", levels.levelItem[j]);
		}
		printf("\n");
	}
}

void SILU::PrintLevels(level_info *levels)
{
    for(int i = 0; i < levels->nlevel; i++)
    {
        for(int j = levels->levelPtr[i]; j <levels->levelPtr[i+1]; j++)
        {
            printf("%ld, ", levels->levelItem[j]);
        }
        printf("\n");
    }
}

int SILU::split_dag(csr_matrix *tri_mat)
{
    if(policy.spmv_method == SPMV_CSR)
	{	
        TopTri_csr.csrRowPtr=(ind_type *)malloc((levels.cum_Rows[policy.split_at_level]+1) * sizeof(ind_type));
		TopTri_csr.csrColIdx=(ind_type *)malloc((levels.cum_Nnzs[policy.split_at_level]) * sizeof(ind_type));
		TopTri_csr.csrVal=(val_type *)malloc((levels.cum_Nnzs[policy.split_at_level]) * sizeof(val_type));
    	TopTri_csr.m = levels.cum_Rows[policy.split_at_level]; 	
		TopTri_csr.n = TopTri_csr.m;
		TopTri_csr.nnz = levels.cum_Nnzs[policy.split_at_level];
		
		int lvi = 1;
		int nnzptr = 0;
		int rowptr = 0;
        int rowcntr = 0;
#ifdef SAVE_MATRIX
	//printf("I WAS HERE\n");
	std::ofstream ia("ia.txt", std::ios::out);
	if(ia.is_open())
	  {
	  }
	else
	  {
	    printf("Failed creating ia file\n");
	    return EXIT_FAILURE;
	  }
	std::ofstream ib("ib.txt", std::ios::out);
	if(ib.is_open())
	  {
	  }
	else
	  {
	    printf("Failed creating ib file\n");
	    return EXIT_FAILURE;
	  }
	std::ofstream val("val.txt", std::ios::out);
	if(val.is_open())
	  {
	  }
	else
	  {
	    printf("Failed creating val file\n");
	    return EXIT_FAILURE;
	  }
#endif
		TopTri_csr.csrRowPtr[rowptr++] = 0;
#ifdef SAVE_MATRIX
		ia << TopTri_csr.csrRowPtr[rowptr-1] << ",";
#endif
    	//printf("Last row number for first split: %d\n", levels.cum_Rows[policy.split_at_level]);
    	// Arrange Top triangle
		while(levels.levelPtr[lvi-1] != levels.cum_Rows[policy.split_at_level])
		{
			for(int i = levels.levelPtr[lvi-1]; i < levels.levelPtr[lvi]; i++)
			{				
			    int node = levels.levelItem[i];
				for(int j = tri_mat->csrRowPtr[node]; j < tri_mat->csrRowPtr[node+1]; j++)
				{
					TopTri_csr.csrVal[nnzptr] = tri_mat->csrVal[j];
					
					TopTri_csr.csrColIdx[nnzptr++] = levels.levelItemNewRowIdx[tri_mat->csrColIdx[j]];
#ifdef SAVE_MATRIX
					ib << TopTri_csr.csrColIdx[nnzptr-1] << ",";
					if(rowcntr < tri_mat->m/2)
					  {
					    val << 300 << ",";					   
					  }
					else if(TopTri_csr.csrColIdx[nnzptr-1] < tri_mat->m/2)
					  {
					    val << 600 << ",";					    
					  }
					else
					  {
					    val << 900 << ",";
					  }
#endif
				}
				rowcntr++;
				
				TopTri_csr.csrRowPtr[rowptr++] = nnzptr;
#ifdef SAVE_MATRIX
				ia << TopTri_csr.csrRowPtr[rowptr-1] << ",";
#endif
			}
			lvi++;
		}
#ifdef SAVE_MATRIX
		ia.close();
		ib.close();
		val.close();
#endif
        ind_type rem_rows = tri_mat->m - levels.cum_Rows[policy.split_at_level];
        ind_type rem_nnz = tri_mat->nnz - nnzptr;
        //printf("Remaining ROWS: %d\r", rem_rows);
	//printf("Ok till here\n");
	
        Mvp_csr.csrRowPtr=(ind_type *)malloc((rem_rows+1) * sizeof(ind_type));
        Mvp_csr.csrColIdx=(ind_type *)malloc((rem_nnz) * sizeof(ind_type));
        Mvp_csr.csrVal=(val_type *)malloc((rem_nnz) * sizeof(val_type));
        BotTri_csr.csrRowPtr=(ind_type *)malloc((rem_rows+1) * sizeof(ind_type));
        BotTri_csr.csrColIdx=(ind_type *)malloc((rem_nnz) * sizeof(ind_type));
        BotTri_csr.csrVal=(val_type *)malloc((rem_nnz) * sizeof(val_type));
        int nnzptrMvp = 0;
        int nnzptrBotTri = 0;
        int rowptrMvp = 0;
        int rowptrBotTri = 0;
        Mvp_csr.csrRowPtr[rowptrMvp++]=0;
        BotTri_csr.csrRowPtr[rowptrBotTri++]=0;
		// Arrange SpMV block and Bottom triangle
		while(levels.levelPtr[lvi-1] != levels.cum_Rows[levels.nlevel-1])
		{
			for(int i = levels.levelPtr[lvi-1]; i < levels.levelPtr[lvi]; i++)
			{				
			    int node = levels.levelItem[i];
				for(int j = tri_mat->csrRowPtr[node]; j < tri_mat->csrRowPtr[node+1]; j++)
				{
					if(levels.levelItemNewRowIdx[tri_mat->csrColIdx[j]] < levels.cum_Rows[policy.split_at_level])
					{
						Mvp_csr.csrVal[nnzptrMvp] = tri_mat->csrVal[j];
						Mvp_csr.csrColIdx[nnzptrMvp++] = levels.levelItemNewRowIdx[tri_mat->csrColIdx[j]];
					}
					else
					{
						BotTri_csr.csrVal[nnzptrBotTri] = tri_mat->csrVal[j];
						BotTri_csr.csrColIdx[nnzptrBotTri++] = levels.levelItemNewRowIdx[tri_mat->csrColIdx[j]] - 
						                                           (levels.cum_Rows[policy.split_at_level]);						
					}	
				}
			    rowcntr++;	
				Mvp_csr.csrRowPtr[rowptrMvp++] = nnzptrMvp;
				BotTri_csr.csrRowPtr[rowptrBotTri++] = nnzptrBotTri;				
			}
						
			lvi++;
		}
		
		Mvp_csr.m = tri_mat->m - levels.cum_Rows[policy.split_at_level];
		Mvp_csr.n = levels.cum_Rows[policy.split_at_level];
		Mvp_csr.nnz=nnzptrMvp;
		BotTri_csr.m = tri_mat->m - levels.cum_Rows[policy.split_at_level];
		BotTri_csr.n = BotTri_csr.m;
		BotTri_csr.nnz=nnzptrBotTri;
		
		// Sort column indices within a row
		/*	
		for (int i = 0; i < TopTri_csr.m; i++)
        {
            quick_sort_key_val_pair<ind_type, val_type>(&TopTri_csr.csrColIdx[TopTri_csr.csrRowPtr[i]],
                                              &TopTri_csr.csrVal[TopTri_csr.csrRowPtr[i]],
                                              TopTri_csr.csrRowPtr[i+1]-TopTri_csr.csrRowPtr[i]);
        }
	printf("Ok till here 2\n");
        for (int i = 0; i < BotTri_csr.m; i++)
        {
            quick_sort_key_val_pair<ind_type, val_type>(&BotTri_csr.csrColIdx[BotTri_csr.csrRowPtr[i]],
                                              &BotTri_csr.csrVal[BotTri_csr.csrRowPtr[i]],
                                              BotTri_csr.csrRowPtr[i+1]-BotTri_csr.csrRowPtr[i]);
        }
	printf("Ok till here 3\n");
        for (int i = 0; i < Mvp_csr.m; i++)
        {
            quick_sort_key_val_pair<ind_type, val_type>(&Mvp_csr.csrColIdx[Mvp_csr.csrRowPtr[i]],
                                              &Mvp_csr.csrVal[Mvp_csr.csrRowPtr[i]],
                                              Mvp_csr.csrRowPtr[i+1]-Mvp_csr.csrRowPtr[i]);
        }
		*/
        //PrintCsrMatrix(&TopTri_csr,'a',10);

        if(policy.sptrsv_method1 == SYNC_FREE)
        {        
            matrix_transpose(&TopTri_csr, &TopTri_csc);
            free(TopTri_csr.csrRowPtr);
            free(TopTri_csr.csrColIdx);
            free(TopTri_csr.csrVal);
            //printf("Cam here\n");    
        }
        if(policy.sptrsv_method2 == SYNC_FREE)
        {
            matrix_transpose(&BotTri_csr, &BotTri_csc);    
            free(BotTri_csr.csrRowPtr);
            free(BotTri_csr.csrColIdx);
            free(BotTri_csr.csrVal);
        }
        if(policy.sptrsv_method1 == SLFC || policy.sptrsv_method1 == ELMC)
        {            
            matrix_transpose(&TopTri_csr, &TopTri_csc);
            findlevels_csr(&TopTri_csr, &levelsTop);
            //printf("I was here\n");
            free(TopTri_csr.csrRowPtr);
            free(TopTri_csr.csrColIdx);
            free(TopTri_csr.csrVal);
        }

        if(policy.sptrsv_method1 == ELMR)
        {
            matrix_transpose(&TopTri_csr, &TopTri_csc);
            findlevels_csr(&TopTri_csr, &levelsTop);   
        }

        if(policy.sptrsv_method2 == SLFC || policy.sptrsv_method2 == ELMC)
        {         
            matrix_transpose(&BotTri_csr, &BotTri_csc);
            findlevels_csr(&BotTri_csr, &levelsBot);    
            free(BotTri_csr.csrRowPtr);
            free(BotTri_csr.csrColIdx);
            free(BotTri_csr.csrVal);
        }
        if(policy.sptrsv_method2 == ELMR)
        {
            matrix_transpose(&BotTri_csr, &BotTri_csc);
            findlevels_csr(&BotTri_csr, &levelsBot);    
        }
        //PrintCsrMatrix(&TopTri_csr,'a');
	}
	//printf("Ok till here 4\n");
	return EXIT_SUCCESS;
}

int SILU::split_dag2(csr_matrix *tri_mat)
{
    if(policy.spmv_method == SPMV_CSR)
    {   

        TopTri_csr.csrRowPtr=(ind_type *)malloc((tri_mat->m+1) * sizeof(ind_type));
        TopTri_csr.csrColIdx=(ind_type *)malloc((tri_mat->nnz) * sizeof(ind_type));
        TopTri_csr.csrVal=(val_type *)malloc((tri_mat->nnz) * sizeof(val_type));
        TopTri_csr.m = tri_mat->m;  
        TopTri_csr.n = TopTri_csr.m;
        TopTri_csr.nnz = tri_mat->nnz;
        
        int lvi = 1;
        int nnzptr = 0;
        int rowptr = 0;
        TopTri_csr.csrRowPtr[rowptr++] = 0;
        
        // Arrange Top triangle
        while(levels.levelPtr[lvi-1] != levels.cum_Rows[policy.split_at_level])
        {            
            for(int i = levels.levelPtr[lvi-1]; i < levels.levelPtr[lvi]; i++)
            {               
                int node = levels.levelItem[i];
                for(int j = tri_mat->csrRowPtr[node]; j < tri_mat->csrRowPtr[node+1]; j++)
                {
                    TopTri_csr.csrVal[nnzptr] = tri_mat->csrVal[j];                 
                    TopTri_csr.csrColIdx[nnzptr++] = levels.levelItemNewRowIdx[tri_mat->csrColIdx[j]];
                }        
                TopTri_csr.csrRowPtr[rowptr++] = nnzptr;
            }
            lvi++;
        }

	
        
        ind_type rem_rows = tri_mat->m - levels.cum_Rows[policy.split_at_level];
        ind_type rem_nnz = tri_mat->nnz - nnzptr;
        
        BotTri_csr.csrRowPtr=(ind_type *)malloc((rem_rows+1) * sizeof(ind_type));
        BotTri_csr.csrColIdx=(ind_type *)malloc((rem_nnz) * sizeof(ind_type));
        BotTri_csr.csrVal=(val_type *)malloc((rem_nnz) * sizeof(val_type));

        int nnzptrBotTri = 0;
        int rowptrBotTri = 0;
        BotTri_csr.csrRowPtr[rowptrBotTri++]=0;

        while(levels.levelPtr[lvi-1] != levels.cum_Rows[levels.nlevel-1])
        {
            for(int i = levels.levelPtr[lvi-1]; i < levels.levelPtr[lvi]; i++)
            {               
                int node = levels.levelItem[i];
                for(int j = tri_mat->csrRowPtr[node]; j < tri_mat->csrRowPtr[node+1]; j++)
                {
                    if(levels.levelItemNewRowIdx[tri_mat->csrColIdx[j]] < levels.cum_Rows[policy.split_at_level])
                    {
                        TopTri_csr.csrVal[nnzptr] = tri_mat->csrVal[j];                 
                        TopTri_csr.csrColIdx[nnzptr++] = levels.levelItemNewRowIdx[tri_mat->csrColIdx[j]];
                    }
                    else
                    {
                        BotTri_csr.csrVal[nnzptrBotTri] = tri_mat->csrVal[j];
                        BotTri_csr.csrColIdx[nnzptrBotTri++] = levels.levelItemNewRowIdx[tri_mat->csrColIdx[j]] - 
                                                                   (levels.cum_Rows[policy.split_at_level]);                        
                    }   
                }                
                TopTri_csr.csrRowPtr[rowptr++] = nnzptr;
                BotTri_csr.csrRowPtr[rowptrBotTri++] = nnzptrBotTri;                
            }                        
            lvi++;
        }
        
        BotTri_csr.m = tri_mat->m - levels.cum_Rows[policy.split_at_level];
        BotTri_csr.n = BotTri_csr.m;
        BotTri_csr.nnz=nnzptrBotTri;
        
        for (int i = 0; i < TopTri_csr.m; i++)
        {
            quick_sort_key_val_pair<ind_type, val_type>(&TopTri_csr.csrColIdx[TopTri_csr.csrRowPtr[i]],
                                              &TopTri_csr.csrVal[TopTri_csr.csrRowPtr[i]],
                                              TopTri_csr.csrRowPtr[i+1]-TopTri_csr.csrRowPtr[i]);
        }

        for (int i = 0; i < BotTri_csr.m; i++)
        {
            quick_sort_key_val_pair<ind_type, val_type>(&BotTri_csr.csrColIdx[BotTri_csr.csrRowPtr[i]],
                                              &BotTri_csr.csrVal[BotTri_csr.csrRowPtr[i]],
                                              BotTri_csr.csrRowPtr[i+1]-BotTri_csr.csrRowPtr[i]);
        }

        for (int i = 0; i < Mvp_csr.m; i++)
        {
            quick_sort_key_val_pair<ind_type, val_type>(&Mvp_csr.csrColIdx[Mvp_csr.csrRowPtr[i]],
                                              &Mvp_csr.csrVal[Mvp_csr.csrRowPtr[i]],
                                              Mvp_csr.csrRowPtr[i+1]-Mvp_csr.csrRowPtr[i]);
        }
 
        if(policy.sptrsv_method1 == SYNC_FREE)
        {        
            matrix_transpose(&TopTri_csr, &TopTri_csc);
            free(TopTri_csr.csrRowPtr);
            free(TopTri_csr.csrColIdx);
            free(TopTri_csr.csrVal);         
        }
        if(policy.sptrsv_method2 == SYNC_FREE)
        {
            matrix_transpose(&BotTri_csr, &BotTri_csc);    
            free(BotTri_csr.csrRowPtr);
            free(BotTri_csr.csrColIdx);
            free(BotTri_csr.csrVal);
        }
        if(policy.sptrsv_method1 == SLFC || policy.sptrsv_method1 == ELMC)
        {            
            matrix_transpose(&TopTri_csr, &TopTri_csc);
            findlevels_csr(&TopTri_csr, &levelsTop);
            //printf("I was here\n");
            free(TopTri_csr.csrRowPtr);
            free(TopTri_csr.csrColIdx);
            free(TopTri_csr.csrVal);
        }

        if(policy.sptrsv_method1 == ELMR)
        {
            matrix_transpose(&TopTri_csr, &TopTri_csc);
            findlevels_csr(&TopTri_csr, &levelsTop);   
        }

        if(policy.sptrsv_method2 == SLFC || policy.sptrsv_method2 == ELMC)
        {         
            matrix_transpose(&BotTri_csr, &BotTri_csc);
            findlevels_csr(&BotTri_csr, &levelsBot);    
            free(BotTri_csr.csrRowPtr);
            free(BotTri_csr.csrColIdx);
            free(BotTri_csr.csrVal);
        }
        if(policy.sptrsv_method2 == ELMR)
        {
            matrix_transpose(&BotTri_csr, &BotTri_csc);
            findlevels_csr(&BotTri_csr, &levelsBot);    
        }
    }

    return EXIT_SUCCESS;
}

int SILU::cut_square_matrix(csr_matrix *in, csr_matrix *out, int start_row, int end_row)
{
    int nnzptr = 0;
    int rowptr = 0;
    out->csrRowPtr[rowptr++]=0; 
    for(int i = start_row; i < end_row; i++)
    {
        for(int j = in->csrRowPtr[i]; j < in->csrRowPtr[i+1]; j++)
        {
            if((in->csrColIdx[j] >= start_row) && (in->csrColIdx[j] < start_row+end_row))
            {
                out->csrColIdx[nnzptr] = in->csrColIdx[j]-start_row;
                out->csrVal[nnzptr++] = in->csrVal[j];
            }
        }
        out->csrRowPtr[rowptr++] = nnzptr;
    }
    out->nnz = nnzptr;

    return EXIT_SUCCESS;
}

int SILU::split_matrix(csr_matrix *mat, csr_matrix *tri_mat)
{
    if(policy.spmv_method == SPMV_CSR)
    {
        if(policy.sptrsv_method1 == PARK)
        {
            Top_csr.csrRowPtr=(ind_type *)malloc((policy.split_at_row+1) * sizeof(ind_type));
            Top_csr.csrColIdx=(ind_type *)malloc(mat->csrRowPtr[policy.split_at_row] * sizeof(ind_type));
            Top_csr.csrVal=(val_type *)malloc(mat->csrRowPtr[policy.split_at_row] * sizeof(val_type));
            Top_csr.m = policy.split_at_row;  
            Top_csr.n = Top_csr.m;
            cut_square_matrix(mat, &Top_csr, 0, policy.split_at_row);
            Top_csr.csrColIdx=(ind_type *)realloc(Top_csr.csrColIdx, Top_csr.nnz * sizeof(ind_type));
            Top_csr.csrVal=(val_type *)realloc(Top_csr.csrVal, Top_csr.nnz * sizeof(ind_type));
        }        

        TopTri_csr.csrRowPtr=(ind_type *)malloc((policy.split_at_row+1) * sizeof(ind_type));
        TopTri_csr.csrColIdx=(ind_type *)malloc(tri_mat->csrRowPtr[policy.split_at_row] * sizeof(ind_type));
        TopTri_csr.csrVal=(val_type *)malloc(tri_mat->csrRowPtr[policy.split_at_row] * sizeof(val_type));
        TopTri_csr.m = policy.split_at_row;  
        TopTri_csr.n = TopTri_csr.m;
        TopTri_csr.nnz = tri_mat->csrRowPtr[policy.split_at_row];
#ifdef SAVE_MATRIX

	std::ofstream ia("ia_mat.txt", std::ios::out);
	if(ia.is_open())
	  {
	  }
	else
	  {
	    printf("Failed creating ia file\n");
	    return EXIT_FAILURE;
	  }
	std::ofstream ib("ib_mat.txt", std::ios::out);
	if(ib.is_open())
	  {
	  }
	else
	  {
	    printf("Failed creating ib file\n");
	    return EXIT_FAILURE;
	  }
	std::ofstream val("val_mat.txt", std::ios::out);
	if(val.is_open())
	  {
	  }
	else
	  {
	    printf("Failed creating val file\n");
	    return EXIT_FAILURE;
	  }
#endif
        for(int i = 0; i < policy.split_at_row; i++)
        {
            TopTri_csr.csrRowPtr[i] = tri_mat->csrRowPtr[i];
#ifdef SAVE_MATRIX
	    ia << TopTri_csr.csrRowPtr[i] << ",";
#endif
            for(int j = tri_mat->csrRowPtr[i]; j < tri_mat->csrRowPtr[i+1]; j++)
            {
                TopTri_csr.csrColIdx[j] = tri_mat->csrColIdx[j];
                TopTri_csr.csrVal[j] = tri_mat->csrVal[j];
#ifdef SAVE_MATRIX
		ib << TopTri_csr.csrColIdx[j] << ",";
		if(i < tri_mat->m/2)
		  {
		    val << 300 << ",";					   
		  }
		else if(TopTri_csr.csrColIdx[j] < tri_mat->m/2)
		  {
		    val << 600 << ",";					    
		  }
		else
		  {
		    val << 900 << ",";					    
		  }
#endif
            }
        }	
        TopTri_csr.csrRowPtr[policy.split_at_row] = tri_mat->csrRowPtr[policy.split_at_row];            
#ifdef SAVE_MATRIX
        ia << TopTri_csr.csrRowPtr[policy.split_at_row] << ",";
	ia.close();
	ib.close();
	val.close();
#endif
	ind_type rem_rows = tri_mat->m - policy.split_at_row;
        ind_type rem_nnz = tri_mat->nnz - tri_mat->csrRowPtr[policy.split_at_row];

        if(policy.sptrsv_method2 == PARK)
        {
            Bot_csr.csrRowPtr=(ind_type *)malloc((rem_rows+1) * sizeof(ind_type));
            Bot_csr.csrColIdx=(ind_type *)malloc(rem_nnz * sizeof(ind_type));
            Bot_csr.csrVal=(val_type *)malloc(rem_nnz * sizeof(val_type));
            Bot_csr.m = rem_rows;  
            Bot_csr.n = Bot_csr.m;
            cut_square_matrix(mat, &Bot_csr, policy.split_at_row, mat->m);
            Bot_csr.csrColIdx=(ind_type *)realloc(Bot_csr.csrColIdx, Bot_csr.nnz * sizeof(ind_type));
            Bot_csr.csrVal=(val_type *)realloc(Bot_csr.csrVal, Bot_csr.nnz * sizeof(ind_type));            
        }

        Mvp_csr.csrRowPtr=(ind_type *)malloc((rem_rows+1) * sizeof(ind_type));
        Mvp_csr.csrColIdx=(ind_type *)malloc((rem_nnz) * sizeof(ind_type));
        Mvp_csr.csrVal=(val_type *)malloc((rem_nnz) * sizeof(val_type));
        BotTri_csr.csrRowPtr=(ind_type *)malloc((rem_rows+1) * sizeof(ind_type));
        BotTri_csr.csrColIdx=(ind_type *)malloc((rem_nnz) * sizeof(ind_type));
        BotTri_csr.csrVal=(val_type *)malloc((rem_nnz) * sizeof(val_type));

        int nnzptrMvp = 0;
        int nnzptrBotTri = 0;
        int rowptrMvp = 0;
        int rowptrBotTri = 0;
        Mvp_csr.csrRowPtr[rowptrMvp++]=0;
        BotTri_csr.csrRowPtr[rowptrBotTri++]=0;

        for(int i = policy.split_at_row; i < tri_mat->m; i++)
        {
            for(int j = tri_mat->csrRowPtr[i]; j < tri_mat->csrRowPtr[i+1]; j++)
            {
                
                if(tri_mat->csrColIdx[j] < policy.split_at_row)
                {
                    Mvp_csr.csrVal[nnzptrMvp] = tri_mat->csrVal[j];
                    Mvp_csr.csrColIdx[nnzptrMvp++] = tri_mat->csrColIdx[j];
                }
                else
                {
                    BotTri_csr.csrVal[nnzptrBotTri] = tri_mat->csrVal[j];
                    BotTri_csr.csrColIdx[nnzptrBotTri++] = tri_mat->csrColIdx[j] - policy.split_at_row;
                }
            }
            Mvp_csr.csrRowPtr[rowptrMvp++] = nnzptrMvp;
            BotTri_csr.csrRowPtr[rowptrBotTri++] = nnzptrBotTri;
        }

        Mvp_csr.m = tri_mat->m - policy.split_at_row;
        Mvp_csr.n = Mvp_csr.m;
        Mvp_csr.nnz = nnzptrMvp;
        BotTri_csr.m = tri_mat->m - policy.split_at_row;
        BotTri_csr.n = BotTri_csr.m;
        BotTri_csr.nnz=nnzptrBotTri;

        if(policy.sptrsv_method1 == MKL)
        {
            findlevels_csr(&TopTri_csr, &levelsTop);
        }
        if(policy.sptrsv_method2 == MKL)
        {
            findlevels_csr(&BotTri_csr, &levelsBot);   
        }
        if(policy.sptrsv_method1 == SYNC_FREE)
        {        
            matrix_transpose(&TopTri_csr, &TopTri_csc);
            free(TopTri_csr.csrRowPtr);
            free(TopTri_csr.csrColIdx);
            free(TopTri_csr.csrVal);
        }
        if(policy.sptrsv_method2 == SYNC_FREE)
        {
            matrix_transpose(&BotTri_csr, &BotTri_csc);    
            free(BotTri_csr.csrRowPtr);
            free(BotTri_csr.csrColIdx);
            free(BotTri_csr.csrVal);
        }
        if(policy.sptrsv_method1 == SLFC || policy.sptrsv_method1 == ELMC)
        {            
            matrix_transpose(&TopTri_csr, &TopTri_csc);
            findlevels_csr(&TopTri_csr, &levelsTop);
            free(TopTri_csr.csrRowPtr);
            free(TopTri_csr.csrColIdx);
            free(TopTri_csr.csrVal);
        }

        if(policy.sptrsv_method1 == ELMR)
        {
            matrix_transpose(&TopTri_csr, &TopTri_csc);
            findlevels_csr(&TopTri_csr, &levelsTop);   
        }

        if(policy.sptrsv_method2 == SLFC || policy.sptrsv_method2 == ELMC)
        {         
            matrix_transpose(&BotTri_csr, &BotTri_csc);
            findlevels_csr(&BotTri_csr, &levelsBot);    
            free(BotTri_csr.csrRowPtr);
            free(BotTri_csr.csrColIdx);
            free(BotTri_csr.csrVal);
        }
        if(policy.sptrsv_method2 == ELMR)
        {
            matrix_transpose(&BotTri_csr, &BotTri_csc);
            findlevels_csr(&BotTri_csr, &levelsBot);    
        }
    }
    printf("Top levels: %d\n", levelsTop.nlevel);
    printf("Bot levels: %d\n", levelsBot.nlevel);

}

int SILU::CPUGPUSplit(csr_matrix *tri_mat)
{
  int retCode = 0;
  int cpu_friendly_lvl_cnt = 0;
  int gpu_friendly_lvl_cnt = 0;
  int gpu_friendly_row_cnt = 0;
  float gpu_friendly_row_perc = 0.0;
  float cpu_friendly_lvl_perc = 0.0;
  // First check levels if suitable for CPU-GPU split
  if(levels.nlevel > MIN_LEVELS)
    {
      // Check minimum matrix size for CPU-GPU split
      if(tri_mat->m > MIN_ROWS)
	{
	  // Calculate PFR and SFL for final CPU-GPU split check
	  for(int i = 0; i < levels.nlevel; i++)
	    {
	      if(levels.num_levelRows[i] < CPU_FRIENDLY_MAX_ROWS_THRESHOLD)
		{
		  cpu_friendly_lvl_cnt += 1;                
		}
	      else
		{
		  gpu_friendly_row_cnt += levels.num_levelRows[i];
		  gpu_friendly_lvl_cnt++;
		}
	    }
	  gpu_friendly_row_perc = (gpu_friendly_row_cnt/(float)tri_mat->m)*100.0;
	  cpu_friendly_lvl_perc = (cpu_friendly_lvl_cnt/(float)levels.nlevel)*100.0;
	  policy.gpu_friendly_row_perc = gpu_friendly_row_perc;
	  policy.cpu_friendly_lvl_perc = cpu_friendly_lvl_perc;
	  policy.cpu_friendly_lvls = cpu_friendly_lvl_cnt;
	  policy.gpu_friendly_lvls = gpu_friendly_lvl_cnt;
	  policy.gpu_friendly_rows = gpu_friendly_row_cnt;
	  printf("CPU levels: %f, GPU rows: %f\n", policy.cpu_friendly_lvl_perc, policy.gpu_friendly_row_perc);
	  
	  if((gpu_friendly_row_perc >= GPU_FRIENDLY_ROWS_MIN_PERC) && 
	     (cpu_friendly_lvl_perc >= CPU_FRIENDLY_LEVELS_MIN_PERC))
	    {
	      printf("CPU-GPU split execution\n");
	      retCode = 1; // CPU-GPU split execution
	    }
	  else
	    {	   
	      printf("Not for CPU-GPU split\n");
	      retCode = 0; // Check for CPU-CPU split
	    }	  
	}      
    }
  else
    {
      retCode = 2; // Check for GPU-GPU split
    }
  return retCode;
}



int SILU::PreProcess(csr_matrix *tri_mat)
{
    int retCode = -1;
    int cpu_friendly_lvl_cnt = 0;
    int gpu_friendly_lvl_cnt = 0;
    int gpu_friendly_row_cnt = 0;
    float gpu_friendly_row_perc = 0.0;
    float cpu_friendly_lvl_perc = 0.0;

    if(tri_mat->m > MIN_ROWS)
    {
        for(int i = 0; i < levels.nlevel; i++)
        {
            if(levels.num_levelRows[i] < CPU_FRIENDLY_MAX_ROWS_THRESHOLD)
            {
                cpu_friendly_lvl_cnt += 1;                
            }
            else
            {
                gpu_friendly_row_cnt += levels.num_levelRows[i];
                gpu_friendly_lvl_cnt++;
            }
        }

        gpu_friendly_row_perc = (gpu_friendly_row_cnt/(float)tri_mat->m)*100.0;
        cpu_friendly_lvl_perc = (cpu_friendly_lvl_cnt/(float)levels.nlevel)*100.0;
        policy.gpu_friendly_row_perc = gpu_friendly_row_perc;
        policy.cpu_friendly_lvl_perc = cpu_friendly_lvl_perc;
        policy.cpu_friendly_lvls = cpu_friendly_lvl_cnt;
        policy.gpu_friendly_lvls = gpu_friendly_lvl_cnt;
        policy.gpu_friendly_rows = gpu_friendly_row_cnt;

        printf("CPU levels: %f, GPU levels: %f\n", policy.cpu_friendly_lvl_perc, policy.gpu_friendly_row_perc);
        if((gpu_friendly_row_perc >= GPU_FRIENDLY_ROWS_MIN_PERC) && 
           (cpu_friendly_lvl_perc >= CPU_FRIENDLY_LEVELS_MIN_PERC) &&
           (levels.nlevel > MIN_LEVELS))
        {
	  retCode = 1; // CPU-GPU split execution
        }
        else
        {
	    
	  if((policy.cpu_friendly_lvl_perc > MIN_LVL_PERC_FOR_CPU_EXEC) &&
	     (policy.gpu_friendly_row_perc < MAX_GPU_FRIENDLY_PERC_FOR_CPU_EXEC))
	    {
	      retCode = 2;   // CPU-CPU HTS
	    }
	  else
	    {
	      retCode = 3; //
	    }
        }
    }
    else
    {
        retCode = 0;
    }

    return retCode;
}

int SILU::DecidePlatform(csr_matrix *tri_mat)
{
    int retCode = 1;
    if((policy.cpu_friendly_lvl_perc > MIN_LVL_PERC_FOR_CPU_EXEC) &&
       (policy.gpu_friendly_row_perc < MAX_GPU_FRIENDLY_PERC_FOR_CPU_EXEC))
    {
        retCode = 1;
    }
    else
    {
        retCode = 2;
    }
    return retCode;
}

//int SILU::GPUSplitORUnified(csr_matrix *tri_mat)

int SILU::GPUGPUSplit(csr_matrix *tri_mat)
{
  int splitGPU = 0;
  int maxRows = levels.num_levelRows[0];
  int max_rows_level_indx = 0;
  int max_cols_level_indx = 0;
  int max_col_len = 0;
  int foundSplit = 0;
  float warp_iters_cols = levels.max_levelColLength[0]/WARP_SIZE;
  int max_col_iters = 0;// warp_iters_cols;
  float warp_iters_rows = levels.max_levelRowLength[0]/WARP_SIZE;
  int max_row_iters = 0; //warp_iters_rows;
  policy.split_at_level=0;

  // Set split level at level with maximum max column length
  for(int i = 0; i < levels.nlevel; i++)
    {
        warp_iters_cols = levels.max_levelColLength[i]/WARP_SIZE;
	if(warp_iters_cols > WARP_COLWISE_ITERATIONS_THRESHOLD && warp_iters_cols > max_col_iters)
        {
	  policy.split_at_level = i;
	  max_col_iters = levels.max_levelColLength[i]/WARP_SIZE;
	  max_cols_level_indx = i;
	  max_col_len=levels.max_levelColLength[i];
	  splitGPU=4; 
        }
    }

    // Check if maximum row length level is beyond selected split point
    // Adjust split level accordingly
  if(splitGPU==4)
    {
      for(int i = policy.split_at_level; i < levels.nlevel; i++)
	{
	  warp_iters_rows = levels.max_levelRowLength[i]/WARP_SIZE;	
	  if(warp_iters_rows > WARP_ROWISE_ITERATIONS_THRESHOLD)
	    {
	      policy.split_at_level = i;
	      float warp_iters_bw_levels = (max_col_len-(levels.cum_Rows[i]
							 -levels.cum_Rows[policy.split_at_level]))/WARP_SIZE;	    
	      if(warp_iters_bw_levels > WARP_COLWISE_ITERATIONS_THRESHOLD)
		{
		  {
		    splitGPU = 4;	
		  }
		}
	      else
		{
		  splitGPU = 0;
		  //break;
		}
	    }            
	}

      for(int i = policy.split_at_level; i < levels.nlevel; i++)
	{
	  warp_iters_cols = levels.max_levelColLength[i]/WARP_SIZE;
	  if(warp_iters_cols > WARP_COLWISE_ITERATIONS_THRESHOLD)
	    {
	      int warp_iters_bw_levels = (max_col_len - 
					  (levels.cum_Rows[i] - levels.cum_Rows[max_cols_level_indx]))/WARP_SIZE;
	      if(warp_iters_cols > warp_iters_bw_levels)
		{
		  policy.split_at_level = i;
		  //max_cols_level_indx = i;    // ??
		  splitGPU = 4;
	      }
	    }
	}
    }
  else
    {
      policy.split_at_level = levels.nlevel-1;
    }

    if((splitGPU == 0) || (policy.split_at_level == levels.nlevel-1))
      {
	splitGPU = 0;
	printf("GPU-GPU DAG slicing failure\n");
      }
    else
      {
	policy.sptrsv_method1 = ELMC;
        policy.sptrsv_method2 = ELMC;
        policy.spmv_platform = SPMV_GPU;
        policy.spmv_method = SPMV_CSR;
	policy.description=4;
      }    
    
    printf("SPLIT AT LEVEL: %d\n", policy.split_at_level);

    return splitGPU;
  
}
int SILU::GPUGPUSplit_old(csr_matrix *tri_mat)
{
    int splitGPU = 0;
    int maxRows = levels.num_levelRows[0];
    int max_rows_level_indx = 0;
    int max_cols_level_indx = 0;
    int foundSplit = 0;
    float warp_iters_cols = levels.max_levelColLength[0]/WARP_SIZE;
    int max_col_iters = 0;// warp_iters_cols;
    float warp_iters_rows = levels.max_levelRowLength[0]/WARP_SIZE;
    int max_row_iters = 0; //warp_iters_rows;
    policy.split_at_level=levels.nlevel;
    // Set split level at level with maximum max column length
    for(int i = 0; i < levels.nlevel; i++)
    {
        //warp_iters_rows = levels.max_levelRowLength[i]/WARP_SIZE;
        //if(warp_iters_rows > WARP_ROWISE_ITERATIONS_THRESHOLD)
        //{
        //    max_row_iters = levels.max_levelRowLength[i]/WARP_SIZE;
         //   max_rows_level_indx = i;
            //printf("HERERERE: %d\n", i);
        //}
        warp_iters_cols = levels.max_levelColLength[i]/WARP_SIZE;
	//printf("warp_iters_cols:%d\n", warp_iters_cols);
        if(warp_iters_cols > WARP_COLWISE_ITERATIONS_THRESHOLD && warp_iters_cols > max_col_iters)
        {
           policy.split_at_level = i;
           max_col_iters = levels.max_levelColLength[i]/WARP_SIZE;
           max_cols_level_indx = i;
           splitGPU=4; 
        }
	//printf("GPU-GPU split at level: %d\n", policy.split_at_level);
    }
  


    // Check if maximum row length level is beyond selected split point
    // Adjust split level accordingly
    for(int i = policy.split_at_level; i < levels.nlevel; i++)
    {
        warp_iters_rows = levels.max_levelRowLength[i]/WARP_SIZE;
        if(warp_iters_rows > WARP_ROWISE_ITERATIONS_THRESHOLD)
        {
            float warp_iters_bw_levels = (levels.cum_Rows[i]-levels.cum_Rows[policy.split_at_level])/WARP_SIZE;
            //printf("CUM ROWS: %f, %d\n", warp_iters_bw_levels,(levels.cum_Rows[i]-levels.cum_Rows[policy.split_at_level]));
            if(warp_iters_bw_levels > WARP_COLWISE_ITERATIONS_THRESHOLD)
            {
                if(warp_iters_bw_levels > 0.5 * max_col_iters)
                {
                    splitGPU = 0;
                    policy.split_at_level = levels.nlevel-1;
                    break;
                }
                else
                {
                    splitGPU = 4;
                    policy.split_at_level = i;
                    //printf("Why not here?\n");
                }
            }
            else
            {
                policy.split_at_level = i;
                splitGPU = 4;
            }
        }
        // if(policy.split_at_level < max_rows_level_indx)
        // {
        //     //printf("Was here\n");
        //     int warp_iters_bw_levels = (levels.cum_Rows[max_rows_level_indx] - levels.cum_Rows[policy.split_at_level])/WARP_SIZE;
        //     printf("%d, ", warp_iters_bw_levels);
        //     if(warp_iters_bw_levels > WARP_COLWISE_ITERATIONS_THRESHOLD)
        //     {
        //         if(warp_iters_bw_levels > 0.5 * max_col_iters)
        //         {
        //             //splitGPU = 0;
        //             break;
        //             policy.split_at_level = levels.nlevel;
        //         }
        //         else
        //         {
        //             policy.split_at_level = max_rows_level_indx; 
        //             splitGPU = 1;
        //         }
        //     }
        //     else
        //     {

        //     }
        // }    
    }
    
    
    if(splitGPU==4)
    {
        for(int i = policy.split_at_level; i < levels.nlevel; i++)
        {
            warp_iters_cols = levels.max_levelColLength[i]/WARP_SIZE;
            if(warp_iters_cols > WARP_COLWISE_ITERATIONS_THRESHOLD)
            {   
                int warp_iters_bw_levels = (levels.cum_Rows[i] - levels.cum_Rows[max_cols_level_indx])/WARP_SIZE;
                if(warp_iters_cols > warp_iters_bw_levels)
                {
                    policy.split_at_level = i;
                }
            }
        }
        //policy.split_at_level = i;
        if(policy.split_at_level == levels.nlevel-1)
        {
            splitGPU = 0;
            printf("Split level happens to be the last level\n");
        }    
        policy.sptrsv_method1 = ELMC;
        policy.sptrsv_method2 = ELMC;
        policy.spmv_platform = SPMV_GPU;
        policy.spmv_method = SPMV_CSR;
	policy.description=4;
	//printf("IWASFDFSDFSDF\n");
    }
    printf("SPLIT AT LEVEL: %d\n", policy.split_at_level);

    return splitGPU;
}

/*
int SILU::GPUSplitORUnified_old(csr_matrix *tri_mat)
{
    int splitGPU = 0;
    int maxRows = levels.num_levelRows[0];
    int maxRowIndx = 0;
    int foundSplit = 0;
    for(int i = 1; i < levels.nlevel; i++)
    {
        if(levels.num_levelRows[i] > maxRows)
        {
            maxRows = levels.num_levelRows[i];
            maxRowIndx = i;
        }    
    }
    for(int i = 0; i < levels.nlevel; i++)
    {
      //float col_len_as_perc_of_mat_size = (levels.max_levelColLength[i]/(float)tri_mat->m) * 100.0;
      //float row_len_as_perc_of_mat_size = (maxRows/(float)tri_mat->m) * 100.0;
      //printf("printing i:%d\n", i);
	int warp_column_wise_iterations_sptrsv = ceil(levels.max_levelColLength[i]/WARP_SIZE);
	if(warp_column_wise_iterations_sptrsv > WARP_COLWISE_ITERATIONS_THRESHOLD)
	  {
	    policy.split_at_level = i;
	    policy.sptrsv_method1 = ELMC;
	    policy.sptrsv_method2 = ELMC;
	    policy.spmv_platform = SPMV_GPU;
	    policy.spmv_method = SPMV_CSR;
	    splitGPU=1;
	    for(int j = i; j < levels.nlevel; j++)
	      {
		int warp_row_wise_iterations_spmv = ceil(levels.max_levelRowLength[j]/WARP_SIZE);
		if(warp_row_wise_iterations_spmv > WARP_ROWISE_ITERATIONS_THRESHOLD)
		  {
		    policy.split_at_level = j;            
		  }
	      }

	  } 
	//int warp_row_wise_iterations_spmv = ceil(levels.max_level
        /*
	if((col_len_as_perc_of_mat_size > MAX_COL_LENGTH_THRESHOLD) &&
           (row_len_as_perc_of_mat_size < MAX_ROW_LENGTH_THRESHOLD))           
        {
            policy.split_at_level = i;
            policy.sptrsv_method1 = ELMC;
            policy.sptrsv_method2 = ELMC;
            policy.spmv_platform = SPMV_GPU;
            policy.spmv_method = SPMV_CSR;
            splitGPU = 1;
            break;
        }
        else
        {
            if((col_len_as_perc_of_mat_size > MAX_COL_LENGTH_THRESHOLD) &&
           (row_len_as_perc_of_mat_size > MAX_ROW_LENGTH_THRESHOLD))
            {
                
                splitGPU = 0;
                policy.sptrsv_method1 = MKL;
                policy.sptrsv_method2 = MKL;
                policy.spmv_platform = SPMV_CPU;
                policy.spmv_method = SPMV_CSR;
                policy.split_at_level = levels.nlevel-1;
                break;
            }
            else
            {
                policy.sptrsv_method1 = ELMC;
                policy.sptrsv_method2 = ELMC;
                policy.spmv_platform = SPMV_GPU;
                policy.spmv_method = SPMV_CSR;
                policy.split_at_level = levels.nlevel-1;
                //splitGPU=0;
            }
            
	    }*/
//      if(splitGPU==1)
//        break;
//    }
//    return splitGPU;
//}



int SILU::FindSplitPoint(csr_matrix *tri_mat)
{
    int retCode = 0;
    int midLevel = (levels.nlevel%2==0)?levels.nlevel/2-1:(levels.nlevel-1)/2-1;
    int start_level = 0;
    int split_level = start_level;
    int end_level = 0;
    int consec_cpu_friendly_levels = 0;
    int cpu_friendly_levels_in_split = 0;
    int gpu_friendly_levels_in_split = 0;
    int gpu_friendly_rows_in_split = 0;
    float gpu_friendly_rows_in_split_perc = 0.0;
    float cpu_friendly_lvl_in_split_perc = 0.0;
    float rows_percentage_mid_level = (levels.cum_Rows[midLevel] / (float)tri_mat->m) * 100.0;
    int inc_j = 1; 
    if(rows_percentage_mid_level >= 0.5)
    {
        // Search split point on the left
        for(int i = start_level; i < levels.nlevel; i++)
        {
            //printf("Level rows: %d\n", levels.num_levelRows[i]);
            if(levels.num_levelRows[i] < CPU_FRIENDLY_MAX_ROWS_THRESHOLD)
            {
                consec_cpu_friendly_levels++;
                //printf("Matdfsdf\n");
                cpu_friendly_levels_in_split++;
                cpu_friendly_lvl_in_split_perc = cpu_friendly_levels_in_split/(float)policy.cpu_friendly_lvls * 100.0;
                if(cpu_friendly_lvl_in_split_perc > MAX_CPU_FRIENDLY_LVLS_IN_SPLIT_PERC)
                {
                    printf("CPU levels in split exceed limit\n");
                    cpu_friendly_levels_in_split -= (consec_cpu_friendly_levels);
                    cpu_friendly_lvl_in_split_perc = cpu_friendly_levels_in_split/(float)policy.cpu_friendly_lvls * 100.0;
                    retCode = 0;
                    break;
                }                
            }
            else
            {
                gpu_friendly_levels_in_split++;
                gpu_friendly_rows_in_split += levels.num_levelRows[i];
                split_level = split_level + consec_cpu_friendly_levels + 1;
                gpu_friendly_rows_in_split_perc = (gpu_friendly_rows_in_split /(float)policy.gpu_friendly_rows) * 100.0;
                //printf("DFD: %d, %d, %d\n",gpu_friendly_levels_in_split, consec_cpu_friendly_levels, split_level);
                consec_cpu_friendly_levels = 0;
                
                if(gpu_friendly_levels_in_split == policy.gpu_friendly_lvls)
                {
                    //retCode = 1;
                    break;
                }
            }
    }
    printf("Split at level: %d, CPU friendly levels: %f, GPU friendly rows:%f\n", split_level-1, cpu_friendly_lvl_in_split_perc, gpu_friendly_rows_in_split_perc);
    if((cpu_friendly_lvl_in_split_perc <= MAX_CPU_FRIENDLY_LVLS_IN_SPLIT_PERC) &&
       (gpu_friendly_rows_in_split_perc > MIN_GPU_FRIENDLY_ROWS_IN_SPLIT_PERC))
    {
      
      policy.description=3;
        retCode = 1;
    }

    }
    else
    {
        // Search split point  on the right
        start_level = midLevel;
        for(int i = start_level; i < levels.nlevel; i++)
        {
            if(levels.num_levelRows[i] >= CPU_FRIENDLY_MAX_ROWS_THRESHOLD)
            {
                inc_j = 0;
            }
            //else
            {
                
                for(int j = i; j < levels.nlevel; j++)
                {                                        
                    if(levels.num_levelRows[i] < CPU_FRIENDLY_MAX_ROWS_THRESHOLD)
                    {
                        if(inc_j == 1)
                        {
                            j++;
                            //split_level = j;
                        }
                        else
                        {
                            consec_cpu_friendly_levels++;
                            cpu_friendly_levels_in_split++;
                            cpu_friendly_lvl_in_split_perc = cpu_friendly_levels_in_split/(float)policy.cpu_friendly_lvls * 100.0;
                        }
                        //consec_cpu_friendly_levels++;
                        //printf("Matdfsdf\n");
                        //cpu_friendly_levels_in_split++;
                        
                        if(cpu_friendly_lvl_in_split_perc > MAX_CPU_FRIENDLY_LVLS_IN_SPLIT_PERC)
                        {
                            printf("CPU levels in split exceed limit (2):%d\n",j);
                            cpu_friendly_levels_in_split-=(consec_cpu_friendly_levels+1);
                            cpu_friendly_lvl_in_split_perc = cpu_friendly_levels_in_split/(float)policy.cpu_friendly_lvls * 100.0;
                            retCode = 0;
                            break;
                        }
                    }
                    else
                    {
                        inc_j = 0;
                        gpu_friendly_levels_in_split++;
                        gpu_friendly_rows_in_split += levels.cum_Rows[j];
                        consec_cpu_friendly_levels = 0;
                        gpu_friendly_rows_in_split_perc = (gpu_friendly_rows_in_split /(float)policy.gpu_friendly_rows) * 100.0;
                        if(gpu_friendly_levels_in_split == policy.gpu_friendly_lvls)
                        {
                            //retCode = 1;
                            break;
                        }

                    }
                }
                if((cpu_friendly_lvl_in_split_perc <= MAX_CPU_FRIENDLY_LVLS_IN_SPLIT_PERC) &&
                   (gpu_friendly_rows_in_split_perc > MIN_GPU_FRIENDLY_ROWS_IN_SPLIT_PERC))
                {
                    //printf("Found level: %d\n", i);
                    split_level = i;
                    retCode = 2;
                    break;
                }                    
            }            
        }

    }
    //printf("retCode: %d\n", retCode);
    if(retCode == 1)
    {
        policy.sptrsv_method1 = ELMC;
        policy.sptrsv_method2 = MKL;
        policy.spmv_method = SPMV_CSR;
        policy.spmv_platform = SPMV_GPU;
        policy.split_at_level = split_level-1;
        printf("Split at level(GPU-CPU): %d\n", policy.split_at_level);
	policy.description=3;
    }
    if(retCode == 2)
    {
        policy.sptrsv_method1 = MKL;
        policy.sptrsv_method2 = ELMC;
        policy.spmv_method = SPMV_CSR;
        policy.spmv_platform = SPMV_GPU;
        policy.split_at_level = split_level-1;
        printf("Split at level (CPU-GPU): %d\n", policy.split_at_level);   
	policy.description=2;
    }
    return retCode;
}

int SILU::HTS_ELMC_LogRes()
{
  int retCode = 0;
  int total_rows = levels.cum_Rows[levels.nlevel-1];
  int total_nnzs = levels.cum_Nnzs[levels.nlevel-1];
  float parallelism = total_rows/(float)(levels.nlevel);
  float mean_row_length = total_nnzs/(float)total_rows;
  int max_col_length = 0;
  double p = 0.0, pred=0.0;
  for(int i = 0; i < levels.nlevel; i++)
    {
      if(levels.max_levelColLength[i] > max_col_length)
	{
	  max_col_length=levels.max_levelColLength[i];
	}
    }
  p = -2.00493726 + (2.44439384e-02 * levels.nlevel) - (2.39147474e-07 * parallelism) 
    - (3.12534312e-02 * mean_row_length) + (9.41692274e-04 * max_col_length); 
  p = exp(-p);
  try
    {
      if(p == HUGE_VAL || p == -HUGE_VAL)
	throw "Overflow exception";
      pred = 1/(1+p);
    }
  catch(const char* e)
    {
      pred=0;
    }

  if (pred > 0.5)
    {
      retCode = 6;  // HTS
    }
  else
    {
      retCode = 1; // ELMC
    }
  
  return retCode;
}

int SILU::MKL_HTS_LogReg()
{
  int total_nnzs = levels.cum_Nnzs[levels.nlevel-1];
  int total_rows = levels.cum_Rows[levels.nlevel-1];
  float parallelism = total_rows/(float)(levels.nlevel);
  float mean_row_length = total_nnzs/(float)total_rows;
  int max_row_length = 0;
  double p = 0.0, pred=0.0;
  int MKL_HTS = 0;
  for(int i = 0; i < levels.nlevel; i++)
    {
      if(levels.max_levelRowLength[i] > max_row_length)
	{
	  max_row_length=levels.max_levelRowLength[i];
	}
    }
  printf("parallelism=%f, mean_row_length=%f, max_row_length=%d\n", parallelism, mean_row_length, max_row_length);
  // Logistic regression for V100-Gold
  p = 1.22003133 - (2.55600460e-02 * parallelism) - (1.06753910e-02 * mean_row_length) - (3.30841013e-05 * max_row_length);
  // Logistic regression for G1080-E5
  //p= 0.33719131 - (0.01352607 * parallelism) - (0.01403581 * mean_row_length) - (0.00017559 * max_row_length);  
  p = exp(-p);
  try
    {
      if(p == HUGE_VAL || p == -HUGE_VAL)
	throw "Overflow exception";
      pred = 1/(1+p);
    }
  catch(const char* e)
    {
      pred=0;
    }
  if (pred > 0.5)
    {
      MKL_HTS = 5;  // MKL
    }
  else
    {
      MKL_HTS = 6; // HTS
    }
  return MKL_HTS;
}

int SILU::CPUCPUSplit(csr_matrix *tri_mat, int hts_mkl)
{
  int retCode = 0; // 5=MKL, 6=HTS, 3=ELMC, unsuitable=4
  if(hts_mkl==1) // HTS-MKL classifier
    {
      printf("Testing if suitable for MKL-HTS\n");
      printf("CPU friendly levels:levels::%0.1f:%d\n", policy.cpu_friendly_lvl_perc, levels.nlevel);
      // Check whether suitable for HTS-MKL classification 
      if(policy.cpu_friendly_lvl_perc > MIN_LVL_PERC_FOR_CPU_EXEC ||
	 levels.nlevel > MAX_LEVELS)
	{
	  // Run HTS-MKL classifier
	  retCode = MKL_HTS_LogReg();
	  printf("Running MKL-HTS logistic regression, %d\n", retCode);
	}
      else
	{
	  printf("Not suitable for MKL-HTS selection\n");
	  retCode=0;
	} 
    }
  else
    { // Run HTS-ELMC classifier
      
      retCode = HTS_ELMC_LogRes();
    }
  return retCode;
}

int SILU::analyse_dag_revised(csr_matrix *tri_mat)
{
#ifndef AUTO_TESTING
  int doSplit = CPUGPUSplit(tri_mat);
  int retCode = 0;
  if(doSplit==1)
    {
      printf("Suitable for split\n");
      int splitStatus = FindSplitPoint(tri_mat);
      // If finding split point is a success
      if(splitStatus==1 || splitStatus==2)
	{
	  printf("Finding split-point for CPU-GPU is a success\n");
	  retCode=3;
	}
      else
	{
	  // Check if suitable for HTS-MKL
	  printf("CPU-CPU test after CPU-GPU split-point selection failure\n");
	  retCode=CPUCPUSplit(tri_mat, 1);
	  if(retCode==0)
	    {
	      // Try GPU-GPU split
	      printf("Trying GPU-GPU split\n");
	      retCode=GPUGPUSplit(tri_mat);
	      printf("GPU-GPU retCode:%d\n", retCode);
	      if(retCode==0)
		{
		  // Try HTS-ELMC
		  printf("Running HTS-ELMC logistic regression\n");
		  retCode=CPUCPUSplit(tri_mat,0);
		  printf("HTS-ELMC retCode:%d\n",retCode);
		  if(retCode==6)
		    {
		      policy.sptrsv_method1 = MKL;   // This is reserved for HTS, MKL as stub
		      policy.spmv_method = SPMV_CSR;
		      policy.spmv_platform = SPMV_CPU;
		      policy.sptrsv_method2 = MKL;  // HTS
		      policy.split_at_level = levels.nlevel-1;
		      policy.description=retCode;
		    }
		  else
		    {
		      policy.sptrsv_method1 = ELMC;   
		      policy.spmv_method = SPMV_CSR;
		      policy.spmv_platform = SPMV_GPU;
		      policy.sptrsv_method2 = ELMC;  
		      policy.split_at_level = levels.nlevel-1;
		      policy.description=retCode;		      
		    }
		} 
	    }
	  else
	    {
	      policy.sptrsv_method1 = MKL;   // This is reserved for HTS, MKL as stub
	      policy.spmv_method = SPMV_CSR;
	      policy.spmv_platform = SPMV_CPU;
	      policy.sptrsv_method2 = MKL;  // HTS
	      policy.split_at_level = levels.nlevel-1;
	      policy.description=retCode;
	    }
	}
    }
  else
    {
      printf("Not suitable for CPU-GPU split, %d\n", doSplit);
      if(doSplit==0) // Run HTS-ELMC classifier  
	{
	  printf("CPU-CPU test after CPU-GPU basic test failure\n");
	  retCode=CPUCPUSplit(tri_mat,1);
	  //printf("retCode:%d\n", retCode);
	  if(retCode==0)
	    {
	      // Try GPU-GPU split
	      printf("Trying GPU-GPU split\n");
	      retCode=GPUGPUSplit(tri_mat);
	      printf("GPU-GPU retCode:%d\n", retCode);
	      if(retCode==0)
		{
		  printf("Running HTS-ELMC logistic regression\n");
		  // Try HTS-ELMC
		  retCode=CPUCPUSplit(tri_mat,0);
		  printf("HTS-ELMC retCode:%d\n",retCode);
		  if(retCode==6)
		    {
		      policy.sptrsv_method1 = MKL;   // This is reserved for HTS, MKL as stub
		      policy.spmv_method = SPMV_CSR;
		      policy.spmv_platform = SPMV_CPU;
		      policy.sptrsv_method2 = MKL;  // HTS
		      policy.split_at_level = levels.nlevel-1;
		      policy.description=retCode;
		    }
		  else
		    {
		      policy.sptrsv_method1 = ELMC;   
		      policy.spmv_method = SPMV_CSR;
		      policy.spmv_platform = SPMV_GPU;
		      policy.sptrsv_method2 = ELMC;  
		      policy.split_at_level = levels.nlevel-1;
		      policy.description=retCode;		      
		    }
		}
	      else
		{    // GPU-GPU split is a success
		  printf("GPU-GPU split successful. Split level: %d\n", policy.split_at_level);
		}
	    }
	  else
	    {
	      if(retCode=5)
		{
		  policy.sptrsv_method1 = MKL;   
		  policy.spmv_method = SPMV_CSR;
		  policy.spmv_platform = SPMV_CPU;
		  policy.sptrsv_method2 = MKL;  
		  policy.split_at_level = levels.nlevel-1;
		  policy.description=retCode;
		}
	      else   // HTS
		{
		  policy.sptrsv_method1 = MKL;   // This is reserved for HTS, MKL as stub
		  policy.spmv_method = SPMV_CSR;
		  policy.spmv_platform = SPMV_CPU;
		  policy.sptrsv_method2 = MKL;  // HTS
		  policy.split_at_level = levels.nlevel-1;
		  policy.description=retCode;
		}
	    }
	}
      if(doSplit==2)  // Number of levels threshold < MIN_LEVELS
	{
	  // Try ELMC GPU-GPU
	  //printf("I WAS HERE\n");
	  retCode=GPUGPUSplit(tri_mat);
	  if(retCode==0)
	    {
	      retCode=CPUCPUSplit(tri_mat,0);
	      if(retCode==6)
		{
		  policy.sptrsv_method1 = MKL;   // This is reserved for HTS, MKL as stub
		  policy.spmv_method = SPMV_CSR;
		  policy.spmv_platform = SPMV_CPU;
		  policy.sptrsv_method2 = MKL;  // HTS
		  policy.split_at_level = levels.nlevel-1;
		  policy.description=retCode;
		}
	      else
		{
		  policy.sptrsv_method1 = ELMC;   
		  policy.spmv_method = SPMV_CSR;
		  policy.spmv_platform = SPMV_GPU;
		  policy.sptrsv_method2 = ELMC;  
		  policy.split_at_level = levels.nlevel-1;
		  policy.description=retCode;		      
		}
	    }
	  else
	    { // GPU-GPU split is a success
	      printf("GPU-GPU split successful. Split level: %d\n", policy.split_at_level);
	    }
	}      
    }

#endif
  policy.rhs = 1;
  policy.splitting_policy = DAG_SPLIT;   // DAG_SPLIT, MATRIX_SPLIT
  if(policy.splitting_policy == DAG_SPLIT)
    {
      if(ELMC_only)
	{
	  policy.sptrsv_method1 = ELMC;
	  policy.spmv_method = SPMV_CSR;
	  policy.spmv_platform = SPMV_GPU;
	  policy.sptrsv_method2 = ELMC;
	  policy.split_at_level = levels.nlevel-1;
	  policy.description=1;
	}

      if(MKL_only)
	{
	  //printf("I was here\n");
	  policy.sptrsv_method1 = MKL;
	  policy.spmv_method = SPMV_CSR;
	  policy.spmv_platform = SPMV_CPU;
	  policy.sptrsv_method2 = MKL;
	  policy.split_at_level = levels.nlevel-1;
	  policy.description=0;    
	}
    }
}

int SILU::analyse_dag(csr_matrix *tri_mat)
{	
    //printf("here\n");
#ifndef AUTO_TESTING
    int doSplit = PreProcess(tri_mat);
    //printf("Was here: %d\n", doSplit);
    if(doSplit==1)
    {
        printf("Suitable for split\n");
        int splitStatus = FindSplitPoint(tri_mat);
        if(splitStatus == 0)
        {
            printf("Later turned out not Suitable for split\n");
            int platform = DecidePlatform(tri_mat);
            if(platform == 1)
            {
                printf("Suitable for execution on CPU(1)\n");
                policy.sptrsv_method1 = MKL;   // This is reserved for HTS, MKL as stub
                policy.spmv_method = SPMV_CSR;
                policy.spmv_platform = SPMV_CPU;
                policy.sptrsv_method2 = MKL;  // HTS
                policy.split_at_level = levels.nlevel-1;
		policy.description=5;
            }
            else
            {
                printf("Suitable for execution on the GPU (1)\n");
                //int gpuExecMode = GPUSplitORUnified(tri_mat);
                int gpuExecMode=GPUGPUSplit(tri_mat);
		if(gpuExecMode == 1)
                {
		  policy.description=4;
                    printf("Split GPU-GPU execution: Split-level: %d\n",policy.split_at_level);
                    //policy.sptrsv_method1 = ELMC;
                    //policy.spmv_method2 = SPMV_CPU;
                    //policy.sptrsv_method2 = MKL
                    //policy.split_at_level = levels.nlevel-1;
                }
                else
                {
                    policy.sptrsv_method1 = ELMC;
                    policy.spmv_method = SPMV_CSR;
                    policy.spmv_platform = SPMV_GPU;
                    policy.sptrsv_method2 = ELMC;
                    policy.split_at_level = levels.nlevel-1;
		    policy.description=1;
                }

            }
        }
    } 
    else
    {
      switch(doSplit)
      {
      case 0:
	policy.sptrsv_method1 = ELMC;
	policy.spmv_method=SPMV_CSR;
	policy.spmv_platform=SPMV_GPU; 
	policy.sptrsv_method2=ELMC;
	policy.split_at_level = levels.nlevel-1;
	policy.description=1;
	break;

      case 2:
	policy.sptrsv_method1 = MKL;    // This is reserved for HTS, MKL as stub
	policy.spmv_method=SPMV_CSR;
	policy.spmv_platform=SPMV_GPU; 
	policy.sptrsv_method2=MKL;
	policy.split_at_level = levels.nlevel-1;
	printf("Policy LOG_REG: %d\n", MKL_HTS_LogReg());
	policy.description=MKL_HTS_LogReg();
	break;

      case 3:
	printf("Testing for GPU split or unified\n");
	//int gpuExecMode = GPUSplitORUnified(tri_mat);
	int gpuExecMode=GPUGPUSplit(tri_mat);
	//printf("gpuExecMode: %d\n", gpuExecMode);
	if(gpuExecMode == 1)
	  {
	    policy.description=4;
	    printf("Split GPU-GPU execution: Split-level: %d\n",policy.split_at_level);
	  }
	else
	  {
	    policy.sptrsv_method1 = ELMC;
	    policy.spmv_method = SPMV_CSR;
	    policy.spmv_platform = SPMV_GPU;
	    policy.sptrsv_method2 = ELMC;
	    policy.split_at_level = levels.nlevel-1;
	    policy.description=1;
	  }
	break;
      }

	/*if(dosplit==2)
	{
	  printf("Suitable for execution on the CPU using HTS\n");
	  policy.sptrsv_method1 = HTS;
	  policy.spmv_method=SPMV_CSR;
	  policy.spmv_platform=SPMV_CPU; 
	  policy.sptrsv_method2=HTS;
	  policy.split_at_level = levels.nlevel-1;
	}  
      //printf("Not Suitable for split\n");
        int platform = DecidePlatform(tri_mat);
        if(platform == 1)
        {
            printf("Suitable for execution on CPU(1)\n");
            policy.sptrsv_method1 = MKL;
            policy.spmv_method = SPMV_CSR;
            policy.spmv_platform = SPMV_CPU;
            policy.sptrsv_method2 = MKL;
            policy.split_at_level = levels.nlevel-1;
        }
        else
        {
            printf("Suitable for execution on the GPU (2)\n");
            int gpuExecMode = GPUSplitORUnified(tri_mat);
            if(gpuExecMode == 1)
            {
                printf("Split GPU-GPU execution: Split-level: %d\n",policy.split_at_level);
            }

        }
	*/

    }
    //printf("Passed through it\n");
    //while(1);

	//policy.sptrsv_method1 = ELMC;
	//policy.sptrsv_method2 = MKL;
	//policy.spmv_method = SPMV_CSR;
	//policy.spmv_platform = SPMV_GPU;

#endif
    policy.rhs = 1;
    policy.splitting_policy = DAG_SPLIT;   // DAG_SPLIT, MATRIX_SPLIT
    
    if(policy.splitting_policy == DAG_SPLIT)
    {
#ifndef AUTO_TESTING
        //policy.split_at_level = levels.nlevel-1;
        //policy.split_at_level = 3;
	//policy.split_at_level = 287;
#endif

#ifdef ELMC_ONLY
    policy.sptrsv_method1 = ELMC;
    policy.spmv_method = SPMV_CSR;
    policy.spmv_platform = SPMV_GPU;
    policy.sptrsv_method2 = ELMC;
    policy.split_at_level = levels.nlevel-1;
    policy.description=1;
#endif
if(ELMC_only)
{
    policy.sptrsv_method1 = ELMC;
    policy.spmv_method = SPMV_CSR;
    policy.spmv_platform = SPMV_GPU;
    policy.sptrsv_method2 = ELMC;
    policy.split_at_level = levels.nlevel-1;
    policy.description=1;
}

if(MKL_only)
{
    printf("I was here\n");
    policy.sptrsv_method1 = MKL;
    policy.spmv_method = SPMV_CSR;
    policy.spmv_platform = SPMV_CPU;
    policy.sptrsv_method2 = MKL;
    policy.split_at_level = levels.nlevel-1;
    policy.description=0;    
} 

#ifdef MKL_ONLY
    policy.sptrsv_method1 = MKL;
    policy.spmv_method = SPMV_CSR;
    policy.spmv_platform = SPMV_CPU;
    policy.sptrsv_method2 = MKL;
    policy.split_at_level = levels.nlevel-1;
    policy.description=0;
#endif 

#ifdef CUSPARSE_ONLY
    policy.sptrsv_method1 = CUSPARSE_V2;
    policy.spmv_method = SPMV_CSR;
    policy.spmv_platform = SPMV_GPU;
    policy.sptrsv_method2 = CUSPARSE_V2;
    policy.split_at_level = levels.nlevel-1;
    policy.description=1;
#endif
        //policy.split_at_level = 1242;
        //policy.split_at_level = 26;
        //policy.sptrsv_method1 = ELMC;
	//policy.sptrsv_method2 = ELMC;
	//policy.spmv_method = SPMV_CSR;
	//policy.spmv_platform = SPMV_GPU;
        //policy.nosplit = 1;
        printf("Split at level=%d, Percent rows=%f%%, Percent nnzs=%f%%\n", 
                policy.split_at_level,
                levels.cum_Rows[policy.split_at_level]/float(tri_mat->m)*100.0,
                levels.cum_Nnzs[policy.split_at_level]/float(tri_mat->nnz)*100.0);           
    }
    else
    {
        policy.split_at_row = tri_mat->m;
        //policy.split_at_row = 1.7*1e7;
        // policy.split_at_row = 0.65 * 1e8;    // For matrix 2803
        //policy.split_at_row = 206500;  // 2376_2265
        //policy.split_at_row = 1270432;  // 2265_2376
        //policy.split_at_row = 69365;
        //policy.split_at_row = 27829;
        //policy.split_at_row = 440020;  // 1893_369 original split point
        //policy.split_at_row = 217918;
        //policy.nosplit = 1;    
        printf("Split at row=%d\n", policy.split_at_row);
    }

    if(policy.sptrsv_method1 == MKL)
    {
        policy.sptrsv1_platform = SPTRSV_CPU;
    }
    else
    {
        policy.sptrsv1_platform = SPTRSV_GPU;
    }
    if(policy.sptrsv_method2 == MKL)
    {
        policy.sptrsv2_platform = SPTRSV_CPU;
    }
    else
    {
        policy.sptrsv2_platform = SPTRSV_GPU;
    }
	
	return EXIT_SUCCESS;	
}

int SILU::findlevels_incremental(csr_matrix *tri_mat, std::string filename)
{
  printf("%s\n", filename);
  levels.level_of_row = (ind_type *)malloc((tri_mat->m) * sizeof(ind_type));
    levels.level_delta = (ind_type *)malloc((tri_mat->m) * sizeof(ind_type));
    levels.parallelism = (float *)malloc((tri_mat->m) * sizeof(ind_type));
    levels.nnz_per_row = (ind_type *)malloc((tri_mat->m) * sizeof(ind_type));
    levels.nnz_per_col = (ind_type *)malloc((tri_mat->m) * sizeof(ind_type));
    levels.avg_row_nnzs = (float *)malloc((tri_mat->m) * sizeof(float));
    levels.avg_col_nnzs = (float *)malloc((tri_mat->m) * sizeof(float));
    levels.avg_level_delta = (float *)malloc((tri_mat->m) * sizeof(float));
    levels.col_center = (float *)malloc((tri_mat->m) * sizeof(float));
    memset(levels.level_of_row, 0, (tri_mat->m) *sizeof(ind_type));
    memset(levels.parallelism, 0, (tri_mat->m) *sizeof(ind_type));

    csc_matrix A_csc;
    matrix_transpose(tri_mat, &A_csc);

    int prev_total_levels = 0;
    int current_total_levels = 0;
    int cum_row_nnzs = 0;
    int cum_col_nnzs = 0;
    float cum_level_delta = 0;
    printf("Rows: %d\n", tri_mat->m);
    for(int i = 0; i < tri_mat->m; i++)
    {
        int max_level = 0;
        levels.nnz_per_col[i] = A_csc.cscColPtr[i+1] - A_csc.cscColPtr[i];
        for(int j = tri_mat->csrRowPtr[i]; j < tri_mat->csrRowPtr[i+1]; j++)
        {
            if(tri_mat->csrColIdx[j] < i)
            {
                if(levels.level_of_row[tri_mat->csrColIdx[j]] > max_level)
                {
                    max_level = levels.level_of_row[tri_mat->csrColIdx[j]];
                }
            }            
        }
        if(current_total_levels < (max_level + 1))
        {
            current_total_levels = (max_level + 1);    
        }
        
        levels.level_of_row[i] = max_level + 1;
        levels.level_delta[i] = (current_total_levels - prev_total_levels) * 100.0;
        cum_level_delta += levels.level_delta[i];
        levels.avg_level_delta[i] = cum_level_delta/float(i+1);
        //levels.level_delta[i] = (current_total_levels); 
        levels.nlevel = current_total_levels;
        prev_total_levels = current_total_levels;
        levels.parallelism[i] = (float)(i+1)/current_total_levels;
        levels.nnz_per_row[i] = tri_mat->csrRowPtr[i+1] - tri_mat->csrRowPtr[i];
        cum_row_nnzs += levels.nnz_per_row[i];
        levels.avg_row_nnzs[i] = cum_row_nnzs/(float)(i+1);
        cum_col_nnzs += levels.nnz_per_col[i];
        levels.avg_col_nnzs[i] = cum_col_nnzs/(float)(i+1);
    }

    for(int i = 0; i < tri_mat->n; i++)
    {
        levels.col_center[i] = 0.0;
        for(int j = A_csc.cscColPtr[i]; j < A_csc.cscColPtr[i+1]; j++)
        {
            levels.col_center[i] += (float)A_csc.cscRowIdx[j]/levels.nnz_per_col[i]; 
        }
    }

    SaveMatrixFeatures(tri_mat,filename);
    // for(int i = 0; i < tri_mat->m; i++)
    // {
    //     printf("row[%d]=%d, delta[%d]=%d, parallelism=%.1f, nnzs_row[%d]=%d, nnz_col[%d]=%d\n",i, levels.level_of_row[i], i, 
    //                                                           levels.level_delta[i], levels.parallelism[i],
    //                                                           i, levels.nnz_per_row[i], i, levels.nnz_per_col[i]);
    // }    

    return EXIT_SUCCESS;
}

int SILU::findlevels_csr(csr_matrix *tri_mat)
{
	if(tri_mat->m != tri_mat->n)
	{
		printf("The matrix is not square. Exiting!\n");
		return EXIT_FAILURE;
	}

    csc_matrix A_csc;
    matrix_transpose(tri_mat, &A_csc);
    
    levels.levelPtr = (ind_type *)malloc((tri_mat->m+1) * sizeof(ind_type));
    levels.levelItem = (ind_type *)malloc((tri_mat->m)*sizeof(ind_type));
    levels.levelItemNewRowIdx = (ind_type *)malloc((tri_mat->m+1)*sizeof(ind_type));
    levels.num_levelRows = (ind_type *)malloc((tri_mat->m)*sizeof(ind_type));
    levels.num_levelNnzs = (ind_type *)malloc((tri_mat->m)*sizeof(ind_type));
    levels.cum_Rows = (ind_type *)malloc((tri_mat->m)*sizeof(ind_type));
    levels.cum_Nnzs = (ind_type *)malloc((tri_mat->m)*sizeof(ind_type));
    levels.num_indegree = (ind_type *)malloc((tri_mat->m)*sizeof(ind_type));
    levels.num_outdegree = (ind_type *)malloc((tri_mat->m)*sizeof(ind_type));
    levels.max_levelRowLength = (ind_type *)malloc((tri_mat->m)*sizeof(ind_type));
    levels.max_levelColLength = (ind_type *)malloc((tri_mat->m)*sizeof(ind_type));

    memset(levels.num_levelRows, 0, (tri_mat->m) *sizeof(ind_type));
    memset(levels.num_levelNnzs, 0, (tri_mat->m) *sizeof(ind_type));
    memset(levels.num_indegree, 0, (tri_mat->m) *sizeof(ind_type));
    memset(levels.num_outdegree, 0, (tri_mat->m) *sizeof(ind_type));
    memset(levels.max_levelRowLength, 0, (tri_mat->m) *sizeof(ind_type));
    memset(levels.max_levelColLength, 0, (tri_mat->m) *sizeof(ind_type));
    
    findlevels(tri_mat, &A_csc);
    
    return EXIT_SUCCESS;
}

int SILU::findlevels_csr(csr_matrix *tri_mat, level_info *levels)
{
    if(tri_mat->m != tri_mat->n)
    {
        printf("The matrix is not square. Exiting!\n");
        return EXIT_FAILURE;
    }

    csc_matrix A_csc;
    matrix_transpose(tri_mat, &A_csc);
    
    levels->levelPtr = (ind_type *)malloc((tri_mat->m+1) * sizeof(ind_type));
    levels->levelItem = (ind_type *)malloc((tri_mat->m)*sizeof(ind_type));
    levels->levelItemNewRowIdx = (ind_type *)malloc((tri_mat->m+1)*sizeof(ind_type));
    levels->num_levelRows = (ind_type *)malloc((tri_mat->m)*sizeof(ind_type));
    levels->num_levelNnzs = (ind_type *)malloc((tri_mat->m)*sizeof(ind_type));
    levels->cum_Rows = (ind_type *)malloc((tri_mat->m)*sizeof(ind_type));
    levels->cum_Nnzs = (ind_type *)malloc((tri_mat->m)*sizeof(ind_type));
    levels->num_indegree = (ind_type *)malloc((tri_mat->m)*sizeof(ind_type));
    levels->num_outdegree = (ind_type *)malloc((tri_mat->m)*sizeof(ind_type));
    levels->max_levelRowLength = (ind_type *)malloc((tri_mat->m)*sizeof(ind_type));
    levels->max_levelColLength = (ind_type *)malloc((tri_mat->m)*sizeof(ind_type));

    memset(levels->num_levelRows, 0, (tri_mat->m) *sizeof(ind_type));
    memset(levels->num_levelNnzs, 0, (tri_mat->m) *sizeof(ind_type));
    memset(levels->num_indegree, 0, (tri_mat->m) *sizeof(ind_type));
    memset(levels->num_outdegree, 0, (tri_mat->m) *sizeof(ind_type));
    memset(levels->max_levelRowLength, 0, (tri_mat->m) *sizeof(ind_type));
    memset(levels->max_levelColLength, 0, (tri_mat->m) *sizeof(ind_type));
    findlevels(tri_mat, &A_csc, levels);
    
    return EXIT_SUCCESS;
}


int SILU::findlevels(csr_matrix *mat_csr, csc_matrix *mat_csc)
{
	int *indegree = (int *) malloc(mat_csr->m * sizeof(ind_type));

	for(int i = 0; i < mat_csr->m; i++)
	{
		indegree[i] = mat_csr->csrRowPtr[i+1] - mat_csr->csrRowPtr[i];
        //printf("ind[%d]=%d\n",i,indegree[i]-1);
	}

	int ptr = 0;

	// Find root items
	levels.levelPtr[0] = 0;
	for(int i = 0; i < mat_csr->m; i++)
	{
		if(indegree[i] == 1)
		{
			levels.levelItem[ptr] = i;
			levels.levelItemNewRowIdx[i] = ptr;
			ptr++;
		}
	}
    
	levels.levelPtr[1] = ptr;
	levels.cum_Rows[0] = ptr;
	levels.cum_Nnzs[0] = ptr;

	int lvi = 1;
	int max_level_row_length = 0;
	int max_level_col_length = 0;
	int rL = 0;
	int cL = 0;
    int test = 0;
	while(levels.levelPtr[lvi-1] != mat_csr->m)
	{
	  max_level_row_length = 0;
	  max_level_col_length = 0;
		for(int i = levels.levelPtr[lvi-1]; i < levels.levelPtr[lvi]; i++)
		{
			//test++;
            int node = levels.levelItem[i];
			for(int j = mat_csc->cscColPtr[node]; j < mat_csc->cscColPtr[node+1]; j++)
			{
				int visit_node = mat_csc->cscRowIdx[j];
				indegree[visit_node]--;
				
				if(indegree[visit_node] == 1)
				{
					levels.levelItem[ptr] = visit_node;
					levels.levelItemNewRowIdx[visit_node] = ptr;
					ptr++;
				}
			}
			
			levels.num_levelNnzs[lvi-1]+=mat_csr->csrRowPtr[node+1]-mat_csr->csrRowPtr[node];

			rL = mat_csr->csrRowPtr[node+1]-mat_csr->csrRowPtr[node]-1;
			cL = mat_csc->cscColPtr[node+1]-mat_csc->cscColPtr[node]-1;
			levels.num_indegree[lvi-1]+=mat_csr->csrRowPtr[node+1]-mat_csr->csrRowPtr[node]-1;
			levels.num_outdegree[lvi-1]+=mat_csc->cscColPtr[node+1]-mat_csc->cscColPtr[node]-1;
			if(rL > max_level_row_length)
			  {
			    max_level_row_length = rL;
			    levels.max_levelRowLength[lvi-1] = max_level_row_length;
			  }
			if(cL > max_level_col_length)
			  {
			    max_level_col_length = cL;
			    levels.max_levelColLength[lvi-1] = max_level_col_length;
			  }
			    
			   levels.num_levelRows[lvi-1]++;
			   if(lvi > 1)
			     {
			       levels.cum_Rows[lvi-1] = levels.cum_Rows[lvi-2] + levels.num_levelRows[lvi-1];
			       levels.cum_Nnzs[lvi-1] = levels.cum_Nnzs[lvi-2] + levels.num_levelNnzs[lvi-1]; 							
			     }
			   
			   }

		  //test = levels.num_levelRows[lvi-1];
		  //printf("%d\n", test);
          lvi++;
			levels.levelPtr[lvi] = ptr;		
		}
        //printf("Rows covered: %d\n", test);

		
		levels.nlevel = lvi-1;
		//printf("total rows: %d\n", test);
        //for(int i = 0; i < levels.nlevel; i++)
        //    printf("%d\n", levels.num_levelNnzs[i]);
        //while(1);
        free(indegree);    
		
		return EXIT_SUCCESS;
	} 

/*
int SILU::findlevels(csr_matrix *mat_csr, csc_matrix *mat_csc)
{
    int *indegree = (int *) malloc(mat_csr->m * sizeof(ind_type));
    int max_level_row_length = 0;
    int max_level_col_length = 0;
    int rL = 0;
    int cL = 0;

    for(int i = 0; i < mat_csr->m; i++)
    {
        indegree[i] = mat_csr->csrRowPtr[i+1] - mat_csr->csrRowPtr[i];
        //printf("ind[%d]=%d\n",i,indegree[i]-1);
    }

    int ptr = 0;

    // Find root items
    levels.levelPtr[0] = 0;
    for(int i = 0; i < mat_csr->m; i++)
    {
        if(indegree[i] == 1)
        {
            rL = mat_csr->csrRowPtr[i+1]-mat_csr->csrRowPtr[i];
            cL = mat_csc->cscColPtr[i+1]-mat_csc->cscColPtr[i];
            if(rL > max_level_row_length)
              {
                max_level_row_length = rL;
                levels.max_levelRowLength[0] = max_level_row_length;
              }
            if(cL > max_level_col_length)
              {
                max_level_col_length = cL;
                levels.max_levelColLength[0] = max_level_col_length;
              }
            levels.num_levelRows[0]++;
            levels.num_levelNnzs[0]+=mat_csr->csrRowPtr[i+1]-mat_csr->csrRowPtr[i];
            levels.levelItem[ptr] = i;
            levels.levelItemNewRowIdx[i] = ptr;
            ptr++;
        }
    }
    
    levels.levelPtr[1] = ptr;
    levels.cum_Rows[0] = ptr;
    levels.cum_Nnzs[0] = ptr;
    //levels.max_levelRowLength[0]=1;

    int lvi = 1;
    
    
    while(levels.levelPtr[lvi-1] != mat_csr->m)
    {
      max_level_row_length = 0;
      max_level_col_length = 0;
        for(int i = levels.levelPtr[lvi-1]; i < levels.levelPtr[lvi]; i++)
        {
            //test++;
            int node = levels.levelItem[i];
            for(int j = mat_csc->cscColPtr[node]; j < mat_csc->cscColPtr[node+1]; j++)
            {
                int visit_node = mat_csc->cscRowIdx[j];
                indegree[visit_node]--;
                
                if(indegree[visit_node] == 1)
                {
                    levels.levelItem[ptr] = visit_node;
                    levels.levelItemNewRowIdx[visit_node] = ptr;
                    ptr++;
                }
            }
            
            levels.num_levelNnzs[lvi]+=mat_csr->csrRowPtr[node+1]-mat_csr->csrRowPtr[node];

            rL = mat_csr->csrRowPtr[node+1]-mat_csr->csrRowPtr[node];
            cL = mat_csc->cscColPtr[node+1]-mat_csc->cscColPtr[node]-1;
            levels.num_indegree[lvi]+=mat_csr->csrRowPtr[node+1]-mat_csr->csrRowPtr[node]-1;
            levels.num_outdegree[lvi]+=mat_csc->cscColPtr[node+1]-mat_csc->cscColPtr[node]-1;
            if(rL > max_level_row_length)
              {
                max_level_row_length = rL;
                levels.max_levelRowLength[lvi] = max_level_row_length;
              }
            if(cL > max_level_col_length)
              {
                max_level_col_length = cL;
                levels.max_levelColLength[lvi] = max_level_col_length;
              }
                
               levels.num_levelRows[lvi]++;
               if(lvi > 1)
                 {
                   levels.cum_Rows[lvi] = levels.cum_Rows[lvi-1] + levels.num_levelRows[lvi];
                   levels.cum_Nnzs[lvi] = levels.cum_Nnzs[lvi-1] + levels.num_levelNnzs[lvi];                           
                 }
               
               }
          
          lvi++;
            levels.levelPtr[lvi] = ptr;     
        }
        //printf("Rows covered: %d\n", test);
        
        
        levels.nlevel = lvi-1;
        
        printf("Printing, levels:%d\n", levels.nlevel);
        for(int i = 0; i < levels.nlevel; i++)
        {
            //printf("level rows\n");
            printf("%d\n ", levels.num_levelRows[i]);
        }
        while(1);
        free(indegree);    
        
        return EXIT_SUCCESS;
    }
    */

int SILU::findlevels(csr_matrix *mat_csr, csc_matrix *mat_csc, level_info *levels)
{
    int *indegree = (int *) malloc(mat_csr->m * sizeof(ind_type));

    for(int i = 0; i < mat_csr->m; i++)
    {
        indegree[i] = mat_csr->csrRowPtr[i+1] - mat_csr->csrRowPtr[i];
    }

    int ptr = 0;

    // Find root items
    levels->levelPtr[0] = 0;
    for(int i = 0; i < mat_csr->m; i++)
    {
        if(indegree[i] == 1)
        {
            levels->levelItem[ptr] = i;
            levels->levelItemNewRowIdx[i] = ptr;
            ptr++;
        }
    }
    
    levels->levelPtr[1] = ptr;
    levels->cum_Rows[0] = ptr;
    levels->cum_Nnzs[0] = ptr;

    int lvi = 1;
    while(levels->levelPtr[lvi-1] != mat_csr->m)
    {
        for(int i = levels->levelPtr[lvi-1]; i < levels->levelPtr[lvi]; i++)
        {
            int node = levels->levelItem[i];
            for(int j = mat_csc->cscColPtr[node]; j < mat_csc->cscColPtr[node+1]; j++)
            {
                int visit_node = mat_csc->cscRowIdx[j];
                indegree[visit_node]--;
                
                if(indegree[visit_node] == 1)
                {
                    levels->levelItem[ptr] = visit_node;
                    levels->levelItemNewRowIdx[visit_node] = ptr;
                    //printf("[%d]=%d\n",visit_node, ptr);
                    ptr++;
                }
            }
            
            levels->num_levelNnzs[lvi-1]+=mat_csr->csrRowPtr[node+1]-mat_csr->csrRowPtr[node];           
            levels->num_levelRows[lvi-1]++;
            if(lvi > 1)
            {
                levels->cum_Rows[lvi-1] = levels->cum_Rows[lvi-2] + levels->num_levelRows[lvi-1];
                levels->cum_Nnzs[lvi-1] = levels->cum_Nnzs[lvi-2] + levels->num_levelNnzs[lvi-1];                          
            }
            
        }

        lvi++;
        levels->levelPtr[lvi] = ptr;     
    }
    
    levels->nlevel = lvi-1;
    free(indegree);    

    return EXIT_SUCCESS;
}

// Set execution policy
// DAG_SPLIT or MATRIX_SPLIT

int SILU::SetExecutionPolicy(int splitting_policy)
{
	//policy.splitting_policy = splitting_policy;
	return EXIT_SUCCESS;
}

void SILU::PrintLevelRows(void)
{
	printf("Matrix levels=%d\n", levels.nlevel);
	for(int i = 0; i < levels.nlevel; i++)
	{
		for(int j = levels.levelPtr[i]; j < levels.levelPtr[i+1]; j++)
		{
			printf("%ld, ", levels.levelItem[j]);
		}
		printf("Rows1=%ld, Nnzs=%ld, cum_Rows=%ld, cum_Nnzs=%ld\n", 
			levels.num_levelRows[i], levels.num_levelNnzs[i], 
			levels.cum_Rows[i], levels.cum_Nnzs[i]);
	}
}

int SILU::matrix_transpose(csr_matrix *in_mat, csc_matrix *out_mat)
{
	
	out_mat->cscRowIdx = (ind_type *)malloc(in_mat->nnz * sizeof(ind_type));
	out_mat->cscColPtr = (ind_type *)malloc((in_mat->n+1) * sizeof(ind_type));
	out_mat->cscVal = (val_type *)malloc((in_mat->nnz) * sizeof(val_type));
    
    out_mat->n = in_mat->n;
    out_mat->m = in_mat->m;
    out_mat->nnz = in_mat->nnz;
    
	memset(out_mat->cscColPtr, 0, (in_mat->n+1) *sizeof(ind_type));

	for(int i = 0; i < in_mat->nnz; i++)
	{
		out_mat->cscColPtr[in_mat->csrColIdx[i]]++;
	}	
	
	exclusive_scan(out_mat->cscColPtr, in_mat->n+1);
        
	ind_type *cscColIncr = (ind_type *)malloc(sizeof(ind_type) * (out_mat->n+1));
    memcpy (cscColIncr, out_mat->cscColPtr, sizeof(ind_type) * (out_mat->n+1));

        
	for(int row = 0; row < in_mat->m; row++)
	{

		for(int j = in_mat->csrRowPtr[row]; j < in_mat->csrRowPtr[row+1]; j++)
		{            
			int col = in_mat->csrColIdx[j];				
			out_mat->cscRowIdx[cscColIncr[col]] =  row;
			out_mat->cscVal[cscColIncr[col]] = in_mat->csrVal[j];
            cscColIncr[col]++;			            
		}
	}	
	free(cscColIncr);
    return EXIT_SUCCESS;	
	
}

int SILU::matrix_transpose(csc_matrix *in_mat, csr_matrix *out_mat)
{
	return EXIT_SUCCESS;
}

int SILU::exclusive_scan(ind_type *ptr, sz_type length)
{
	if(length == 0 || length == 1)
        return EXIT_SUCCESS;
    ind_type prev_val, new_val;

    prev_val = ptr[0];
    ptr[0] = 0;

    for(int i = 1; i < length; i++)
    {
    	new_val = ptr[i];
    	ptr[i] = prev_val + ptr[i-1];
    	prev_val = new_val;
    }

    return EXIT_SUCCESS;
}

int SILU::SerRef(val_type *x, val_type *x_ref, val_type *b, csc_matrix *inmat, int rhs, int ref_mode, int substitution)
{	
	
    SetRefVars(x_ref, b, inmat, rhs, ref_mode);
}

int SILU::SetRefVars(val_type *x_ref, val_type *b, csc_matrix *inmat, int rhs, int ref_mode)
{
	
   	if(ref_mode == RANDOM_REF)
	{
		for ( int i = 0; i < inmat->n; i++)
        	for (int j = 0; j < rhs; j++)
            {
            	x_ref[i * rhs + j] = rand() % 10 + 1; 
            }	

	        for (int i = 0; i < inmat->m * rhs; i++)
    	    	b[i] = 0;

		    for (int i = 0; i < inmat->n; i++)
		    {
		        for (int j = inmat->cscColPtr[i]; j < inmat->cscColPtr[i+1]; j++)
		        {
		            int rowid = inmat->cscRowIdx[j]; 
		            for (int k = 0; k < rhs; k++)
		            {
		                b[rowid * rhs + k] += inmat->cscVal[j] * x_ref[i * rhs + k];
		            }
		        }
		    }    	
	}    

	return EXIT_SUCCESS;
}

int SILU::ValidateResult(val_type *x, val_type *x_ref, sz_type n, int rhs)
{
	double accuracy = 1e-4;
    double ref = 0.0;
    double res = 0.0;

    for (int i = 0; i < n * rhs; i++)
    {
        ref += abs(x_ref[i]);
        if(policy.splitting_policy == DAG_SPLIT)
            res += abs(x[i] - x_ref[levels.levelItem[i]]);
        else
            res += abs(x[i] - x_ref[i]);    
    }
    res = ref == 0 ? res : res / ref;

    if (abs(res) < accuracy)
        printf("SILU test passed! |x-xref|/|xref| = %8.4e\n", res);
    else
        printf("SILU test NOT passed! |x-xref|/|xref| = %8.4e\n", res);

    return EXIT_SUCCESS;
}

int SILU::Test(csr_matrix *in_mat, val_type *x, val_type *b)
{
    struct matrix_descr testdescrA;
    sparse_matrix_t testcsrA;

    testdescrA.type = SPARSE_MATRIX_TYPE_GENERAL;
    testdescrA.mode = SPARSE_FILL_MODE_UPPER; 
    testdescrA.diag = SPARSE_DIAG_NON_UNIT;    
    sp_status = mkl_sparse_d_create_csr(&testcsrA, SPARSE_INDEX_BASE_ZERO, in_mat->m,
                                            in_mat->n, in_mat->csrRowPtr, in_mat->csrRowPtr + 1,
                                            in_mat->csrColIdx, in_mat->csrVal);
    if(sp_status == SPARSE_STATUS_SUCCESS)
    {
        printf("MKL create matrix was success\n");
    }
    sp_status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, testcsrA, testdescrA, x, 0.0, b);
    if(sp_status == SPARSE_STATUS_SUCCESS)
    {
        printf("MKL matrix SpMV was success\n");
    }

    return EXIT_SUCCESS;
}

void SILU::PrintMethodsInfo()
{   
    printf("SPTRSV Method 1: ");
    if(policy.sptrsv_method1 == SYNC_FREE)
    {
        printf("Sync Free\n");
    }
    if(policy.sptrsv_method1 == CUSPARSE_V2)
    {
        printf("cuSPARSE v2\n");
    }
    if(policy.sptrsv_method1 == MKL)
    {
        printf("MKL\n");
    }
    if(policy.sptrsv_method1 == SLFC)
    {
        printf("SLFC\n");
    }
    if(policy.sptrsv_method1 == ELMR)
    {
        printf("ELMR\n");
    }
    if(policy.sptrsv_method1 == ELMC)
    {
        printf("ELMC\n");
    }

    printf("SpMV Platform: ");
    if(policy.spmv_platform == SPMV_CPU)
    {
        printf("CPU\n");
    }
    if(policy.spmv_platform == SPMV_GPU)
    {
        printf("GPU\n");
    }

    printf("SPTRSV Method 2: ");
    if(policy.sptrsv_method2 == SYNC_FREE)
    {
        printf("Sync Free\n");
    }
    if(policy.sptrsv_method2 == CUSPARSE_V2)
    {
        printf("cuSPARSE v2\n");
    }
    if(policy.sptrsv_method2 == MKL)
    {
        printf("MKL\n");
    }
    if(policy.sptrsv_method2 == SLFC)
    {
        printf("SLFC\n");
    }
    if(policy.sptrsv_method2 == ELMR)
    {
        printf("ELMR\n");
    }
    if(policy.sptrsv_method2 == ELMC)
    {
        printf("ELMC\n");
    }
}

float SILU::getSpTRSV1_Time()
{
    return time_sptrsv1;
}
float SILU::getSpTRSV2_Time()
{
    return time_sptrsv2;
}
float SILU::getSpMV_Time()
{
    return time_spmv;
}
float SILU::getd2h_Time()
{
    return time_d2h;
}
float SILU::getd2h_sptrsv1_Time()
{
    return time_d2h_sptrsv1;
}
float SILU::getd2h_spmv_Time()
{
    return time_d2h_spmv;
}
float SILU::geth2d_Time()
{
    return time_h2d;
}
float SILU::geth2d_sptrsv1_Time()
{
    return time_h2d_sptrsv1;
}
float SILU::geth2d_spmv_Time()
{
    return time_h2d_spmv;
}
float SILU::getd2h_final_Time()
{
    return time_d2h_final;
}
float SILU::getlevelAna_Time()
{
    return time_levels;
}
float SILU::getDagAna_Time()
{
    return time_dag_analysis;
}
float SILU::getSplit_Time()
{
    return time_split;
}
float SILU::getDataStruct_Time()
{
    return time_data_structures;
}
float SILU::getAlgoAna_Time()
{
    return time_algo_analysis;
}

int SILU::getSpTRSV1_platform()
{
    return policy.sptrsv1_platform;
}
int SILU::getSpMV_platform()
{
    return policy.spmv_platform;
}
int SILU::getSpTRSV2_platform()
{
    return policy.sptrsv2_platform;
}

int SILU::getSpTRSV1_method()
{
    return policy.sptrsv_method1;
}

int SILU::getSpTRSV2_method()
{
    return policy.sptrsv_method2;
}

int SILU::getSplit_policy()
{
    return policy.splitting_policy;
}

int SILU::getSplit_level()
{
    return policy.split_at_level;
}
int SILU::getSplit_row()
{
    return policy.split_at_row;
}

int SILU::getSplit_cumrow()
{
    return levels.cum_Rows[policy.split_at_level];
}

int SILU::getLevels()
{
    return levels.nlevel;
}


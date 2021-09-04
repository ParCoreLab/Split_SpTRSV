// SILU.h
// Header file for Split Execution Framework Class
// Created: 03-12-2019
// Author: Najeeb Ahmad

#ifndef __HEADER_SILU__
#define __HEADER_SILU__

#include "matrix.h"
#include "levels.h"
#include "policy.h"
#include "common.h"
#include "device_ds.h"

#include <cusparse_v2.h>
//#include <mkl.h>
#include <mkl_spblas.h>
#include <string>
#include <time.h>
#include <chrono>

using namespace std::chrono;

class SILU
{

protected:
	level_info levels;
	level_info levelsTop;
	level_info levelsBot;

	//hyb_csc_csr_csc_matrix hyb_crc;        // Hybrid column (CSC), row(CSR), column (CSR)
	csc_matrix TopTri_csc;
  	csc_matrix Mvp_csc;
  	csc_matrix BotTri_csc;
  	csr_matrix TopTri_csr;
  	csr_matrix Mvp_csr;
  	csr_matrix BotTri_csr;

  	// Park's
  	csr_matrix Top_csr;
  	csr_matrix Bot_csr;
	
	int GPU_device_id;
	device_data_csc device_csc_matrix_top;
	device_data_csr device_csr_matrix_top;
	device_data_csc device_csc_matrix_bottom;
	device_data_csr device_csr_matrix_bottom;
	device_data_csc device_csc_matrix_mvp;
	device_data_csr device_csr_matrix_mvp;
    // SyncFree
	int *d_graphInDegree_1;
	int *d_graphInDegree_2;
	val_type *d_left_sum_1;
	// Solution
	val_type *d_x;
	val_type *d_x_1;
	// RHS
	val_type *d_b;
	val_type *d_b_1;
	char *ready_1;
	char *ready_2;
    // CUSPARSE
	cusparseHandle_t cusparseHandle;
	cusparseStatus_t cusparseStatus;
	cusparseMatDescr_t descr;
	cusparseMatDescr_t descrL;
	cusparseMatDescr_t descrL1;
	double doubleone;
	double minusdoubleone;
	double doublezero;
	
	void *pBuffer;
	int pBufferSizeInBytesL;
	int pBufferSize;
	csrsv2Info_t infoA;

	void *pBuffer1;
	int pBufferSizeInBytesL1;
	int pBufferSize1;
	csrsv2Info_t infoA1;

	// MKL
	sparse_status_t sp_status;
	sparse_matrix_t csrA, csrLU, csrLU1;
	struct matrix_descr MKLdescrA, MKLdescrL, MKLdescrL1;
	MKL_INT *ipar;
	double  *dpar;
	MKL_INT ivar, ierr;
	int m;
	int n;
	int n_iter;
	int n_iter_LU;
        int n_iter_LU1;


        //HTS
    //typedef HTS<int, int, double> ihts; 
    //typedef typename ihts::Real Real;
    //typename ihts::CrsMatrix* T;
    //std::vector<int> ir;
    //std::vector<int> jc, p, q;
    //std::vector<double> v;
    //std::vector<double> r;


	// Timing
	struct timeval t1, t2;
	float time_sptrsv1;
	float time_sptrsv2;
	float time_spmv;
	float time_d2h;
	float time_d2h_sptrsv1;
	float time_d2h_spmv;
	float time_d2h_final;
	float time_h2d;
	float time_h2d_sptrsv1;
	float time_h2d_spmv;
	
	
	float time_levels;
	float time_dag_analysis;
	float time_split;
	float time_data_structures;
	float time_algo_analysis;
	float time_factorize;
	cudaEvent_t start, stop;
	high_resolution_clock::time_point hrt1;
	high_resolution_clock::time_point hrt2;

	int timing;
	//val_type *x_ref;

	int SyncFreeAnalyzer(device_data_csc *device_csc_matrix, int *inDegree);
	int SyncFreeExecutor(device_data_csc *device_csc_matrix, val_type *d_x, val_type *d_b, 
		                 int *inDegree, int direction);
	int SLFCExecutor(device_data_csc *device_csc_matrix, val_type *d_x, 
                       val_type *diag, ind_type *jlev, int *inDegree);
	int ELMRExecutor(device_data_csr *device_csc_matrix, val_type *d_x, 
                       val_type *d_b, char *ready);
	int ELMCExecutor(device_data_csc *device_csc_matrix, val_type *d_x, 
                       val_type *d_b, ind_type *count);
	int HTSExecutor(device_data_csr *device_csr_matrix);

	// Smart ILU initialization routine
	// @return Initialization status

	int findlevels(csr_matrix *mat_csr, csc_matrix *mat_csc);
	int findlevels_csr(csr_matrix *tri_mat);
	int findlevels(csr_matrix *mat_csr, csc_matrix *mat_csc, level_info *li);
	int findlevels_csr(csr_matrix *tri_mat, level_info *li);
	int matrix_transpose(csc_matrix *input, csr_matrix *output);
	int exclusive_scan(ind_type *ptr, sz_type length);
	int analyse_dag(csr_matrix *tri_mat);
	int analyse_dag_revised(csr_matrix *tri_mat);
	int split_dag(csr_matrix *tri_mat);
	int split_dag2(csr_matrix *tri_mat);
	int split_matrix(csr_matrix *mat, csr_matrix *tri_mat);
	int setup_execution_devices(void);
	int setup_device_datastructures(val_type *b);
	int algorithm_analysis();
	int cut_square_matrix(csr_matrix *in, csr_matrix *out, int start_row,int end_row);
	int end_profile_CPU(float *elapsed_time);
	int end_profile_GPU(float *elapsed_time);
	int start_profile_CPU();
	int start_profile_GPU();
	int free_memory();
	int free_memory_1();
	void init_ptrs();

public:
	execution_policy policy;
	int ELMC_only;
    int MKL_only;
	SILU();
	~SILU();
	int Init();

	// Smart ILU Factorization routine
	// @return Factorization status
	// @param  A  input csr matrix
	// @param  L  output lower triangular hybrid matrix
	// @param  U  output upper triangular hybrid matrix
	// @param  method ILU0 factorization method
	//         Options: dummy: Upper and lower parts of A used as ILU0 factors
	//                  MKL: MKL ILU0 factorization is used
	//                  cuSPARSE: cuSPARSE ILU0 factorization is used
	int SetRefVars(val_type *x_ref, val_type *b, csc_matrix *inmat, int rhs, int ref_mode);
	int SerRef(val_type *x, val_type *x_ref, val_type *b, csc_matrix *inmat, 
		             int rhs, int ref_mode, int substitution);
	int ValidateResult(val_type *x, val_type *x_ref, sz_type n, int rhs);
	int Factorize(csr_matrix *A, csr_matrix *L, csr_matrix *U, std::string method);

	// Smart ILU Analysis routine
	// @return Analysis status
	int Analyze(csr_matrix *mat, csr_matrix *tri_mat, val_type *b);

	int SetExecutionPolicy(int policy);
	int matrix_transpose(csr_matrix *input, csc_matrix *output);

	void PrintLevelRows(void);
	//void PrintCsrMatrix(csr_matrix *mat);
	void PrintCsrMatrix(csr_matrix *mat, char rows='a', int norows=20);
	void PrintCscMatrix(csc_matrix *mat);
	void PrintLevels();
	void PrintCumNnzs();
	void PrintLevels(level_info *levels);
	void PrintMethodsInfo();
	int SpTRSV_Time_Breakup();
	int Print_SubMatrix_stats();
	float getSpTRSV1_Time();
	float getSpTRSV2_Time();
	float getSpMV_Time();
	float getd2h_Time();
	float geth2d_Time();
	float getd2h_final_Time();
	float getlevelAna_Time();
	float getDagAna_Time();
	float getSplit_Time();
	float getDataStruct_Time();
	float getAlgoAna_Time();
	float getd2h_sptrsv1_Time();
    float getd2h_spmv_Time();
    float geth2d_sptrsv1_Time();
    float geth2d_spmv_Time();
    int getSpTRSV1_platform();
    int getSpMV_platform();
    int getSpTRSV2_platform();
    int getSpTRSV1_method();
    int getSpTRSV2_method();
    int getSplit_policy();
    int getSplit_level();
    int getSplit_row();
    int getSplit_cumrow();
    int getLevels();
    int SaveLevelFeatures(int mat_id);
    int SILU::SavePolicy(std::string filename, int mat_id);
    int SaveMatrixFeatures(csr_matrix *tri_mat, std::string filename);
    int findlevels_incremental(csr_matrix *tri_mat, std::string filename);
    int PreProcess(csr_matrix *tri_mat);
    int GPUGPUSplit(csr_matrix *tri_mat);
    int GPUGPUSplit_old(csr_matrix *tri_mat);
    int CPUGPUSplit(csr_matrix *tri_mat);
    int CPUCPUSplit(csr_matrix *tri_mat, int hts_mkl);
    int DecidePlatform(csr_matrix *tri_mat);
    int GPUSplitORUnified(csr_matrix *tri_mat);
    int FindSplitPoint(csr_matrix *tri_mat);
    int MKL_HTS_LogReg();
    int HTS_ELMC_LogRes();
	// Smart ILU triangular solve routine
	// @return triangular solve status
	int trsv(val_type *rhs, val_type *x, int direction);
	int Test(csr_matrix *in_mat, val_type *x, val_type *b);	
};

#endif


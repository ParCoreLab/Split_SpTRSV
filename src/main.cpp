// Split-Execution Framework driver code
// Created: 03-12-2019
// Author: Najeeb Ahmad


#include "common.h"
#include "matrix.h"
#include "util.h"
#include "SILU.h"
#include "quicksort.h"
#include "policy.h"
#include "mmio.h"
#include <cuda_runtime.h>

using namespace std;
using namespace std::chrono;

int main(int argc, char **argv)
{
    if(argc < 3)
    {
    	std::cout << "Error: Missing arguments\n";
    	std::cout << "Usage: " << argv[0] << " matrix_id mode\n";
    	return EXIT_FAILURE;
    }

    csr_matrix matrix;
    string matrix_name;
    string mode_str(argv[2]);
    int mode = stoi(mode_str);
    ofstream outfile("results.txt", ios::app);
    ofstream outfile1;
    if(mode==1) // hybrid
    {
        outfile1.open("results_for_revision_hybrid.txt", ios::app);
    }
    if(mode==2)
    {    // ELMC}
        outfile1.open("results_for_revision_ELMC.txt", ios::app);
    }
    if(mode==3)
    {
        outfile1.open("results_for_revision_MKL.txt", ios::app);
    }     
    
    
    if(outfile.is_open())
    {
    }
    else
    {
        cout << "Failed to open results.txt file!" << endl;
        return EXIT_FAILURE;
    }
    if(outfile1.is_open())
    {
    }
    else
    {
        cout << "Failed to open results_for_revision_XXX file!" << endl;
        return EXIT_FAILURE;
    }

    if (argc == 4)     // If matrix is to be read from .mtx file
    {
        printf("Reading .mtx file\n");
        int retCode = 0;
	    matrix_name = argv[4];
        retCode = mm_read_unsymmetric_sparse(argv[3], &matrix.m, &matrix.n, &matrix.nnz,
                &matrix.csrVal, &matrix.csrRowPtr, &matrix.csrColIdx);
        coo2csr_in(matrix.m, matrix.nnz, matrix.csrVal, matrix.csrRowPtr, matrix.csrColIdx);
        printf("Done reading file\n");
        printf("Matrix Rows: %d\n", matrix.m);
        printf("Matrix Cols: %d\n", matrix.n);
        
        if(retCode == -1)
        {
            cout << "Error reading input .mtx file\n";
            return EXIT_FAILURE;
        }
    }
    else
    {
        printf("Downloading matrix from Suitsparse\n");
        string matrix_id(argv[1]), output_folder(argv[2]);
        int mat_id = stoi(matrix_id);

        uf_matrix_reader reader;
        matrix_name = reader.get_matrix_name(mat_id);
        if (does_file_exist(path_join(output_folder, matrix_name))) 
        {
            std::cout << "Already calculated for " << matrix_name << "\n";
            return EXIT_SUCCESS;
        }
        matrix = reader.get_matrix_csr_libufget(mat_id);

    }
       
    float sptrsv1_time_new = 0.0;
    float sptrsv2_time_new = 0.0;
    float spmv_time_new = 0.0;
    float d2h_time = 0.0;
    float d2h_sptrsv1_time_new = 0.0;
    float d2h_spmv_time_new = 0.0;
    float h2d_time = 0.0;
    float h2d_sptrsv1_time_new = 0.0;
    float h2d_spmv_time_new = 0.0;
    float d2h_time_final_new = 0.0;
    float levels_time_new = 0.0;
    float dag_analysis_time_new = 0.0;
    float split_time_new = 0.0;
    float data_structures_time_new = 0.0;
    float algo_analysis_time_new = 0.0;
    float total_sptrsv_time_new = 0.0;
    float total_analysis = 0.0;
    float total_comm = 0.0;
    int SpTRSV1_method = 0;
    int SpTRSV2_method = 0;
    int SpMV_platform = 0;
    int split_level = 0;
    int split_row = 0;
    int split_policy = 0;
    int ITERATIONS = 10;
    int policy = 0;
    int split_level_new = 0;
    
#ifdef AUTO_TESTING
    // For manual testing
    int algo1 = stoi(argv[3]);
    int algo2 = stoi(argv[4]);
    int spmv_plat = stoi(argv[5]);
    int splitting_level = stoi(argv[6]);
    printf("algo1: %d\n", algo1);
    printf("algo2: %d\n", algo2);
    printf("spmv_plat: %d\n", spmv_plat);
    printf("splitting_level: %d\n", splitting_level);
#endif

    printf("Starting iterations\n");
    for(int i = 0; i < ITERATIONS; i++)
    {
        csr_matrix L, U;
        csc_matrix L_csc;
        val_type *y;
        SILU silu_test;

#ifdef AUTO_TESTING
    silu_test.policy.sptrsv_method1 = algo1;
    silu_test.policy.sptrsv_method2 = algo2;
    silu_test.policy.spmv_method = SPMV_CSR;
    silu_test.policy.spmv_platform = spmv_plat;
    silu_test.policy.split_at_level = splitting_level;
#endif
        switch(mode)
        {
            case 1:
            silu_test.ELMC_only = 0;
            silu_test.MKL_only = 0;
            break;
            case 2:
            silu_test.ELMC_only = 1;
            break;
            case 3:
            silu_test.MKL_only = 1;
            break;
        }
        // Dummy ILU factorization. In dummy factorization, we are only
	// intersted in the sprsity pattern, not values        
        silu_test.Factorize(&matrix, &L, &U, "dummy");
        int rhs = 1;
        val_type *b;
	cudaError_t status = cudaMallocHost((void **)&b, sizeof(val_type)*matrix.n*rhs);
	val_type *x = (val_type *)malloc(sizeof(val_type) * matrix.n * rhs);
        val_type *x_ref = (val_type *)malloc(sizeof(val_type) * matrix.n * rhs);
        silu_test.matrix_transpose(&L, &L_csc);
        silu_test.SerRef(x, x_ref, b, &L_csc, rhs, RANDOM_REF, SUBSTITUTION_FORWARD);
        
	// Analysis phase
        if(silu_test.Analyze(&matrix, &L, b) == EXIT_SUCCESS)
        {        
            printf("Levels: %d\n", silu_test.getLevels());
            silu_test.PrintMethodsInfo();
            y = (val_type *)malloc(sizeof(val_type) * matrix.m);
	    memset(y, 0, (matrix.m) *sizeof(ind_type));
            // Perform sparse triangular solve (split/unified) as per policy
	    // decided in Analysis phase
	    silu_test.trsv(b, y, SUBSTITUTION_FORWARD);        
            silu_test.ValidateResult(y, x_ref,matrix.n, 1);
            printf("-------------------------------\n");                
        }
        else
        {
            printf("Errors encountered. Exiting program.\n");
            return EXIT_FAILURE;
        }
        
	cudaFree(b);
        free(x);
        free(x_ref);
        free(y);
	if(i == ITERATIONS - 1)
        {
            SpTRSV1_method = silu_test.getSpTRSV1_method();
            SpTRSV2_method = silu_test.getSpMV_platform();
            SpMV_platform = silu_test.getSpTRSV2_method();
            split_policy = silu_test.getSplit_policy(); 
            if(silu_test.getSplit_policy() == DAG_SPLIT)
            {
                split_level = silu_test.getSplit_level();
                split_row = silu_test.getSplit_cumrow();
            }
            else
            {
                split_row = silu_test.getSplit_row();
            }
            
	    string matrix_id(argv[1]);
	    int mat_id = stoi(matrix_id);
	    silu_test.SaveLevelFeatures(mat_id);
	    silu_test.SavePolicy("exec_policy", mat_id);
	    silu_test.Print_SubMatrix_stats();
	    
        }    
	// For timing
	switch(silu_test.policy.description)
	  {
	  case 0:   // CPU only
	    sptrsv1_time_new = silu_test.getSpTRSV1_Time();
	    total_sptrsv_time_new += sptrsv1_time_new;
	    total_comm = 0.0;
	    break;
	  case 1:   // GPU only
	    sptrsv1_time_new = silu_test.getSpTRSV1_Time();
	    total_sptrsv_time_new += sptrsv1_time_new;
	    total_comm = 0.0;
	    break;
	  case 2:  // CPU-GPU
	    sptrsv1_time_new = silu_test.getSpTRSV1_Time();
	    sptrsv2_time_new = silu_test.getSpTRSV2_Time();
	    spmv_time_new = silu_test.getSpMV_Time();
	    h2d_sptrsv1_time_new = silu_test.geth2d_sptrsv1_Time();
	    total_sptrsv_time_new += sptrsv1_time_new + spmv_time_new + sptrsv2_time_new + h2d_sptrsv1_time_new;
	    total_comm += h2d_sptrsv1_time_new;                
	    break;
	  case 3: // GPU-CPU
	    sptrsv1_time_new = silu_test.getSpTRSV1_Time();
	    sptrsv2_time_new = silu_test.getSpTRSV2_Time();
	    spmv_time_new = silu_test.getSpMV_Time();
	    d2h_spmv_time_new = silu_test.getd2h_spmv_Time();
	    total_sptrsv_time_new += sptrsv1_time_new + spmv_time_new + sptrsv2_time_new + d2h_spmv_time_new;
	    total_comm += d2h_spmv_time_new;
	    break;
	  case 4: // GPU-GPU
	    sptrsv1_time_new = silu_test.getSpTRSV1_Time();
	    sptrsv2_time_new = silu_test.getSpTRSV2_Time();
	    spmv_time_new = silu_test.getSpMV_Time();
	    total_sptrsv_time_new += sptrsv1_time_new + spmv_time_new + sptrsv2_time_new;
	    total_comm = 0.0;
	    break;
	  case 5: // MKL
	    sptrsv1_time_new = silu_test.getSpTRSV1_Time();
	    total_sptrsv_time_new += sptrsv1_time_new;
	    total_comm = 0.0;
	    break;
	  case 6: // HTS
	    sptrsv1_time_new = silu_test.getSpTRSV1_Time();
	    total_sptrsv_time_new += sptrsv1_time_new;
	    total_comm = 0.0;
	    break;	 
        }

        d2h_time_final_new += silu_test.getd2h_final_Time();
        levels_time_new = silu_test.getlevelAna_Time();
        dag_analysis_time_new = silu_test.getDagAna_Time();
        split_time_new = silu_test.getSplit_Time();
        data_structures_time_new = silu_test.getDataStruct_Time();
        algo_analysis_time_new = silu_test.getAlgoAna_Time();
        total_analysis += levels_time_new + dag_analysis_time_new + split_time_new + data_structures_time_new + algo_analysis_time_new;
        policy = silu_test.policy.description;
        split_level_new = silu_test.policy.split_at_level;
    }
    // Times
    printf("SPTRSV: %f, COMM: %f, D2H_FINAL:%f, ANALYSIS: %f\n", 
        total_sptrsv_time_new/ITERATIONS,total_comm/ITERATIONS, 
        d2h_time_final_new/ITERATIONS, total_analysis/ITERATIONS);
    //total_sptrsv_time /= ITERATIONS;    
    outfile1 << argv[1] << "," << total_sptrsv_time_new/ITERATIONS << ","
             << total_comm/ITERATIONS <<  "," << d2h_time_final_new/ITERATIONS << "," 
             <<  total_analysis/ITERATIONS << "," << policy << "," 
             << split_level_new << endl; 
    outfile1.close();
    printf("--------------------------------------------\n");
    printf("Results stored in file\n");
    printf("All times in milli-seconds\n");
    outfile << "SpTRSV1 method: " << SpTRSV1_method << endl;
    outfile << "SpMV platform: " << SpTRSV2_method << endl;
    outfile << "SpTRSV2 method: " << SpMV_platform << endl;
    outfile << "Split policy: " << split_policy << endl;
    outfile << "Split level: " << split_level << endl;
    outfile << "Split row: " << split_row << endl;    
    outfile.close();
    
    return EXIT_SUCCESS;
}

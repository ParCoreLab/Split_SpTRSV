// matrix.h
#pragma once

#include "common.h"
#include <libufget.h>

struct csr_matrix 
{
  ind_type * csrRowPtr;
  ind_type * csrColIdx;
  val_type * csrVal;
  sz_type m;
  sz_type n;
  sz_type nnz;
  std::string name;

  csr_matrix();
  ~csr_matrix();
  csr_matrix(ind_type * csrRowPtr,
             ind_type * csrColIdx,
             val_type * csrVal,
             sz_type m,
             sz_type n,
             sz_type nnz,
             const std::string &name);
};

struct csc_matrix
{
  ind_type * cscRowIdx;
  ind_type * cscColPtr;
  val_type * cscVal;
  sz_type m;
  sz_type n;
  sz_type nnz;
  std::string name;

  csc_matrix();
  csc_matrix(ind_type * cscRowIdx,
             ind_type * cscColPtr,
             val_type * cscVal,
             sz_type m,
             sz_type n,
             sz_type nnz,
             const std::string &name);
};

struct hyb_csr_matrix
{
  csr_matrix TopTri;
  csr_matrix Mvp;
  csr_matrix BotTri;
  /*
   ind_type * csrRowPtr_TopTri;
   ind_type * csrColIdx_TopTri;
   val_type * csrVal_TopTri;
   sz_type m_TopTri;
   sz_type n_TopTri;
   sz_type nnz_TopTri;

   ind_type * csrRowPtr_Mvp;
   ind_type * csrColIdx_Mvp;
   val_type * csrVal_Mvp;
   sz_type m_Mvp;
   sz_type n_Mvp;
   sz_type nnz_Mvp;

   ind_type * csrRowPtr_BotTri;
   ind_type * csrColIdx_BotTri;
   val_type * csrVal_BotTri;
   sz_type m_BotTri;
   sz_type n_BotTri;
   sz_type nnz_BotTri;
   */

   //hyb_csr_matrix(){}
   /*hyb_csr_matrix(ind_type * const csrRowPtr_TopTri,
             ind_type * const csrColIdx_TopTri,
             val_type * const csrVal_TopTri,
             const sz_type m_TopTri,
             const sz_type n_TopTri,
             const sz_type nnz_TopTri,
             ind_type * const csrRowPtr_Mvp,
             ind_type * const csrColIdx_Mvp,
             val_type * const csrVal_Mvp,
             const sz_type m_Mvp,
             const sz_type n_Mvp,
             const sz_type nnz_Mvp,
             ind_type * const csrRowPtr_BotTri,
             ind_type * const csrColIdx_BotTri,
             val_type * const csrVal_BotTri,
             const sz_type m_BotTri,
             const sz_type n_BotTri,
             const sz_type nnz_BotTri);*/
};

struct hyb_csc_csr_csc_matrix
{
  csc_matrix TopTri;
  csr_matrix Mvp;
  csc_matrix BotTri;
};

class uf_matrix_reader 
{
 private:
  uf_collection_t* collection;
 public:
  uf_matrix_reader();
  ~uf_matrix_reader();

  csr_matrix get_matrix_csr_libufget(int id);
  std::string get_matrix_name(int id);
};

struct lu_matrices 
{
  csr_matrix l;
  csr_matrix u;

  lu_matrices(const csr_matrix &l, const csr_matrix &u);
};

lu_matrices get_ilu0_decomposition(const csr_matrix &matrix);


// matrix.cpp
#include "matrix.h"
#include "common.h"
csr_matrix::csr_matrix(ind_type * const csrRowPtr,
             		   ind_type * const csrColIdx,
             		   val_type * const csrVal,
             		   const sz_type m,
             	     const sz_type n,
             		   const sz_type nnz,
             		   const std::string &name)
    : csrRowPtr(csrRowPtr), csrColIdx(csrColIdx), csrVal(csrVal), m(m), n(n), nnz(nnz), name(name) {}

csc_matrix::csc_matrix(ind_type * const cscRowIdx,
                   ind_type * const cscColPtr,
                   val_type * const cscVal,
                   const sz_type m,
                   const sz_type n,
                   const sz_type nnz,
                   const std::string &name)
    : cscRowIdx(cscRowIdx), cscColPtr(cscColPtr), cscVal(cscVal), m(m), n(n), nnz(nnz), name(name) {}

uf_matrix_reader::uf_matrix_reader()
    : collection(uf_collection_init1(UF_COLLECTION_VERBOSE)) {

  std::cout << "Database contains " << uf_collection_num_matrices(collection) << " matrices.\n";

}

csr_matrix::csr_matrix()
{

}

csc_matrix::csc_matrix()
{
  
}

csr_matrix::~csr_matrix()
{
  //free(csrRowPtr);
  //free()
}

csr_matrix uf_matrix_reader::get_matrix_csr_libufget(int id) {
  printf("ID: %d\n", id);
  uf_matrix_t *mat = uf_collection_get_by_id(collection, id);
  uf_matrix_get(mat);

  uf_field_t field;
  sz_type m, n, nnz;
  ind_type *csrRowPtr, *csrColIdx;
  val_type *csrVal;

  
  int retCode = uf_matrix_coord_int32(&field, &m, &n, &nnz, &csrRowPtr, &csrColIdx, (void **) &csrVal, mat);
  
  if(retCode != 0)
  {
    printf("Error reading the matrix archive file!\n");
    printf("Function return code: %d\n", retCode);
    exit(1);
  }

  uf_matrix_coord_to_csr_int32(field, m, n, nnz, &csrRowPtr, &csrColIdx, (void **) &csrVal);

  printf("matrix name = %s\n", mat->name);
  printf("Rows: %10d\t Cols: %10d\t nnz: %12d\n", (int) m, (int) n, (int) nnz);

  //csrRowPtr = (ind_type *)malloc((m+1) * sizeof(ind_type));
  //csrColIdx = (ind_type *)malloc(nnz * sizeof(ind_type));
  //csrVal    = (val_type *)malloc(nnz * sizeof(val_type));

  uf_matrix_free(mat);

  return csr_matrix{
      csrRowPtr,
      csrColIdx,
      csrVal,
      m,
      n,
      nnz,
      mat->name
  };
}

uf_matrix_reader::~uf_matrix_reader() {
  uf_collection_finalize(collection);
}

std::string uf_matrix_reader::get_matrix_name(int id) {
  uf_matrix_t *mat = uf_collection_get_by_id(collection, id);
  std::string name(mat->name);
  uf_matrix_free(mat);
  return name;
}
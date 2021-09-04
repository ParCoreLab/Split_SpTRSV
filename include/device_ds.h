#ifndef __HEADER_DEVICEDS__
#define __HEADER_DEVICEDS__

#include "common.h"

struct device_data_csc
{
	ind_type *d_cscColPtr;
    ind_type *d_cscRowIdx;
    val_type *d_cscVal;
    sz_type m;
    sz_type n;
    sz_type nnz;
    //val_type *d_b;
    val_type *d_x;
    val_type *d_left_sum;
    val_type *diag;
    ind_type *jlev;
};

struct device_data_csr
{
	ind_type *d_csrRowPtr;
    ind_type *d_csrColIdx;
    val_type *d_csrVal;
    sz_type m;
    sz_type n;
    sz_type nnz;
    //val_type *d_b;
    val_type *d_x;
    val_type *d_left_sum;
    val_type *diag;
    ind_type *jlev;
};

#endif
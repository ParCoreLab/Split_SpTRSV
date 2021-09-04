/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright (C) Martin Koehler, 2015
 */

#ifndef LIBUFGET_H
#define LIBUFGET_H

#ifdef __cplusplus
extern "C" {
#endif
    #include <stdint.h>

    #define LIBUFGET_MAJOR 1
    #define LIBUFGET_MINOR 0
    #define LIBUFGET_PATCH 2

    #define LIBUFGET_VERSION (((LIBUFGET_MAJOR*100)+LIBUFGET_MINOR)*100+LIBUFGET_PATCH)


    /* Forward definition for internal usage  */
    struct _uf_collection_t;
    struct _uf_query_t;

    /**
     * \brief Opaque object referencing the collection.
     *
     * The \b uf_collection_t object contains the information about the currently used collection.
     * It contents are not vissible to the API because of internal changes in later versions.
     * Therefore, the object is only accessed via a pointer. Never use \b uf_collection_t as a
     * non-pointer variable.
     *
     * \ingroup collection
     */
    typedef struct _uf_collection_t uf_collection_t;

    /**
     * \brief Opaque object referencing a query.
     *
     * The \b uf_query_t object is referencing a query to the collection. The API only requires
     * a forward reference to the internal data structure. Therefore, the user should only use
     * the \b uf_query_t object as a pointer.
     *
     * \ingroup query
     */
    typedef struct _uf_query_t uf_query_t;

    /**
     * \brief Representation of a matrix entry in the collection.
     *
     * The uf_matrix_t collection represents a matrix in the collection with all the information from
     * the data base.
     *
     * \ingroup matrix
     * */
    typedef struct _uf_matrix_t {
        int id;                        /**<\brief ID of the Matrix from the collection. This is the same ID as on
                                          the UF collection website. */
        char group_name[1024];         /**<\brief The name of the group where the matrix belongs to. */
        char name[1024];               /**<\brief The name of the matrix. The full name of a matrix should always
                                          be given as \b group_name/name. */
        int64_t rows;                  /**<\brief The number of rows in the matrix.  */
        int64_t cols;                  /**<\brief The number of columns in the matrix.  */
        int64_t nnz;                   /**<\brief The number of numerically nonzero entries in the matrix. */
        int64_t nzero;                 /**<\brief The number of explicit entries present in the matrix that are provided
                                          by the matrix author but which are numerically zero.  */
        double pattern_symmetry;       /**<\brief Let S be the binary pattern of the explicit nonzeros in the matrix.
                                          Let pmatched be the number of matched offdiagonal entries, where both S(i,j)
                                          and S(j,i) are one, with i not equal to j. Let nzoffdiag be the number of
                                          offdiagonal entries in S. Then pattern_symmetry is the ratio of pmatched/nzoffdiag.
                                          Note that if S(i,j) and S(j,i) are both one, then this pair of entries is counted
                                          twice in both pmatched and nzoffdiag. If the matrix is rectangular, this statistic
                                          is zero. If there are no offdiagonal entries, the statistic is 1.   */
        double numerical_symmetry;     /**<\brief Let xmatched be the number of matched offdiagonal entries, where A(i,j)
                                          is equal to the complex conjugate of A(j,i) and where i and j are not equal.
                                          Then numerical_symmetry is the ration xmatched / nzoffdiag (or 1 if nzoffdiag
                                          is zero). This statistic is zero for rectangular matrices. Note that this statistic
                                          measures how close a matrix is to being Hermitian (A=A' in MATLAB). For complex
                                          symmetric matrices (A=A.' in MATLAB), this ratio will be less than one (unless
                                          all offdiagonal entries are real). */
        int isBinary;                  /**<\brief  1 if the matrix is binary (all entries are 0 or 1), 0 otherwise. */
        int isReal;                    /**<\brief  1 if the matrix is real, 0 otherwise. */
        int isComplex;                 /**<\brief  1 if the matrix is complex, 0 otherwise. */
        int64_t nnzdiag;               /**<\brief The number of numerically nonzero entries on the diagonal
                                          (nnz (diag (Problem.A)) in MATLAB notation).*/
        int posdef;                    /**<\brief 1 if the matrix is known to be symmetric positive definite (or Hermitian
                                          positive definite for the complex case), 0 if it is known not to be, -1 if it is
                                          symmetric (or Hermitian) but with unknown positive-definiteness. If the statistic
                                          is unknown (-1) this may be revised in subsequent versions of the index. */
        int64_t amd_lnz;               /**<\brief Let C=S+S' where S=spones(A) is the binary pattern of A. Then amd_lnz
                                          is number of nonzeros in the Cholesky factorization of the matrix C(p,p)
                                          (assuming C is positive definite) where p=amd(C) is the AMD fill-reducing ordering.
                                          This statistic is -2 for rectangular matrices or for graph problems.
                                          It is -1 if it is not computed. This figure gives an estimate of the memory
                                          requirements for the direct solution of Ax=b if A is a square matrix. */
        int64_t amd_flops;             /**<\brief The floating-point operation count for computing the Cholesky
                                          factorization of C(p,p). See \ref amd_lnz. */
        int64_t amd_vnz;               /**<\brief The number of entries in a Householder-vector representation of
                                          the Q factor of R, after a COLAMD fill-reducing ordering. This is an upper
                                          bound on L for the LU factorization of A. */
        int64_t amd_rnz;               /**<\brief The number of entries in R for the QR factorization of A, after a COLAMD
                                          fill-reducing ordering. This is an upper bound on U for the LU factorization of A. */
        int64_t nblocks;               /**<\brief The number of blocks from the Dulmage-Mendelsohn Decomposition. */
        int64_t sprank;                /**<\brief The structural rank of the matrix, which is the number of maximual
                                          number of nonzero entries that can be permuted to the diagonal.
                                          This statistic is not computed for problems that represent graphs,
                                          since in those cases the diagonal of the matrix is often not relevant
                                          (self-edges are often ignored). See also \ref nblocks. */
        char RBtype[4];                /**<\brief The Rutherford Boeing type of the matrix. RBtype is a a 3-letter
                                          lower-case string.  */
        int cholcand;                  /**<\brief 1 if the matrix is symmetric (Hermitian if complex) and if all
                                          the diagonal entries are postive and real. Zero otherwise. If 1, the matrix
                                          is a candidate for a Cholesky factorization, */
        int64_t ncc;                   /**<\brief The number of of strongly-connected components in the graph of A
                                          (if A is square) or in the bipartite graph if A is rectangular. The diagonal
                                          is ignored. */
        int isND;                      /**<\brief  1 if the problem comes from a 2D/3D discretization, zero otherwise.
                                          This determination is not a property of the matrix, but a qualitative assesment
                                          of the kind of problem the matrix represents. */
        int isGraph;                   /**<\brief 1 if the problem is best considered as a graph rather than a system
                                          of equations, zero otherwise. This determination is not a property of the matrix,
                                          but a qualitative assesment of the kind of problem the matrix represents. */
        int64_t lowerbandwidth;        /**<\brief Bandwidth of the left lower triangle. */
        int64_t upperbandwidth;        /**<\brief Bandwidth of the upper right lower triangle. */
        int64_t rcm_lowerbandwidth;    /**<\brief Bandwidth of the left lower triangle after RCM reordering. */
        int64_t rcm_upperbandwidth;    /**<\brief Bandwidth of the upper right lower triangle after RCM reordering. */
        double xmin_real;              /**<\brief Undocumented information from the database. (Real Part)*/
        double xmin_imag;              /**<\brief Undocumented information from the database. (Imaginary Part)*/
        double xmax_real;              /**<\brief Undocumented information from the database. (Real Part)*/
        double xmax_imag;              /**<\brief Undocumented information from the database. (Imaginary Part)*/

        /* Internal data */
        uf_collection_t *col;           /**< \brief \b Internal.*/
        char *localpath;                /**< \brief \b Internal. */
    } uf_matrix_t;

    /**
     * \brief Enumerator to identify the value type of the matrix entries.
     *
     * The \b uf_field_t enumerator is used to identify if a matrix is real, binary, or complex. It is used by the
     * convert routines to select which typecast is applied to the values. For easier usage binary matrices will
     * be stored as realed value double precision arrays with only 0 and 1 entries.
     *
     * \ingroup matrix
     */
    typedef enum {
        UF_REAL = 1,     /**< \brief The entries of the matrix are real.  */
        UF_PATTERN = 2,  /**< \brief The entries of the matrix are only binary values (0 and 1). */
        UF_COMPLEX = 3   /**< \brief The entries of the matrix are complex values.  */
    } uf_field_t;

    /* Initialization and finalization  */

    /**
     * \brief Default flag for the initialization of the collection.
     *
     * The UF_COLLECTION_DEFAULT macro sets the default flags for the initialization of the library. This means
     * that the library does not perform verbose output, the database is only refreshed if it is older than one
     * month and all downloaded matrices will be cached.
     * \sa uf_collection_init
     *
     * \ingroup collection */
    #define UF_COLLECTION_DEFAULT       0x00

    /**
     * \brief Flag to force the refresh of the database.
     *
     * The UF_COLLECTION_FORCE_REFRESH flag forces the refresh of the database. The flag is orthogonal to the other
     * flags which means that they could be combined using the binary-or operation.
     * \sa uf_collection_init
     * \ingroup collection
     */
    #define UF_COLLECTION_FORCE_REFRESH 0x01

    /**
     * \brief Flag to turn on verbose outputs.
     *
     * The UF_COLLECTION_VERBOSE flag enables some output for the library, like the download progress and similar stuff.
     * The flag is orthogonal to the other flags which means that they could be combine using binary-or operation.
     * \sa uf_collection_init
     * \ingroup collection
     * */
    #define UF_COLLECTION_VERBOSE       0x02

    /**
     * \brief Flag to enable the cleanup of downloaded matrices.
     *
     * If the UF_COLLECTION_CLEANUP flag is set the library removes \b all downloaded matrices from the cache after the
     * collection is deinitialized.
     * The flag is orthogonal to the other flags which means that they could be combine using binary-or operation.
     * \sa uf_collection_init
     * \ingroup collection
     */
    #define UF_COLLECTION_CLEANUP       0x04

    /* Collection Access  */

    /**
     * \brief Default initialization routine for the library.
     * \return A pointer to a initialized \ref uf_collection_t object or NULL in case of an error.
     *
     * The uf_collection_init function is the default initialization routine for the matrix collection. It
     * specifies the \ref UF_COLLECTION_DEFAULT flag and reads the following environment variables to sets
     * the corresponding flag if their value is larger than zero:
     *
     * \li \b UF_COLLECTION_VERBOSE, sets the \ref UF_COLLECTION_VERBOSE flag,
     * \li \b UF_COLLECTION_CLEANUP, sets the \ref UF_COLLECTION_CLEANUP flag,
     * \li \b UF_COLLECTION_BASEURL, sets an alternative base url passed to \ref uf_collection_init2,
     * \li \b UF_COLLECTION_CACHE_DIR, sets an alternative cache directory passed to \ref uf_collection_init2.
     *
     * The default base url used by libufget is http://www.cise.ufl.edu/research/sparse/. The default cache
     * directory is \b $HOME/\b .ufget/.
     *
     * \sa uf_collection_init1
     * \sa uf_collection_init2
     * \sa uf_collection_finalize
     *
     * \ingroup collection
     */
    uf_collection_t * uf_collection_init();

    /**
     * \brief Extended  initialization routine for the library.
     * \param[in] flags  A set of the initialization flags combine with binary or.
     * \return A pointer to a initialized \ref uf_collection_t object or NULL in case of an error.
     *
     * The uf_collection_init1 function is the extended initialization routine for the matrix collection. The different
     * initialization flags like \ref UF_COLLECTION_VERBOSE, \ref UF_COLLECTION_CLEANUP, or  \ref UF_COLLECTION_FORCE_REFRESH
     * can be given. Furthermore, it reads the following environment variables:
     *
     * \li \b UF_COLLECTION_BASEURL, sets an alternative base url passed to \ref uf_collection_init2,
     * \li \b UF_COLLECTION_CACHE_DIR, sets an alternative cache directory passed to \ref uf_collection_init2.
     *
     * The default base url used by libufget is http://www.cise.ufl.edu/research/sparse/. The default cache
     * directory is \b $HOME/\b .ufget/.
     *
     * \sa uf_collection_init
     * \sa uf_collection_init2
     * \sa uf_collection_finalize
     *
     * \ingroup collection
     */
    uf_collection_t * uf_collection_init1(unsigned long flags);


    /**
     * \brief Fine-grained initialization routine for the library.
     * \param[in] baseurl  Base url of the collection to use.
     * \param[in] cache_dir Cache directory to cache the downloaded matrices to.
     * \param[in] flags  A set of the initialization flags combine with binary or.
     * \return A pointer to a initialized \ref uf_collection_t object or NULL in case of an error.
     *
     * The uf_collection_init2 function is the fine-grained initialization routine for the matrix collection. It allows
     * to specify the flags like \ref uf_collection_init1 but additionally gets the base url and the cache directory
     * as parameters. For this reason it does not read any environment variable nor it has default values for them.
     * \sa uf_collection_init1
     * \sa uf_collection_init
     * \sa uf_collection_finalize
     *
     * \ingroup collection
     */
    uf_collection_t * uf_collection_init2(const char *baseurl, const char *cache_dir,  unsigned long flags);

    /**
     * \brief Close a given collection.
     * \param[in] collection Collection to close.
     *
     * The uf_collection_finalize function closes a given collection. It closes the underlying database and if the
     * \ref UF_COLLECTION_CLEANUP flag was specified when initializing the collection the downloaded matrices will
     * be removed from the cache.
     * \sa uf_collection_init
     * \ingroup collection
     */
    void uf_collection_finalize(uf_collection_t * collection);


    /**
     * \brief Cleanup the downloaded matrices.
     * \param[in] collection Collection where the downloaded matrices should be removed.
     *
     * The uf_collection_cleanup function triggers the clean up of a collection independently from if the
     * \ref UF_COLLECTION_CLEANUP flag is set or not. After calling this functions all downloaded matrices
     * belonging to the collection are removed.
     * \sa uf_collection_finalize
     * \ingroup collection
     */
    void uf_collection_cleanup( uf_collection_t * collection);

    /**
     * \brief Return the number of matrices in the collection.
     * \param[in] collection Collection to return the number of matrices from.
     * \return The number of matrices in the collection.
     *
     * The uf_collection_num_matrices function returns the number of matrices in the given collection. In case of an error
     * a value < 0 is returned.
     *
     * \ingroup collection
     */
    int uf_collection_num_matrices(uf_collection_t *collection);

    /* Query  */
    /**
     * \brief Retrieve a matrix by its ID.
     * \param[in] collection Collection where to retrieve the matrix from.
     * \param[in] id  ID of the matrix.
     * \return The matrix object containing the meta information of the matrix.
     *
     * The uf_collection_get_by_id function returns the meta information of a given matrix. The matrix is specified
     * by its unique ID as stored in th UF sparse collection. The matrix is \b not downloaded by this function call.
     * Each returned matrix needs to be freed because the function does not return a pointer
     * to internal data.
     * \sa uf_matrix_get
     * \ingroup query
     * */
    uf_matrix_t * uf_collection_get_by_id(uf_collection_t * collection, int id);

    /**
     * \brief Retrieve a matrix by its name.
     * \param[in] collection Collection where to retrieve the matrix from.
     * \param[in] group_name Name of the group the matrix belongs to.
     * \param[in] name       Name of the matrix.
     * \return The matrix object containing the meta information of the matrix.
     *
     * The uf_collection_get_by_name function returns the meta information of a given matrix. The matrix is specified by
     * its group name and its name inside the group. The matrix is \b not downloaded by this function call.
     * Each returned matrix needs to be freed because the function does not return a pointer
     * to internal data.
     * \sa uf_matrix_get
     *
     * \ingroup query
     * */
    uf_matrix_t * uf_collection_get_by_name(uf_collection_t * collection, const char * group_name, const char *name);

    /**
     * \brief Use a SQL-like query to select matrices.
     * \param[in] collection Collection where to select the matrices from.
     * \param[in] sql   SQL query that can be used inside a WHERE-clause.
     * \return A query to iterate over the matrices matched by the SQL query.
     *
     * The uf_query_sql function allows to create an iterator over a set of matrices. Thereby, the set of matrices is
     * selected via a SQL query. The given string is placed as WHERE clause in the following SQL query:
     * \code{c}
     *  SELECT * FROM matrices WHERE <the_sql_string_passed_to_the_function>;
     * \endcode
     * which is performed on the collection's database. The possible column names are the same as the components of
     * the \ref uf_matrix_t structure where also the description of all of them can be found. The \ref uf_query_next
     * function is used to iterate over the matrices selected by the query. If the query fails, e.g. the SQL query is malicious
     * the returned query object is NULL.
     *
     * ### Examples
     * \li Select all real rectangular matrices:
     * \code{c}
     * uf_query_t * q = uf_query_sql(col, "rows!=cols AND isReal==1");
     * \endcode
     * \li Select all symmetric positive definite matrices:
     * \code{c}
     * uf_query_t * q = uf_query_sql(col, "rows==cols AND isReal==1 AND posdef==1 AND numerical_symmetry==1.0");
     * \endcode
     * \li Select all complex square matrices with more than 10000 non zero elements:
     * \code{c}
     * uf_query_t * q = uf_query_sql(col, "rows==cols AND isComplex==1 AND nnz>10000");
     * \endcode
     *
     * \sa uf_query_next
     * \sa uf_query_free
     *
     * \ingroup query
     * */
    uf_query_t *uf_query_sql(uf_collection_t *collection, const char *sql);

    /**
     * \brief Iterate over a given query.
     * \param[in] query Query to iterate over.
     * \return A \ref uf_matrix_t object containing the next matrix or NULL when the query is empty or finished.
     *
     * The uf_query_next function iterates over a given query. Each call to the function will return the next
     * matrix matched by the query. If no matrices are left or the query does not match any matrix NULL will be
     * returned. The returned \ref uf_matrix_t object contains only the meta information of the matrix as in the
     * return value of \ref uf_collection_get_by_id and \ref uf_collection_get_by_name. The matrices are downloaded
     * using \ref uf_matrix_get. Each returned matrix needs to be freed because the function does not return a pointer
     * to internal data.
     *
     * The following example will download all spd matrices selected by \ref uf_query_spd :
     * \code{.c}
     * uf_query_t *qry = NULL;
     * uf_matrix_t *mat = NULL;
     *
     * qry = uf_query_spd(col);
     *
     * while ( (mat = uf_query_next(qry)) != NULL) {
     *   uf_matrix_get(mat);
     *   uf_matrix_free(mat);
     * }
     *
     * uf_query_free(qry);
     * \endcode
     *
     * \ingroup query
     * */
    uf_matrix_t * uf_query_next(uf_query_t *query);

    /**
     * \brief Free a query.
     * \param[in] query Query to free.
     *
     * The uf_query_free function frees a given query.
     *
     * \ingroup query
     * */
    void  uf_query_free(uf_query_t * query);

    /* Predefined Queries  */
    /**
     * \brief Select all symmetric positive definite matrices.
     * \param[in] collection Collection to select the matrices from.
     * \return A query object for iterating over all spd matrices.
     *
     * The uf_query_spd function is a wrapper around:
     * \code{c}
     * uf_query_sql(col, "rows==cols AND isReal==1 AND posdef==1 AND numerical_symmetry==1.0");
     * \endcode
     *
     * \sa uf_query_sql
     * \ingroup query
     * */
    uf_query_t * uf_query_spd(uf_collection_t *collection);

     /**
     * \brief Select all square matrices with a given value type.
     * \param[in] collection Collection to select the matrices from.
     * \param[in] field Desire value type.
     * \return A query object for iterating over all square matrices.
     *
     * The uf_query_spd function is a wrapper around one of the following function calls depending on the
     * value of field:
     * \li if field == UF_REAL:
     * \code{c}
     * uf_query_sql(col, "rows==cols AND isReal==1");
     * \endcode
     * \li if field == UF_COMPLEX:
     * \code{c}
     * uf_query_sql(col, "rows==cols AND isComplex==1");
     * \endcode
     * \li if field == UF_PATTERN:
     * \code{c}
     * uf_query_sql(col, "rows==cols AND isBinary==1");
     * \endcode
     *
     * \sa uf_query_sql
     * \ingroup query
     */
    uf_query_t * uf_query_square(uf_collection_t *collection, uf_field_t field);

    /**
     * \brief Select all rectangular matrices with a given value type.
     * \param[in] collection Collection to select the matrices from.
     * \param[in] field Desire value type.
     * \return A query object for iterating over all rectangular matrices.
     *
     * The uf_query_spd function is a wrapper around one of the following function calls depending on the
     * value of field:
     * \li if field == UF_REAL:
     * \code{c}
     * uf_query_sql(col, "rows!=cols AND isReal==1");
     * \endcode
     * \li if field == UF_COMPLEX:
     * \code{c}
     * uf_query_sql(col, "rows!=cols AND isComplex==1");
     * \endcode
     * \li if field == UF_PATTERN:
     * \code{c}
     * uf_query_sql(col, "rows!=cols AND isBinary==1");
     * \endcode
     *
     * \sa uf_query_sql
     * \ingroup query
     */
    uf_query_t * uf_query_rectangular(uf_collection_t *collection, uf_field_t field);



    /* Matrix Get Functions */
    /**
     * \brief  Download a given matrix into the cache.
     * \param[in] matrix Matrix entry containing the matrix meta data.
     * \return 0 on success, a non zero value otherwise.
     *
     * The uf_matrix_get function checks if a matrix addressed by its entry form the database
     * is already in the cache. If the matrix is not yet in the cache it is downloaded from
     * the collection.
     *
     * \ingroup matrix
     * */
    int  uf_matrix_get(uf_matrix_t * matrix);

    /**
     * \brief Frees a matrix entry object.
     * \param[in] matrix   Matrix entry to free.
     *
     * The uf_matrix_free function frees a matrix entry. It have to be used on every \ref uf_matrix_t object
     * returned by any uf_collection_* or uf_query_* function in order to avoid memory leaks.
     *
     * \ingroup matrix
     * */
    void uf_matrix_free(uf_matrix_t * matrix);

    /**
     * \brief Print the meta information contained in a matrix entry from the database.
     * \param[in] matrix Matrix entry from the database  to print
     * \param[in] brief Indicator if the output should be brief or not.
     *
     * The uf_matrix_print function prints the meta information of the matrix entry to
     * the standard output. If the brief indicator is true (!=0) the output only contains
     * the most important information.
     *
     * \ingroup matrix
     * */
    void uf_matrix_print(uf_matrix_t * matrix, int brief);

    /**
     * \brief Read a matrix into coordinate storage with 64 bit integers.
     * \param[out] field  Value field of the read matrix.
     * \param[out] nrows  Number of rows of the read matrix.
     * \param[out] ncols  Number of columns of the read matrix.
     * \param[out] nnz    Number of non zeros elements of the read matrix.
     * \param[out] rowptr Array containing the row indices.
     * \param[out] colptr Array containing the column indices.
     * \param[out] values Array containing the values belonging to rowptr and colptr.
     * \param[in]  matrix Matrix entry from the database.
     * \return 0 on success, a non zero error otherwise.
     *
     * The uf_matrix_coord_int64 function takes a matrix entry from the database and return the matrix in
     * coordinate storage. That means the matrix consists of three array of length nnz and the matrix
     * has only non zero entries at \f$ A(rowptr[i], colptr[i]) = values[i]\f$. The arrays
     * rowptr, colptr, and values ares allocate by the function as malloc compatible memory location. Therefore,
     * the have to be freed using the \b free function after usage. The values array is either allocated as
     * an array of sizeof(double)*nnz in the case of a real or a binary matrix or of sizeof(double complex)*nnz
     * in case of a complex matrix. The field output parameters contains again the information of the read matrix
     * is real, binary, or complex valued. In order to avoid compiler warning or errors the function should be call like
     * \code{c}
     *  uf_field_t field;
     *  int64_t nrows, ncols, nnz;
     *  int64_t *rowptr, *colptr;
     *  double *values;
     *  uf_matrix_coord_int64(&field, &nrows, &ncols, &nnz, &rowptr, &colptr, (void **) &values, mat);
     * \endcode
     *
     * \attention
     * The function uses \b int64_t values for all integers. These are normally compatible to the standard \b long
     * data type of C on 64 bit POSIX compatible systems. For 32 bit integers values use
     * \ref uf_matrix_coord_int32 instead.
     * \remark
     * \li If the matrix was not downloaded before using \ref uf_matrix_get it will be downloaded by the call to
     * this function.
     * \li C counting is used, i.e., the row and column indices start with zero.
     *
     * \sa uf_matrix_coord_int32
     *
     * \ingroup matrix
     * */
    int uf_matrix_coord_int64(uf_field_t *field, int64_t *nrows, int64_t *ncols, int64_t *nnz, int64_t **rowptr, int64_t **colptr, void **values, uf_matrix_t *matrix);


    /**
     * \brief Read a matrix into coordinate storage with 32 bit integers.
     * \param[out] field  Value field of the read matrix.
     * \param[out] nrows  Number of rows of the read matrix.
     * \param[out] ncols  Number of columns of the read matrix.
     * \param[out] nnz    Number of non zeros elements of the read matrix.
     * \param[out] rowptr Array containing the row indices.
     * \param[out] colptr Array containing the column indices.
     * \param[out] values Array containing the values belonging to rowptr and colptr.
     * \param[in]  matrix Matrix entry from the database.
     * \return 0 on success, a non zero error otherwise.
     *
     * The uf_matrix_coord_int32 function takes a matrix entry from the database and return the matrix in
     * coordinate storage. That means the matrix consists of three array of length nnz and the matrix
     * has only non zero entries at \f$ A(rowptr[i], colptr[i]) = values[i]\f$. The arrays
     * rowptr, colptr, and values ares allocate by the function as malloc compatible memory location. Therefore,
     * the have to be freed using the \b free function after usage. The values array is either allocated as
     * an array of sizeof(double)*nnz in the case of a real or a binary matrix or of sizeof(double complex)*nnz
     * in case of a complex matrix. The field output parameters contains again the information of the read matrix
     * is real, binary, or complex valued. In order to avoid compiler warning or errors the function should be call like
     * \code{c}
     *  uf_field_t field;
     *  int32_t nrows, ncols, nnz;
     *  int32_t *rowptr, *colptr;
     *  double *values;
     *  uf_matrix_coord_int32(&field, &nrows, &ncols, &nnz, &rowptr, &colptr, (void **) &values, mat);
     * \endcode
     *
     * \attention
     * The function uses \b int32_t values for all integers. These are normally compatible to the standard \b int
     * data type of C. For 64 bit integers values use \ref uf_matrix_coord_int32 instead.
     * \remark
     * \li If the matrix was not downloaded before using \ref uf_matrix_get it will be downloaded by the call to
     * this function.
     * \li C counting is used, i.e., the row and column indices start with zero.
     *
     * \sa uf_matrix_coord_int64
     *
     * \ingroup matrix
     * */
    int uf_matrix_coord_int32(uf_field_t *field, int32_t *nrows, int32_t *ncols, int32_t *nnz, int32_t **rowptr, int32_t **colptr, void **values, uf_matrix_t *matrix);

    /* Conversion between Coord -> {CSR,CSC}  */
    /**
     * \brief Convert the coordinate storage to Compressed Sparse Row(CSR) storage with 32 bit integers.
     * \param[in] field     Value field of the matrix.
     * \param[in] nrows     Number of rows of the matrix.
     * \param[in] ncols     Number of columns of the matrix.
     * \param[in] nnz       Number of non zero elements of the matrix.
     * \param[in,out] rowptr Array containing the row indices on input and the row pointer on output.
     * \param[in,out] colptr Array containing the column indices on input and the column pointer on output.
     * \param[in,out] values Array containing the values on input and the resorted values on output.
     * \return 0 on success, a non zero error otherwise.
     *
     * The uf_matrix_coord_to_csr_int32 function converts a given matrix from coordinate storage to compressed
     * sparse row storage. Thereby, two different algorithms can be used. Either a fast but memory consuming one
     * or a slow but inplace algorithm. By default the fast one is used but it can be change via the
     * \ref uf_matrix_conv_memsave function. The input arrays need to be double pointers because they may be reallocated
     * during the conversion depending on the chosen algorithm and the structure of the matrix. If the arrays are
     * reallocated the new pointer are compatible to malloc as in \ref uf_matrix_coord_int32 and must be freed
     * using \b free after usage.
     *
     * \attention
     * The function uses \b int32_t values for all integers. These are normally compatible to the standard \b int
     * data type of C. For 64 bit integers values use \ref uf_matrix_coord_to_csr_int64 instead.
     *
     * \sa uf_matrix_coord_to_csr_int64
     * \sa uf_matrix_coord_to_csc_int32
     * \sa uf_matrix_coord_to_csc_int64
     *
     * \ingroup matrix
     * */
    int uf_matrix_coord_to_csr_int32(uf_field_t field, int32_t nrows, int32_t ncols, int32_t nnz, int32_t **rowptr, int32_t **colptr, void **values);

    /**
     * \brief Convert the coordinate storage to Compressed Sparse Row(CSR) storage with 64 bit integers.
     * \param[in] field     Value field of the matrix.
     * \param[in] nrows     Number of rows of the matrix.
     * \param[in] ncols     Number of columns of the matrix.
     * \param[in] nnz       Number of non zero elements of the matrix.
     * \param[in,out] rowptr Array containing the row indices on input and the row pointer on output.
     * \param[in,out] colptr Array containing the column indices on input and the column pointer on output.
     * \param[in,out] values Array containing the values on input and the resorted values on output.
     * \return 0 on success, a non zero error otherwise.
     *
     * The uf_matrix_coord_to_csr_int64 function converts a given matrix from coordinate storage to compressed
     * sparse row storage. Thereby, two different algorithms can be used. Either a fast but memory consuming one
     * or a slow but inplace algorithm. By default the fast one is used but it can be change via the
     * \ref uf_matrix_conv_memsave function. The input arrays need to be double pointers because they may be reallocated
     * during the conversion depending on the chosen algorithm and the structure of the matrix. If the arrays are
     * reallocated the new pointer are compatible to malloc as in \ref uf_matrix_coord_int64 and must be freed
     * using \b free after usage.
     *
     * \attention
     * The function uses \b int64_t values for all integers. These are normally compatible to the standard \b long
     * data type of C on 64 bit POSIX compatible systems. For 32 bit integers values use
     * \ref uf_matrix_coord_to_csr_int32 instead.
     *
     * \sa uf_matrix_coord_to_csr_int32
     * \sa uf_matrix_coord_to_csc_int32
     * \sa uf_matrix_coord_to_csc_int64
     *
     * \ingroup matrix
     * */
    int uf_matrix_coord_to_csr_int64(uf_field_t field, int64_t nrows, int64_t ncols, int64_t nnz, int64_t **rowptr, int64_t **colptr, void **values);

    /**
     * \brief Convert the coordinate storage to Compressed Sparse Column (CSC) storage with 32 bit integers.
     * \param[in] field     Value field of the matrix.
     * \param[in] nrows     Number of rows of the matrix.
     * \param[in] ncols     Number of columns of the matrix.
     * \param[in] nnz       Number of non zero elements of the matrix.
     * \param[in,out] rowptr Array containing the row indices on input and the row pointer on output.
     * \param[in,out] colptr Array containing the column indices on input and the column pointer on output.
     * \param[in,out] values Array containing the values on input and the resorted values on output.
     * \return 0 on success, a non zero error otherwise.
     *
     * The uf_matrix_coord_to_csc_int32 function converts a given matrix from coordinate storage to compressed
     * sparse column storage. Thereby, two different algorithms can be used. Either a fast but memory consuming one
     * or a slow but inplace algorithm. By default the fast one is used but it can be change via the
     * \ref uf_matrix_conv_memsave function. The input arrays need to be double pointers because they may be reallocated
     * during the conversion depending on the chosen algorithm and the structure of the matrix. If the arrays are
     * reallocated the new pointer are compatible to malloc as in \ref uf_matrix_coord_int32 and must be freed
     * using \b free after usage.
     *
     * \attention
     * The function uses \b int32_t values for all integers. These are normally compatible to the standard \b int
     * data type of C. For 64 bit integers values use \ref uf_matrix_coord_to_csc_int64 instead.
     *
     * \sa uf_matrix_coord_to_csc_int64
     * \sa uf_matrix_coord_to_csr_int32
     * \sa uf_matrix_coord_to_csr_int64
     *
     * \ingroup matrix
     * */
    int uf_matrix_coord_to_csc_int32(uf_field_t field, int32_t nrows, int32_t ncols, int32_t nnz, int32_t **rowptr, int32_t **colptr, void **values);

    /**
     * \brief Convert the coordinate storage to Compressed Sparse Column(CSC) storage with 64 bit integers.
     * \param[in] field     Value field of the matrix.
     * \param[in] nrows     Number of rows of the matrix.
     * \param[in] ncols     Number of columns of the matrix.
     * \param[in] nnz       Number of non zero elements of the matrix.
     * \param[in,out] rowptr Array containing the row indices on input and the row pointer on output.
     * \param[in,out] colptr Array containing the column indices on input and the column pointer on output.
     * \param[in,out] values Array containing the values on input and the resorted values on output.
     * \return 0 on success, a non zero error otherwise.
     *
     * The uf_matrix_coord_to_csc_int64 function converts a given matrix from coordinate storage to compressed
     * sparse column. Thereby, two different algorithms can be used. Either a fast but memory consuming one
     * or a slow but inplace algorithm. By default the fast one is used but it can be change via the
     * \ref uf_matrix_conv_memsave function. The input arrays need to be double pointers because they may be reallocated
     * during the conversion depending on the chosen algorithm and the structure of the matrix. If the arrays are
     * reallocated the new pointer are compatible to malloc as in \ref uf_matrix_coord_int64 and must be freed
     * using \b free after usage.
     *
     * \attention
     * The function uses \b int64_t values for all integers. These are normally compatible to the standard \b long
     * data type of C on 64 bit POSIX compatible systems. For 32 bit integers values use
     * \ref uf_matrix_coord_to_csr_int32 instead.
     *
     * \sa uf_matrix_coord_to_csr_int32
     * \sa uf_matrix_coord_to_csr_int64
     * \sa uf_matrix_coord_to_csc_int64
     *
     * \ingroup matrix
     * */
    int uf_matrix_coord_to_csc_int64(uf_field_t field, int64_t nrows, int64_t ncols, int64_t nnz, int64_t **rowptr, int64_t **colptr, void **values);

    /**
     * \brief Select the algorithm to convert between Coordinate and Compressed storage.
     * \param[in] save  Boolean value to enable or disable the memory saving conversion algorithm.
     * \return The old value of the setting.
     *
     * The uf_matrix_conv_memsave function selects the algorithm which is used in \ref uf_matrix_coord_to_csc_int32,
     * \ref uf_matrix_coord_to_csc_int64, \ref uf_matrix_coord_to_csr_int32, and \ref uf_matrix_coord_to_csr_int64.
     * By default this algorithm is a fast but memory consuming bucket sort which is linear in the number of
     * non zero elements. If the save parameter is true the algorithm changes to an inplace Quicksort/Insertionsort
     * implementation which is in the O(n log n) runtime class but does not need any additional memory.
     * The return value contains the previous setting of this option before it is set to the new value.
     *
     * \ingroup matrix
     * */
    int uf_matrix_conv_memsave(int save);


#ifdef __cplusplus
};
#endif
#endif /* end of include guard: LIBUFGET_H */



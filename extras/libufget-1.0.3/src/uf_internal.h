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
 * Copyright (C) Martin Koehler, 2015, 2017
 */


#ifndef UF_INTERNAL_H

#define UF_INTERNAL_H

#define UF_INDEX_FILE "files/ss_index.mat"
#define UF_INDEX_URL "http://sparse-files.engr.tamu.edu"
#define UF_CONTAINER_STRUCT "ss_index"
#define UF_SQLITE_TABLE_VERSION 3

#define UFGET_COPYRIGHT_STRING "Copyright (C) 2015, 2016, 2017, 2020 Martin Koehler"

struct _uf_collection_t {
    sqlite3 *db;
    char *baseurl;
    char *cache_dir;
    time_t stamp;
    unsigned long flags;
};
// typedef struct _uf_collection_t uf_collection_t;

struct _uf_query_t {
    uf_collection_t *col;
    sqlite3_stmt *stmt;
    int done;
};
// typedef struct _uf_query_t uf_query_t;

int uf_matrix_set_localpath(uf_matrix_t *mat, const char * localpath);
void uf_matrix_from_stmt(uf_matrix_t *matrix, sqlite3_stmt *res);

int sql_get_long(sqlite3 * db, const char * sql, ...);
int sql_table_exists(sqlite3 * db, const char *table_name);
int sql_get_table(sqlite3 *db, char ***resultp, int *nrow, int *ncol, const char *sql, ...);
int sql_execute(sqlite3 *db, const char* sql, ...);

extern int uf_get_verbose_level;

#define	DPRINTF( level, text, args... )	do { if ( uf_get_verbose_level >= (level)) {fprintf(stderr, text, ## args); } } while(0)


#define MAX(A,B) (((A)<(B))?(B):(A))

/* Column in the SQL table  */
#define COL_ID     0
#define COL_GROUP  1
#define COL_NAME   2
#define COL_NROWS  3
#define COL_NCOLS  4
#define COL_NNZ 5
#define COL_NZERO 6
#define COL_PATTERN_SYMMETRY 7
#define COL_NUMERICAL_SYMMETRY 8
#define COL_ISBINARY 9
#define COL_ISREAL 10
#define COL_ISCOMPLEX 11
#define COL_NNZDIAG 12
#define COL_POSDEF 13
#define COL_AMD_LNZ 14
#define COL_AMD_FLOPS 15
#define COL_AMD_VNZ 16
#define COL_AMD_RNZ 17
#define COL_NBLOCKS 18
#define COL_SPRANK 19
#define COL_RBTYPE 20
#define COL_CHOLCAND 21
#define COL_NCC 22
#define COL_ISND 23
#define COL_ISGRAPH 24
#define COL_LOWERBANDWIDTH 25
#define COL_UPPERBANDWIDTH 26
#define COL_RCM_LOWERBANDWIDTH 27
#define COL_RCM_UPPERBANDWIDTH 28
#define COL_XMIN_REAL 29
#define COL_XMIN_IMAG 30
#define COL_XMAX_REAL 31
#define COL_XMAX_IMAG 32
#define COL_LOCALPATH 33


/* typedef enum {
    MTX_REAL = 0,
    MTX_COMPLEX = 1,
    MTX_PATTERN = 2
} mtx_datatype_t;


#define MTX_IS_REAL(m)		((*m).data_type == MTX_REAL)
#define MTX_IS_INTEGER(m)		((*m).data_type == MTX_REAL)
#define MTX_IS_COMPLEX(m)	((*m).data_type == MTX_COMPLEX)
 */


typedef enum {
	MTX_UNKNOWN = 0,
	MTX_DENSE = 1,
	MTX_COORD = 2
} mtx_storage_t;

#define MTX_IS_DENSE(c)	((*c).store_type == MTX_DENSE)
#define MTX_IS_COORD(c)	((*c).store_type == MTX_COORD)

typedef enum {
	MTX_GENERAL = 0,
	MTX_SYMMETRIC,
	MTX_SKEWSYMMETRIC,
	MTX_HERMITIAN
} mtx_symmetry_t;

#define MTX_IS_GENERAL(m)	((*m).symmetry == MTX_GENERAL)
#define MTX_IS_UNSYMMETRIC(m)	((*m).symmetry == MTX_GENERAL)
#define MTX_IS_SYMMETRIC(m)	((*m).symmetry == MTX_SYMMETRIC)
#define MTX_IS_SKEWSYMMETRIC(m)	((*m).symmetry == MTX_SKEWSYMMETRIC)
#define MTX_IS_HERMITIAN(m)	((*m).symmetry == MTX_HERMITIAN)

typedef struct {
	int64_t cols;
	int64_t rows;
	int64_t nnz;
	mtx_storage_t store_type;
	mtx_symmetry_t symmetry;
	uf_field_t data_type;
} mtx_info_t ;

// struct _uf_io_file_t;
// mtx_info_t * uf_matrix_read_info(struct _uf_io_file_t* file);

#endif /* end of include guard: UF_INTERNAL_H */

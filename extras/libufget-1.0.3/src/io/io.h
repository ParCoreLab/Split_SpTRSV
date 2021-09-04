/* 
* CSCUTILS - A collection of various software routines uses in CSC projects
* Copyright (C) 2015 Martin Koehler
* 
* This library is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published
* by the Free Software Foundation; either version 2.1 of the License, or
* (at your option) any later version.
* 
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU Lesser General Public License for more details.
* 
* You should have received a copy of the GNU Lesser General Public License
* along with this library; if not, see <http://www.gnu.org/licenses/>.
* 
*/ 

#ifndef UF_IO_H
#define UF_IO_H

#ifdef __cplusplus
extern "C" {
#endif
#if __GNUC__ >= 4 
#define CSC_ATTR_SCANF(pos1, pos2)   __attribute__ ((format (scanf, pos1, pos2)));
#else 
#define CSC_ATTR_SCANF(pos1, pos2)  
#endif 


	#include <stdio.h>
	#include <stdlib.h>
	#include <stdarg.h>
	#include <sys/stat.h>

    typedef enum {
		UF_IO_FILE_ERROR = -1,       /**< Indicates that an error occurred during the compresseion detection. */
		UF_IO_FILE_UNCOMPRESSED = 0, /**< Indicates that the file is uncompressed. */
		UF_IO_FILE_GZIP = 1,         /**< Indicates that the file is compressed using gzip. */
		UF_IO_FILE_BZIP2 = 2,        /**< Indicates that the file is compressed using bzip2. */      
		UF_IO_FILE_XZ = 3,            /**< Indicates that the file is compressed using xz(lzma). */
        UF_IO_FILE_HANDLE= 4          /**< Indicates that an external Handle is used. */
	} uf_io_compress_io_type_t;

	typedef enum {
		UF_IO_FILE_READ = 0,   /**< Specifies that the file is opened to read.  */
		UF_IO_FILE_WRITE = 1,  /**< Specifies that the file is opened to write. */
	} uf_io_mode_t;

    typedef struct _compressed_io_handler_t {
		char extension[10];
		uf_io_compress_io_type_t type;
		char magic[10];
		int (*open)(void **data, const char * path, const uf_io_mode_t mode, struct _compressed_io_handler_t *handler);
		int (*close)(void **data);
		size_t (*read)(void *data, void *buf, size_t len);
		size_t (*write)(void *data, void *buf, size_t len);
	} _compressed_io_handler;


	extern  _compressed_io_handler compress_methods[];
	extern int _compressed_io_handler_len;
	extern _compressed_io_handler compress_fallback; 

	void uf_io_init ();

	#define UF_IO_FILE_BUFFER_SIZE 4096

	typedef struct __uf_io_file_t {
		void *handler;                        /**< Handle of the underlying io handler. */
		void *data;                           /**< Auxiliary data for the io handler. */
		char buffer[UF_IO_FILE_BUFFER_SIZE]; /**< Buffer cache of size \ref UF_IO_FILE_BUFFER_SIZE. */
		size_t pos;                           /**< Current position in the buffer cache. */
		size_t avail;                         /**< Available bytes in the buffer cache (Reading) */
		uf_io_mode_t mode;                     /**< Access mode.  */
		int  eof; 	                      /**< Indicate the end of the file */
	} uf_io_file_t;

	
	uf_io_compress_io_type_t  uf_io_is_compressed (const char *path);
	uf_io_file_t* uf_io_open(const char * path, uf_io_mode_t mode);
    int uf_io_flush (uf_io_file_t  *f );
	int uf_io_close (uf_io_file_t  *f );
	int uf_io_eof(uf_io_file_t *f); 
	int uf_io_getc (uf_io_file_t  * f );
	char*   uf_io_gets( char *buf, int len, uf_io_file_t * file );
	ssize_t uf_io_getline( char **buf, size_t *len, uf_io_file_t *file);
	int     uf_io_scanf(uf_io_file_t * f, const char *fmt, ...) CSC_ATTR_SCANF(2,3);
	size_t	uf_io_read(void *ptr, size_t selem, size_t nelem, uf_io_file_t *f);
	int uf_io_putc ( int c, uf_io_file_t *f  );
	int uf_io_puts ( const char *buf, uf_io_file_t * f);
	int uf_io_printf( uf_io_file_t * f, const char * fmt, ...);
	size_t uf_io_write( void * ptr, size_t selem, size_t nelem, uf_io_file_t *f);

	int uf_file_mkdir(const char *dir, mode_t m); 
	int uf_file_exist(const char *path);  


#ifdef __cplusplus
};
#endif



#endif





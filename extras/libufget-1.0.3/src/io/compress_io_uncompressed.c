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

#if _FILE_OFFSET_BITS == 64
#ifndef _LARGEFILE64_SOURCE
#define _LARGEFILE64_SOURCE
#endif
#endif


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

/*-----------------------------------------------------------------------------
 * Memory Mapped IO under Linux   
 *-----------------------------------------------------------------------------*/
#ifdef MMAP_IO 
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h> 
#endif 


#include "io.h"

static int uncompressed_open  (void **data, const char * path, const uf_io_mode_t mode, _compressed_io_handler * handler);
static int uncompressed_open2 (void **data, const char * path, const uf_io_mode_t mode, _compressed_io_handler * handler);
static int uncompressed_close(void **data);
static size_t uncompressed_read(void *data, void *buf, size_t len);
static size_t uncompressed_write(void *data, void *buf, size_t len);



/* Internal structure to represent an uncompressed file. */
typedef struct {
	FILE *fp;
	char * path;
	uf_io_mode_t mode;
}_uncompressed_file;

// Internal structure to represent an uncompressed file. 
typedef struct {
	int fd; 
	size_t size; 
	size_t pos; 
	char * file; 
	char * path; 
	uf_io_mode_t mode;
	int rw; 
}_uncompressed_mmap_file;  

/* Registers the access to uncompressed files. */ 
_compressed_io_handler io_register_uncompressed(){
	_compressed_io_handler handler;
	strcpy(handler.extension, "");
	handler.type = UF_IO_FILE_UNCOMPRESSED;
	strcpy(handler.magic,"");
	handler.open = uncompressed_open;
	handler.close = uncompressed_close;
	handler.write = uncompressed_write;
	handler.read  = uncompressed_read;
	return handler;
}

_compressed_io_handler io_register_no_gzip(){
	_compressed_io_handler handler;
	strcpy(handler.extension, ".gz");
	handler.type = UF_IO_FILE_UNCOMPRESSED;
        strcpy(handler.magic,  "\x1F\x8B\x08"); 
	handler.open = uncompressed_open2;
	handler.close = uncompressed_close;
	handler.write = uncompressed_write;
	handler.read  = NULL;

	return handler;
}

_compressed_io_handler io_register_no_bzip2(){
	_compressed_io_handler handler;
	strcpy(handler.extension, ".bz2");
	handler.type = UF_IO_FILE_UNCOMPRESSED;
	strcpy(handler.magic, "BZh"); 
	handler.open = uncompressed_open2;
	handler.close = uncompressed_close;
	handler.write = uncompressed_write;
	handler.read  = NULL;
	return handler;
}

_compressed_io_handler io_register_no_xz(){
	_compressed_io_handler handler;
	 const char HEADER_MAGIC[6] = { 0xFD, '7', 'z', 'X', 'Z', 0x00 };
	strcpy(handler.extension, ".xz");
	handler.type = UF_IO_FILE_UNCOMPRESSED;
	strcpy(handler.magic,(char *) HEADER_MAGIC);
	handler.open = uncompressed_open2;
	handler.close = uncompressed_close;
	handler.write = uncompressed_write;
	handler.read  = NULL;
	return handler;
}

/* Local function which opens an uncompressed file. */
static int uncompressed_open (void **data, const char * path, const uf_io_mode_t mode, _compressed_io_handler * handler) {
	FILE *fp;
	_uncompressed_file *f;
	int err;

	char _mode[4];
	if ( mode == UF_IO_FILE_WRITE ) {
		strcpy(_mode, "w");
	} else if ( mode == UF_IO_FILE_READ ) {
		strcpy(_mode, "r");
	} else {
		fprintf(stderr, "Wrong mode argument. This must either be UF_IO_FILE_READ or UF_IO_FILE_WRITE\n");
		return -1;
	}
	
	if ( !(fp = fopen(path, _mode))) {
		err = errno;
		fprintf(stderr, "opening file: %s failed, errno: %03d - %s\n", path, err, strerror(err));
		return -1;
	}

	f = (_uncompressed_file * ) malloc( sizeof(_uncompressed_file));
	f->fp = fp;
	f->path = strdup(path);
	f->mode = mode;
	*data = (void *) f;
	return 0;
}

/* Local function which opens an uncompressed file. */ 
static int uncompressed_open2 (void **data, const char * path, const uf_io_mode_t mode, _compressed_io_handler * handler){
	FILE *fp;
	_uncompressed_file *f;
	int err;
	char *pos;
	char *internal_fn;

	char _mode[4];
	if ( mode == UF_IO_FILE_WRITE ) {
		strcpy(_mode, "w");
	} else if ( mode == UF_IO_FILE_READ ) {
		strcpy(_mode, "r");
	} else {
		fprintf(stderr, "Wrong mode argument. This must either be UF_IO_FILE_READ or UF_IO_FILE_WRITE\n");
		return -1;
	}
	/*-----------------------------------------------------------------------------
	 *  dectect extension of the archiver
	 *-----------------------------------------------------------------------------*/
	internal_fn = strdup(path);
	if (strlen(handler->extension) >  0 ) {
		pos = strstr(internal_fn, handler->extension);
	} else {
		pos = NULL;
	}
	if ( pos != NULL ) {
		*pos = '\0';
		fprintf(stderr, "Support for \"%s\" compression not available, truncate to \"%s\".", handler->extension, internal_fn);
	}


	if ( !(fp = fopen(internal_fn, _mode))) {
		free(internal_fn);
		err = errno;
		fprintf(stderr, "opening file: %s failed, errno: %03d - %s\n", path, err, strerror(err));

		return -1;
	}

	f = (_uncompressed_file * ) malloc( sizeof(_uncompressed_file));
	f->fp = fp;
	f->path = strdup(path);
	f->mode = mode;
	*data = (void *) f;
	free(internal_fn);
	return 0;
}


/* Local function to close and free an uncompressed file */ 
static int uncompressed_close(void **data) {
	_uncompressed_file *f;

	if ( data == NULL ) {
		fprintf(stderr, "Error: data == NULL\n");
		return -1;
	}
	if ( *data == NULL ) {
		fprintf(stderr, "Error: *data == NULL\n");
		return -1;
	}

	f = (_uncompressed_file*) *data;
	fclose(f->fp);
	free(f->path);
	free(f);
	*data = NULL;
	return 0;
}


/* Local function to read a block from an uncompressed file */ 
static size_t uncompressed_read(void *data, void *buf, size_t len) {
	_uncompressed_file *f;
	if ( data == NULL) return -1;
	if ( buf == NULL) return -1;
	f = (_uncompressed_file*) data;
	return fread(buf,1,len,f->fp);
}

/* Local function to write a block to an uncompressed file. */ 
static size_t uncompressed_write(void *data, void *buf, size_t len){
	_uncompressed_file *f;
	if ( data == NULL) return -1;
	if ( buf == NULL) return -1;
	f = (_uncompressed_file*) data;
	return (size_t) fwrite(buf,1,len,f->fp);
}


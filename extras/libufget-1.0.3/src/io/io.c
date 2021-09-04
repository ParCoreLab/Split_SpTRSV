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

#if _FILE_OFFSET_BITS == 64 && !defined(_LARGEFILE64_SOURCE)
#define _LARGEFILE64_SOURCE
#endif

#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdarg.h>
#include <stdint.h>
#include <malloc.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "io.h"


#define MIN(A,B) (((A)<(B))?(A):(B))
#define MAX(A,B) (((A)>(B))?(A):(B))

extern void uf_io_init();

char *repl_str(const char *str, const char *old, const char *new) {

    /* Adjust each of the below values to suit your needs. */

    /* Increment positions cache size initially by this number. */
    size_t cache_sz_inc = 16;
    /* Thereafter, each time capacity needs to be increased,
     * multiply the increment by this factor. */
    const size_t cache_sz_inc_factor = 3;
    /* But never increment capacity by more than this number. */
    const size_t cache_sz_inc_max = 1048576;

    char *pret, *ret = NULL;
    const char *pstr2, *pstr = str;
    size_t i, count = 0;
    ptrdiff_t *pos_cache = NULL;
    size_t cache_sz = 0;
    size_t cpylen, orglen, retlen, newlen, oldlen = strlen(old);

    /* Find all matches and cache their positions. */
    while ((pstr2 = strstr(pstr, old)) != NULL) {
        count++;

        /* Increase the cache size when necessary. */
        if (cache_sz < count) {
            cache_sz += cache_sz_inc;
            pos_cache = realloc(pos_cache, sizeof(*pos_cache) * cache_sz);
            if (pos_cache == NULL) {
                goto end_repl_str;
            }
            cache_sz_inc *= cache_sz_inc_factor;
            if (cache_sz_inc > cache_sz_inc_max) {
                cache_sz_inc = cache_sz_inc_max;
            }
        }

        pos_cache[count-1] = pstr2 - str;
        pstr = pstr2 + oldlen;
    }

    orglen = pstr - str + strlen(pstr);

    /* Allocate memory for the post-replacement string. */
    if (count > 0) {
        newlen = strlen(new);
        retlen = orglen + (newlen - oldlen) * count;
    } else  retlen = orglen;
    ret = malloc(retlen + 1);
    if (ret == NULL) {
        goto end_repl_str;
    }

    if (count == 0) {
        /* If no matches, then just duplicate the string. */
        strcpy(ret, str);
    } else {
        /* Otherwise, duplicate the string whilst performing
         * the replacements using the position cache. */
        pret = ret;
        memcpy(pret, str, pos_cache[0]);
        pret += pos_cache[0];
        for (i = 0; i < count; i++) {
            memcpy(pret, new, newlen);
            pret += newlen;
            pstr = str + pos_cache[i] + oldlen;
            cpylen = (i == count-1 ? orglen : pos_cache[i+1]) - pos_cache[i] - oldlen;
            memcpy(pret, pstr, cpylen);
            pret += cpylen;
        }
        ret[retlen] = '\0';
    }

end_repl_str:
    /* Free the cache and return the post-replacement string,
     * which will be NULL in the event of an error. */
    free(pos_cache);
    return ret;
}

/* Detect the compression type of a file  */
uf_io_compress_io_type_t  uf_io_is_compressed (const char *path)
{
	FILE *fp;
	char head[10];
	size_t read;
	int i;

	/*-----------------------------------------------------------------------------
	 *  open the file
	 *-----------------------------------------------------------------------------*/
	fp = fopen (path, "r");
	if (!fp) {
		for (i = 0 ; i < _compressed_io_handler_len;  i++) {
			if ( strlen(compress_methods[i].extension)  > 0  && strstr(path, compress_methods[i].extension) != NULL){
				return compress_methods[i].type;
			}
		}
		return UF_IO_FILE_UNCOMPRESSED;
	} else {
		/*-----------------------------------------------------------------------------
		 *  read the header
		 *-----------------------------------------------------------------------------*/
		read=fread(head, sizeof ( char ) ,  10 , fp);
		fclose(fp);

		for ( i = 0; i < _compressed_io_handler_len; i++){
			if ( strlen(compress_methods[i].magic) > 0 && strncmp(head, compress_methods[i].magic,MIN(strlen(compress_methods[i].magic),read)) == 0 ) {
				return compress_methods[i].type;
			}
		}
		return UF_IO_FILE_UNCOMPRESSED;
	}
	return UF_IO_FILE_ERROR;
}

/* Get the compression handler  */
static _compressed_io_handler * get_compress_handler(const char *path) {
	FILE *fp;
	int err = 0;
	char head[10];
	size_t read;
	int i;
	_compressed_io_handler *handler = NULL;

	/*-----------------------------------------------------------------------------
	 *  open the file
	 *-----------------------------------------------------------------------------*/
	fp = fopen (path, "r");
	if (!fp) {
		err = errno;
		fprintf(stderr, "opening file: %s failed, errno: %03d - %s\n", path, err, strerror(err));
		return NULL;
	}

	/*-----------------------------------------------------------------------------
	 *  read the header
	 *-----------------------------------------------------------------------------*/
	read=fread(head, sizeof ( char ) ,  10 , fp);
	fclose(fp);

	for ( i = 0; i < _compressed_io_handler_len; i++){

		if ( strlen(compress_methods[i].magic) > 0  && strncmp(head, compress_methods[i].magic,MIN(strlen(compress_methods[i].magic),read)) == 0 ) {
			handler =  &(compress_methods[i]);
			break;
		}
	}
	if ( handler == NULL ) {
		handler = &compress_fallback;
	}
	return handler;
}


/* Open a file   */
uf_io_file_t* uf_io_open(const char * path, uf_io_mode_t mode)
{
	int rw = 0;
	int ret = 0;
	int i = 0;
	uf_io_file_t  *f;

	if ( path == NULL) {
		fprintf(stderr, "path points to NULL\n");
	}


	if ( mode != UF_IO_FILE_READ && mode != UF_IO_FILE_WRITE ) {
		fprintf(stderr, "The mode parameter must either be UF_IO_FILE_WRITE or UF_IO_FILE_READ.");
		return NULL;
	}

	uf_io_init ();

	/*-----------------------------------------------------------------------------
	 *  detect write mode
	 *-----------------------------------------------------------------------------*/
	if ( mode == UF_IO_FILE_WRITE ) {
		rw = 1;
	}

	f = (uf_io_file_t *) malloc ( sizeof(uf_io_file_t));
	if ( f == NULL ) {
		fprintf(stderr, "Allocation of the mess_file object failed.\n");
		return NULL;
	}
	/*-----------------------------------------------------------------------------
	 *  Open the file
	 *-----------------------------------------------------------------------------*/
	if ( rw == 0) {
		_compressed_io_handler *handler;
		f->pos = 0;
		f->avail = 0;
		f->data  = NULL;
		f->mode = UF_IO_FILE_READ;
		f->eof = 0;
		handler = get_compress_handler(path);
		if ( !handler ) {
			fprintf(stderr, "Can not determine file IO hanlder. \n");
			free(f);
			return NULL;
		}

		if ( handler->read == NULL) {
			fprintf(stderr, "File \"%s\" cannot be opened for reading. The \"%s\" archive extension is not supported.\n", path, handler->extension);
			free(f);
			return NULL;
		}
		f->handler =(void *) handler;
		ret = handler->open(&(f->data), path, mode, handler);
		if ( ret != 0) {
		    fprintf(stderr, "Opening %s failed.\n", path);
			free(f);
			return NULL;
		}
	} else {
		_compressed_io_handler *handler = NULL;
		_compressed_io_handler *fallback = NULL;
		for (i = 0 ; i < _compressed_io_handler_len;  i++) {
			if ( strlen(compress_methods[i].extension)  > 0  && strstr(path, compress_methods[i].extension) != NULL){
				handler = &(compress_methods[i]);
				break;
			}
			if ( strlen(compress_methods[i].extension) == 0) {
				fallback = &(compress_methods[i]);
			}
		}
		if ( handler == NULL) handler= fallback;
		if ( handler == NULL) {
			fprintf(stderr, "No suitable IO handler for %s found.\n", path);
			free(f);
			return NULL;
		}
		if ( handler->write == NULL ) {
			fprintf(stderr, "The select archive format \"%s\" does not support writing.\n", path );
			free(f);
			return NULL;
		}

		f->mode=UF_IO_FILE_WRITE;
		f->pos = 0;
		f->avail = 0;
		f->data  = NULL;
		f->eof = 0;
		f->handler =(void *) handler;
		ret = handler->open(&(f->data), path, mode, handler);
		if ( ret != 0) {
			fprintf(stderr, "Opening %s failed.\n", path);

			free(f);
			return NULL;
		}
	}
	return f;
}

/* Flush the buffer cache   */
int uf_io_flush (uf_io_file_t  *f )
{
	size_t w =0;
	void *htmp;
	struct _compressed_io_handler_t * handler;

	if (!f) return -1;
	if (f->mode != UF_IO_FILE_WRITE) return -1;

	htmp = (void *) (f->handler);
	handler = (struct _compressed_io_handler_t *) htmp;

	if (f->pos == 0 ) return 0;
	w = handler->write(f->data, f->buffer, f->pos);
	if ( w != f->pos ) {
		return -1;
	}
	f->pos = 0;
	return 0;
}

/* Close the file */
int uf_io_close (uf_io_file_t  *f )
{
	_compressed_io_handler *h;
	if ( f == NULL ) {
	    fprintf(stderr, "File f points to NULLs\n");
		return -1;
	}

	if ( f->mode == UF_IO_FILE_WRITE) {
		uf_io_flush(f);
	}
	h = (_compressed_io_handler *) f->handler;
	h->close(&(f->data));
	free(f);
	return 0;

}

/* Check for EOF  */
int uf_io_eof(uf_io_file_t *f){
	if (f == NULL )
		return 1;
	if (f->eof)
		return 1;
	return 0;
}

/* Read a character from a file    */
int     uf_io_getc (uf_io_file_t  * f )
{
	void *htmp;
	struct _compressed_io_handler_t * handler;

	if (!f) return EOF;
	if (f->mode != UF_IO_FILE_READ) return -1;

	htmp = (void *) (f->handler);
	handler = (struct _compressed_io_handler_t *) htmp;

	if ( handler->read == NULL ) {
		fprintf(stderr, "The used io-handler does not provide any read function.\n");
		return -1;
	}

	if ( f -> pos == f->avail ) {
		f->avail = handler->read(f->data, f->buffer, UF_IO_FILE_BUFFER_SIZE);
        // printf("f->avail: %ld\n", f->avail);
		if ( f->avail <= 0) {
			f->eof = 1;
			return -1;
		}
		f->pos = 0;
	}
    // printf("avail: %ld\npos = %ld\n", f->avail,f->pos);
    if ( f->avail <= 0 ) {
        return -1;
    }
	if ( f -> pos < f->avail) {
		return (int) f->buffer[f->pos++];
	}
	return -1;
}

/* Read a string from the file */
char*   uf_io_gets( char *buf, int len, uf_io_file_t * f)
{
	char *buf_sic = buf;
	int c;

	if ( f == NULL ) {
		fprintf(stderr, "File f points to NULLs\n");
		return NULL;
	}
	if ( buf == NULL) {
		fprintf(stderr, "Buffer buf points to NULL\n");
		return NULL;
	}

	if ( len <= 0) return NULL;
	memset(buf, 0, len);
	while ( --len > 0 ) {
		c = uf_io_getc(f);
		if ( c < 0 ) {
			f->eof = 1;
			return buf_sic;
		}
		if ( c == '\n' || c=='\r') {
			if ( c == '\r') c = uf_io_getc(f);
			*buf = '\0';
			return buf_sic;
		}
		*buf = (char) c;
		buf++;
	}
	return buf_sic;

}

#define CHUNK_SIZE 128
/* Replacement for getline */
ssize_t uf_io_getline( char **buf, size_t *len, uf_io_file_t *file)
{
	int idx = 0;
	int c;

	if ( buf  == NULL || len == NULL || file == NULL)  {
		fprintf(stderr, "Either buf, len or file points to NULL");
		return -1;
	}

    if (file->eof) {
        *len = 0;
        return 0;
    }

	if (*buf == NULL)
	{
		*buf = malloc (CHUNK_SIZE);
		if (*buf == NULL) {
			return -1;
		}
		*len = CHUNK_SIZE;
	}

	while (1)
	{
        c = uf_io_getc (file);
        if ( c < 0 )
            break;
		if (idx >= *len)
		{
			*buf = realloc (*buf, *len + CHUNK_SIZE);
			if (*buf == NULL) {
				return -1;
			}
			*len += CHUNK_SIZE;
		}

		(*buf)[idx++] = c;

		if (c == '\n')
			break;
	}
    if (c < 0) {
        file->eof =1;
    }

	if (idx >= *len)
	{
			*buf = realloc (*buf, *len + 1);
			if (*buf == NULL) {
				return -1;
			}
			*len += 1;
	}
	/* Null terminate the buffer.  */
	(*buf)[idx++] = 0;

	return (c == EOF && (idx - 1) == 0) ? -1 : idx - 1;
}


#ifdef HAVE_VSSCANF
/* Fscanf replacement */
int     uf_io_scanf(uf_io_file_t * f, const char *fmt, ...)
{
	int ret = 0 ;
	char *buf = NULL;
	size_t len = 0 ;
	ssize_t iret = 0;
	va_list args;

	iret = uf_io_getline(&buf, &len, f);
	if ( iret < 0 ) {
		if (buf != NULL)
			free(buf);
		return -1;
	}
	va_start(args,fmt);
	ret = vsscanf(buf, fmt, args);
	va_end(args);

	free(buf);
	return ret;
}
#endif


/* Binary Read   */
size_t	uf_io_read(void *ptr, size_t selem, size_t nelem, uf_io_file_t *f)
{
	void *htmp;
	struct _compressed_io_handler_t * handler;
	size_t nbytes = selem * nelem;
	size_t remain_buffer, pos_read = 0;
	char *_ptr = ptr;

	if (!f) return 0;
	if ( selem == 0 ) return 0;
	if ( nelem == 0 ) return 0;
	if (f->mode != UF_IO_FILE_READ) return 0;

	htmp = (void *) (f->handler);
	handler = (struct _compressed_io_handler_t *) htmp;

	if ( handler->read == NULL ) {
		fprintf(stderr, "The used io-handler does not provide any read function.\n");
		return -1;
	}

	remain_buffer = f->avail-f->pos;

	if ( remain_buffer > nbytes ) {
		memcpy(_ptr, &(f->buffer[f->pos]), nbytes);
		f->pos += nbytes;
		pos_read += nbytes;
	} else {
		while (nbytes > 0 ) {
			size_t len;
			remain_buffer = f->avail-f->pos;
			len = MIN(remain_buffer, nbytes);
			memcpy(&_ptr[pos_read], &(f->buffer[f->pos]), len);
			pos_read += len;
			nbytes -= len;
			f->pos += len;
			if ( f -> pos == f->avail ) {
				f->avail = handler->read(f->data, f->buffer, UF_IO_FILE_BUFFER_SIZE);
				if ( f->avail <= 0) {
					f->eof = 1;
				}
				f->pos = 0;
			}
		}

	}
	return (pos_read/selem);

}

/* Write a character */
int uf_io_putc ( int c, uf_io_file_t *f  )
{
	void *htmp;
	struct _compressed_io_handler_t * handler;

	size_t w = 0;
	if (!f) return EOF;
	if (f->mode != UF_IO_FILE_WRITE) return EOF;

	htmp = (void *) (f->handler);
	handler = (struct _compressed_io_handler_t *) htmp;

	if ( handler->write == NULL ) {
		fprintf(stderr, "The used io-handler does not provide any write function.\n");
		return -1;
	}

	f->buffer[f->pos++] = (char) c;
	if ( f -> pos == UF_IO_FILE_BUFFER_SIZE ) {
		f->pos = 0;
		w = handler->write(f->data, f->buffer, UF_IO_FILE_BUFFER_SIZE);
		if ( w != UF_IO_FILE_BUFFER_SIZE ) {
			return EOF;
		}
	}
	return c;
}

/* Write a string  */
int uf_io_puts ( const char *buf, uf_io_file_t * f)
{
	_compressed_io_handler * h;
	size_t len, free_in_buffer;
	size_t to_write;
	size_t w;
	size_t written = 0;

	if ( f == NULL ) return -1;
	if ( buf == NULL) return -1;
	if ( f->mode != UF_IO_FILE_WRITE )  return -1;


	h = ( _compressed_io_handler *) f->handler;
	if ( h->write == NULL ) {
		fprintf(stderr, "The used io-handler does not provide any write function.\n");
		return -1;
	}


	len = strlen(buf);

	free_in_buffer = UF_IO_FILE_BUFFER_SIZE - f->pos;
	if ( len < free_in_buffer ) {
		memcpy(&f->buffer[f->pos],buf, len);
		f->pos += len;
		written = len;
	} else {
		while (len > 0 ) {
			free_in_buffer = UF_IO_FILE_BUFFER_SIZE - f->pos;
			to_write = MIN (free_in_buffer, len);
			memcpy(&f->buffer[f->pos],buf, to_write );
			f->pos += to_write;
			if ( f -> pos == UF_IO_FILE_BUFFER_SIZE ) {
				f->pos = 0;
				w = h->write(f->data, f->buffer, UF_IO_FILE_BUFFER_SIZE);
				if ( w != UF_IO_FILE_BUFFER_SIZE ) {
					fprintf(stderr, "Writing Buffer to File failed. Only %ld of %ld bytes written.\n", (long) w, (long) UF_IO_FILE_BUFFER_SIZE);
					return -1;
				}
			}
			buf += to_write;
			len -= to_write;
			written += to_write;
		}
	}
	return written;
}

/* Formated printing   */
#ifdef HAVE_VSNPRINTF
static char * make_message(const char *fmt, va_list ap)
{
	int n;
	int size = 100;     /* Guess we need no more than 100 bytes. */
	char *p, *np;

	if ((p = malloc(size)) == NULL)
		return NULL;

	while (1) {

		/* Try to print in the allocated space. */
		// va_start(ap, fmt);
		n = vsnprintf(p, size, fmt, ap);
		// va_end(ap);

		/* If that worked, return the string. */
		if (n > -1 && n < size)
			return p;

		/* Else try again with more space. */

		if (n > -1)    /* glibc 2.1 */
			size = n+1; /* precisely what is needed */
		else           /* glibc 2.0 */
			size *= 2;  /* twice the old size */

		if ((np = realloc (p, size)) == NULL) {
			free(p);
			return NULL;
		} else {
			p = np;
		}
	}
	return NULL;
}

int uf_io_printf( uf_io_file_t * f, const char * fmt, ...)
{
	va_list args;
	char *str;
	int ret = 0;

	if ( f == NULL) return -1;
	if ( fmt == NULL) return -1;

	va_start(args, fmt);
	str = make_message(fmt, args);
	va_end(args);
	if ( str == NULL )
		return -1;
	ret = uf_io_puts(str,f );
	free(str);
	return ret;
}


#endif

/* Binary Write */
size_t uf_io_write( void * ptr, size_t selem, size_t nelem, uf_io_file_t *f)
{
	_compressed_io_handler * h;
	size_t len, free_in_buffer;
	size_t to_write;
	size_t w;
	size_t written = 0;
	char *buf = ptr;

	if ( f == NULL ) return 0;
	if ( buf == NULL) return 0;
	if ( f->mode != UF_IO_FILE_WRITE )  return 0;
	if ( selem == 0 ) return 0;


	h = ( _compressed_io_handler *) f->handler;
	if ( h->write == NULL ) {
		fprintf(stderr, "The used io-handler does not provide any write function.\n");
		return 0;
	}
	len = selem * nelem;

	free_in_buffer = UF_IO_FILE_BUFFER_SIZE - f->pos;
	if ( len < free_in_buffer ) {
		memcpy(&f->buffer[f->pos],buf, len);
		f->pos += len;
		written = len;
	} else {
		while (len > 0 ) {
			free_in_buffer = UF_IO_FILE_BUFFER_SIZE - f->pos;
			to_write = MIN (free_in_buffer, len);
			memcpy(&f->buffer[f->pos],buf, to_write );
			f->pos += to_write;
			if ( f -> pos == UF_IO_FILE_BUFFER_SIZE ) {
				f->pos = 0;
				w = h->write(f->data, f->buffer, UF_IO_FILE_BUFFER_SIZE);
				if ( w != UF_IO_FILE_BUFFER_SIZE ) {
					fprintf(stderr, "Writing Buffer to File failed. Only %ld of %ld bytes written.\n", (long) w, (long) UF_IO_FILE_BUFFER_SIZE);
					return -1;
				}
			}
			buf += to_write;
			len -= to_write;
			written += to_write;
		}
	}
	return written/selem;

}

#include <errno.h>

int uf_file_mkdir(const char *dir, mode_t m)
{
    char *tmp = NULL;
    char *p = NULL;
    struct stat dir_stat;
    size_t len;
    int ret ;

    tmp = repl_str(dir, "//","/");
    len = strlen(tmp);
    if(tmp[len - 1] == '/')
        tmp[len - 1] = 0;
    for(p = tmp + 1; *p; p++)
        if(*p == '/') {
            *p = 0;
            // printf("tmp: %s\n", tmp);
            ret = stat(tmp, &dir_stat);
            if (ret == 0 && S_ISDIR(dir_stat.st_mode) ) {
                    *p = '/';
                    continue;
            }
            ret = mkdir(tmp, m);
            if ( ret != 0 ) {
                fprintf(stderr, "Failed to create %s\n", tmp);
                free(tmp);
                return ret;
            }
            *p = '/';

        }
    // printf("tmp: %s\n", tmp);
    ret = stat(tmp, &dir_stat);
    if (ret == 0 && S_ISDIR(dir_stat.st_mode) ) {
        free(tmp);
        return 0;
    }

    ret = mkdir(tmp, m);
    free(tmp);
    return ret;
}

int uf_file_exist(const char *path) {
	if( access( path, F_OK ) != -1 ) {
		return -1;
	} else {
		return 0;
	}
}

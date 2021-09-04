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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include <archive.h>
#include <archive_entry.h>

int main(int argc, char **argv)
{
   struct archive *a;
   struct archive *ext;
   struct archive_entry         *entry;
   int r;
   size_t size; 
   int64_t offset; 
   const void *buffer; 
   char *cbuffer; 

   a = archive_read_new();
   archive_read_support_format_raw(a);
   archive_read_support_format_all(a);

   archive_read_support_filter_all(a); 
  
   if ((r = archive_read_open_filename(a, "test.gz", 10240))){
       archive_read_free(a);
       fprintf(stderr, "Failed to open %s\n", archive_error_string(a));
       return -1; 
   }
   
   while (1 ) {
       int64_t s; 
        r = archive_read_next_header(a, &entry); 
        if ( r == ARCHIVE_EOF) 
            break;
        if ( r != ARCHIVE_OK) {
            fprintf(stderr, "ERROR: %s\n", archive_error_string(a));
            break; 
        }

        // s = archive_entry_size(entry); 
        if ( strcmp(archive_entry_pathname(entry), internal_filename) == 0) {
            fprintf(stderr, "Found %s in file.\n", internal_filename);
        }

        printf("FILE: %s \t %ld\n", archive_entry_pathname(entry), (long) s ); 
        /* printf("Content\n");
        while ( ( r = archive_read_data_block(a, &buffer, &size, &offset)) == ARCHIVE_OK) {
            printf("size = %lu offset = %ld\n", size, offset);
            cbuffer = (char*) buffer; 
            printf("buffer: %s\n", cbuffer);
        } */

   }

   archive_read_close(a);
   archive_read_free(a);
   return 0;
}

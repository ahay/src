#ifndef __jpeg_h
#define __jpeg_h

#include <stdio.h>
#include <jpeglib.h>
#include <rsf.h>

GLOBAL(void)
write_JPEG_file (JSAMPLE * image_buffer, 
		 int image_height, int image_width, int numcol);
GLOBAL(void) read_JPEG_file (sf_file out);

#endif

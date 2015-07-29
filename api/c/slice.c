/* Read slices from files */
/*
  Copyright (C) 2006 Colorado School of Mines
  Copyright (C) 2004 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE
#endif
#include <sys/types.h>
#include <unistd.h>

#include <stdio.h>
#include <stdlib.h>

#include "file.h"
/*^*/

#include "slice.h"
#include "alloc.h"
#include "error.h"

#ifndef _sf_slice_h

typedef struct sf_Slice *sf_slice;
/* abstract data type */
/*^*/

typedef struct sf_FSlice *sf_fslice;
/* abstract data type */
/*^*/

#endif

/*------------------------------------------------------------*/
struct sf_Slice {
    off_t start;
    sf_file file;
    int n12, n3;
};

struct sf_FSlice {
    FILE* file;
    char* name;
    off_t n1;
    int   n2;
};
/*------------------------------------------------------------*/

sf_slice sf_slice_init(sf_file file, int n1, int n2, int n3)
/*< initialize a sliceable file object >*/
{
    sf_slice sl;

    sl = (sf_slice) sf_alloc(1,sizeof(*sl));
    sl->file = file;
    sl->n12 = n1*n2;
    sl->n3 = n3;

    sf_unpipe(sl->file,n1*n2*n3*sizeof(float));
    sl->start = sf_tell(sl->file);
    
    return sl;
}

/*------------------------------------------------------------*/
void sf_slice_get(sf_slice sl, int i3, float* data)
/*< get a slice at level i3 >*/
{
    sf_seek(sl->file,sl->start+i3*(sl->n12)*sizeof(float),SEEK_SET);
    sf_floatread(data,sl->n12,sl->file);
}

/*------------------------------------------------------------*/
void sf_slice_put(sf_slice sl, int i3, float* data)
/*< put a slice at level i3 >*/
{
    sf_seek(sl->file,sl->start+i3*(sl->n12)*sizeof(float),SEEK_SET);
    sf_floatwrite(data,sl->n12,sl->file);
}

/*------------------------------------------------------------*/
void sf_cslice_get(sf_slice sl, int i3, sf_complex* data)
/*< get a slice at level i3 >*/
{
    sf_seek(sl->file,sl->start+i3*(sl->n12)*sizeof(sf_complex),SEEK_SET);
    sf_complexread(data,sl->n12,sl->file);
}

/*------------------------------------------------------------*/
void sf_cslice_put(sf_slice sl, int i3, sf_complex* data)
/*< put a slice at level i3 >*/
{
    sf_seek(sl->file,sl->start+i3*(sl->n12)*sizeof(sf_complex),SEEK_SET);
    sf_complexwrite(data,sl->n12,sl->file);
}

/*------------------------------------------------------------*/
sf_fslice sf_fslice_init(int n1, int n2, size_t size)
/*< initialize a sliceable file object >*/
{
    sf_fslice sl;

    sl = (sf_fslice) sf_alloc(1,sizeof(*sl));
    sl->file = sf_tempfile(&(sl->name), "w+b");
    sl->n1 = n1*size;
    sl->n2 = n2;
    
    return sl;
}

/*------------------------------------------------------------*/
void sf_fslice_get(sf_fslice sl, int i2, void* data)
/*< get a slice at level i2 >*/
{
    extern int fseeko(FILE *stream, off_t offset, int whence);

    if (0 > fseeko(sl->file,i2*(sl->n1),SEEK_SET))
	sf_error ("%s: seek problem:",__FILE__);
    if (sl->n1 != fread(data,1,sl->n1,sl->file))
	sf_error ("%s: read problem:",__FILE__);
}

/*------------------------------------------------------------*/
void sf_fslice_put(sf_fslice sl, int i2, void* data)
/*< put a slice at level i2 >*/
{
    extern int fseeko(FILE *stream, off_t offset, int whence);

    if (0 > fseeko(sl->file,i2*(sl->n1),SEEK_SET))
	sf_error ("%s: seek problem:",__FILE__);
    if (sl-> n1 != fwrite(data,1,sl->n1,sl->file))
	sf_error ("%s: write problem:",__FILE__);
}

/*------------------------------------------------------------*/
void sf_fslice_close(sf_fslice sl)
/*< remove the file and free allocated storage >*/
{
    unlink(sl->name);
    free(sl);
}

/*------------------------------------------------------------*/
void sf_fslice_load(sf_file in, sf_fslice sl, sf_datatype type)
/*< load the contents of in into a slice file >*/
{
    off_t nleft, n;
    char buf[BUFSIZ];
    extern int fseeko(FILE *stream, off_t offset, int whence);

    /* data size */
    n = (sl->n2)*(sl->n1);

    /* rewind */
    if (0 != fseeko(sl->file,0,SEEK_SET))
	sf_error("%s: seeking error:",__FILE__);

    for (nleft = BUFSIZ; n > 0; n -= nleft) {
	if (nleft > n) nleft=n;
	switch (type) {
	    case SF_FLOAT:
		sf_floatread  ((float*) buf,
			       nleft/sizeof(float),in);
		break;
	    case SF_COMPLEX:
		sf_complexread((sf_complex*) buf,
			       nleft/sizeof(sf_complex),in);
		break;
	    default:
		sf_error("%s: unsupported type %d",__FILE__,type);
		break;
	}
	if (nleft != fwrite (buf,1,nleft,sl->file))
	    sf_error("%s: writing error:",__FILE__);
    }
}

/*------------------------------------------------------------*/
void sf_fslice_dump(sf_file out, sf_fslice sl, sf_datatype type)
/*< dump the contents of a slice file to out >*/
{
    off_t nleft, n;
    char buf[BUFSIZ];
    extern int fseeko(FILE *stream, off_t offset, int whence);

    /* data size */
    n = (sl->n2)*(sl->n1);

    /* rewind */
    if (0 != fseeko(sl->file,0,SEEK_SET))
	sf_error("%s: seeking error:",__FILE__);

    for (nleft = BUFSIZ; n > 0; n -= nleft) {
	if (nleft > n) nleft=n;
	if (nleft != fread (buf,1,nleft,sl->file))
	    sf_error("%s: reading error:",__FILE__);
	switch (type) {
	    case SF_FLOAT:
		sf_floatwrite  ((float*) buf,
				nleft/sizeof(float),out);
		break;
	    case SF_COMPLEX:
		sf_complexwrite((sf_complex*) buf,
				nleft/sizeof(sf_complex),out);
		break;
	    default:
		sf_error("%s: unsupported type %d",__FILE__,type);
		break;
	}
    }
}

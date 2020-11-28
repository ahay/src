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

#include <rsf.h>
/*^*/

#include "slice.h"

#ifndef _slice_h

typedef struct Slice *slice;
/* abstract data type */
/*^*/

typedef struct FSlice *fslice;
/* abstract data type */
/*^*/

#endif

struct Slice {
    off_t start;
    sf_file file;
    int n12, n3;
};

struct FSlice {
    FILE* file;
    char* name;
    off_t n1;
    int   n2;
};


slice slice_init(sf_file file, int n1, int n2, int n3)
/*< initialize a sliceable file object >*/
{
    slice sl;

    sl = (slice) sf_alloc(1,sizeof(*sl));
    sl->file = file;
    sl->n12 = n1*n2;
    sl->n3 = n3;

    sf_unpipe(sl->file,n1*n2*n3*sizeof(float));
    sl->start = sf_tell(sl->file);
    
    return sl;
}

void slice_get(slice sl, int i3, float* data)
/*< get a slice at level i3 >*/
{
    sf_seek(sl->file,sl->start+i3*(sl->n12)*sizeof(float),SEEK_SET);
    sf_floatread(data,sl->n12,sl->file);
}

void slice_put(slice sl, int i3, float* data)
/*< put a slice at level i3 >*/
{
    sf_seek(sl->file,sl->start+i3*(sl->n12)*sizeof(float),SEEK_SET);
    sf_floatwrite(data,sl->n12,sl->file);
}

void cslice_get(slice sl, int i3, sf_complex* data)
/*< get a slice at level i3 >*/
{
    sf_seek(sl->file,sl->start+i3*(sl->n12)*sizeof(sf_complex),SEEK_SET);
    sf_complexread(data,sl->n12,sl->file);
}

void cslice_put(slice sl, int i3, sf_complex* data)
/*< put a slice at level i3 >*/
{
    sf_seek(sl->file,sl->start+i3*(sl->n12)*sizeof(sf_complex),SEEK_SET);
    sf_complexwrite(data,sl->n12,sl->file);
}

/*------------------------------------------------------------*/

fslice fslice_init(int n1, int n2, size_t size)
/*< initialize a sliceable file object >*/
{
    fslice sl;

    sl = (fslice) sf_alloc(1,sizeof(*sl));
    sl->file = sf_tempfile(&(sl->name), "w+b");
    sl->n1 = n1*size;
    sl->n2 = n2;
    
    return sl;
}

void fslice_get(fslice sl, int i2, void* data)
/*< get a slice at level i2 >*/
{
    extern int fseeko(FILE *stream, off_t offset, int whence);

    fseeko(sl->file,i2*(sl->n1),SEEK_SET);
    if (!fread(data,sl->n1,1,sl->file))
        abort();
}

void fslice_put(fslice sl, int i2, void* data)
/*< put a slice at level i2 >*/
{
    extern int fseeko(FILE *stream, off_t offset, int whence);

    fseeko(sl->file,i2*(sl->n1),SEEK_SET);
    fwrite(data,sl->n1,1,sl->file);
}

void fslice_close(fslice sl)
/*< remove the file and free allocated storage >*/
{
    unlink(sl->name);
    free(sl);
}

void fslice_load(sf_file in, fslice sl, sf_datatype type)
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

void fslice_dump(sf_file out, fslice sl, sf_datatype type)
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

/*------------------------------------------------------------*/

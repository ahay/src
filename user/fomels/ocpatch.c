/* Out-of-core patching */
/*
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
/*^*/

#include <rsf.h>
/*^*/

#include "ocpatch.h"

static off_t n1, n2, *tab;
static int np;
static FILE *table;

void ocpatch_init(int dim     /* number of dimensions */, 
		  int nw      /* total patch size */, 
		  int np1     /* total number of patches */, 
		  int* npatch /* number of patches [dim] */, 
		  int* nwall  /* data size [dim] */, 
		  int* nwind  /* patch size [dim] */)
/*< initialize >*/
{
    int i, ip, i2, ff[SF_MAX_DIM], gg[SF_MAX_DIM], t2;
    off_t t;

    n1 = nwind[0];
    n2 = nw/n1;
    np = np1;

    table = tmpfile();
    if (NULL == table) 
	sf_error("%s: trouble with temporary storage",__FILE__);
    tab = sf_largeintalloc(n2);

    for (ip=0; ip < np; ip++) {
	t = ip;
	for(i = 0; i < dim; i++) {
	    t2 = t%npatch[i];
	    t /= npatch[i];

	    if(npatch[i] == 1) {
		ff[i] = 0;
	    } else if (t2 == npatch[i]-1) {
		ff[i] = nwall[i] - nwind[i];
	    } else {	    
		ff[i] = t2*(nwall[i] - nwind[i])/(npatch[i] - 1.0);
	    }
	}	
	t = sf_cart2line (dim-1, nwall+1, ff+1);
	tab[0] = (t*nwall[0] + ff[0])*sizeof(float);
	for (i2=1; i2 < n2; i2++) {
	    t = i2;
	    for (i = 1; i < dim-1; i++) {
		/* cartesian coordinates in window */
		gg[i] = ff[i] + t%nwind[i];
		t /= nwind[i];
	    }
	    t += ff[dim-1];
	    for (i = dim-2; i >= 1; i--) {
		/* line coordinates in input */
		t = t*nwall[i] + gg[i];
	    }
	    tab[i2] = (t*nwall[0] + ff[0])*sizeof(float); 
	}
	if (n2 != fwrite(tab,sizeof(off_t),n2,table)) 
	    sf_error("%s: trouble writing table",__FILE__);
    }
}

void ocpatch_close(void)
/*< free allocated storage >*/
{
    free(tab);
    fclose(table);
}

void ocpatch_lop (int ip      /* patch number */, 
		  bool adj    /* pull or push */, 
		  FILE *wall  /* data */, 
		  float* wind /* patch */)
/*< apply patching >*/
{
    off_t i2;
    extern int fseeko(FILE *stream, off_t offset, int whence);
    
    if (0 != fseeko(table,ip*n2*sizeof(off_t),SEEK_SET))
	sf_error("%s: table seeking error:",__FILE__);
    if (n2 != fread(tab,sizeof(off_t),n2,table))
	sf_error("%s: trouble reading table",__FILE__);

    for (i2=0; i2 < n2; i2++, wind += n1) {
	if (0 != fseeko(wall,tab[i2],SEEK_SET))
	    sf_error("%s: seeking error:",__FILE__);

	if (adj) {
	    if (n1 != fwrite(wind,sizeof(float),n1,wall))
		sf_error("%s: writing error:",__FILE__);
	} else {
	    if (n1 != fread(wind,sizeof(float),n1,wall))
		sf_error("%s: reading error:",__FILE__);
	}
    }
}

void ocpatch_flop (int ip       /* patch number */, 
		   bool adj     /* pull or push */, 
		   sf_file wall /* data */, 
		   float* wind  /* patch */)
/*< apply patching with RSF files >*/
{
    off_t i2;
    extern int fseeko(FILE *stream, off_t offset, int whence);

    if (0 != fseeko(table,ip*n2*sizeof(off_t),SEEK_SET))
	sf_error("%s: table seeking error:",__FILE__);
    if (n2 != fread(tab,sizeof(off_t),n2,table))
	sf_error("%s: trouble reading table",__FILE__);
    
    for (i2=0; i2 < n2; i2++, wind += n1) {
	sf_seek(wall,tab[i2],SEEK_SET);
	
	if (adj) {
	    sf_floatwrite(wind,n1,wall);
	} else {
	    sf_floatread(wind,n1,wall);
	}
    }
}

/* Additional operations with RSF files. */
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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <limits.h>

#include "files.h"
#include "getpar.h"
#include "simtab.h"
#include "error.h"
#include "alloc.h"
#include "_defs.h"

#include "_bool.h"
#include "file.h"
/*^*/

int sf_filedims (sf_file file, /*@out@*/ int *n) 
/*< Find file dimensions.
--- 
Outputs the number of dimensions dim and a dimension array n[dim] >*/
{
    int i, dim;
    char key[3];

    dim = 1;
    for (i=0; i < SF_MAX_DIM; i++) {
	(void) snprintf(key,3,"n%d",i+1);
	if (!sf_histint(file,key,n+i)) {
	    n[i]=1;
	} else if (n[i] > 1) {
	    dim=i+1;
	}
    }
    return dim;
}

int sf_largefiledims (sf_file file, /*@out@*/ off_t *n) 
/*< Find file dimensions.
--- 
Outputs the number of dimensions dim and a dimension array n[dim] >*/
{
    int i, dim;
    char key[3];

    dim = 1;
    for (i=0; i < SF_MAX_DIM; i++) {
	(void) snprintf(key,3,"n%d",i+1);
	if (!sf_histlargeint(file,key,n+i)) {
	    n[i]=1;
	} else if (n[i] > 1) {
	    dim=i+1;
	}
    }
    return dim;
}

int sf_memsize()
/*< Returns memory size by:
  1. checking RSFMEMSIZE environmental variable
  2. using hard-coded "def" constant
  >*/
{
    char *memenv;
    long longsize;
    int memsize;
    const int def=100; /* default value (Mbytes) */

    if (NULL != (memenv = getenv("RSFMEMSIZE"))) {
	longsize = strtol(memenv,NULL,10);
	if (ERANGE == errno || longsize < 0 || longsize > INT_MAX) 
	    sf_error("wrong value in RSFMEMSIZE environmental variable");
	memsize = longsize;
    } else {
	memsize = def;
    }
    return memsize;
}

off_t sf_filesize (sf_file file) 
/*< Find file size (product of all dimensions) >*/
{    
    return sf_leftsize (file, 0);
}

off_t sf_leftsize (sf_file file, int dim) 
/*< Find file size for dimensions greater than dim >*/
{
    off_t ni, size;
    char key[3];

    for (size=1; dim < SF_MAX_DIM; dim++, size *= ni) {
	(void) snprintf(key,3,"n%d",dim+1);
	if (!sf_histlargeint(file,key,&ni)) break;
    }
    return size;
}

void sf_cp(sf_file in, sf_file out)
/*< Copy file in to file out >*/
{
    int esize;
    off_t nsiz, nbuf;
    char *buf;
    
    nsiz = sf_bytes (in);
    if (nsiz < 0) { /* reading from "stdin" */
	nsiz = sf_filesize(in);
	esize = sf_esize(in);
	nsiz *= esize;
    }

    sf_fileflush(out,in);
    sf_setformat(in,"raw");
    sf_setformat(out,"raw");

    nbuf = sf_bufsiz(in);
    buf = sf_charalloc(SF_MIN(nbuf,nsiz));

    while (nsiz > 0) {
	if (nbuf > nsiz) nbuf=nsiz;
	sf_charread (buf,nbuf,in);
	sf_charwrite (buf,nbuf,out);
	nsiz -= nbuf;
    }

    free(buf);
}

void sf_rm(const char* filename, bool force, bool verb, bool inquire)
/*< Remove an RSF file.
---
force, verb, and inquire flags should behave similar to the corresponding flags in the Unix "rm" command. >*/
{
    int c, c2;
    char cc, *in;
    FILE *file=NULL, *query=NULL;
    sf_simtab tab;
    struct stat buf;
    mode_t mod;
    const int tabsize=10;

    tab = sf_simtab_init (tabsize);
    query = fopen ("/dev/tty","w+");
    if (inquire) {
	if (NULL == query) sf_error ("%s: Cannot open terminal",__FILE__);
	setbuf (query,NULL);
	fprintf (query,"sf_rm: Remove '%s'? ",filename);
	c2 = c = getc (query);
	while (c2 != EOF && (cc = (char) c2) != '\n' && cc != '\0') 
	    c2 = getc(query);
	cc = (char) c;
	if ('y' != cc && 'Y' != cc) return;
    }
    if (verb) sf_warning("sf_rm: Removing header %s",filename);
    file = fopen (filename,"r");
    if (NULL == file) sf_error ("%s: Cannot open file %s:",__FILE__,filename);
    sf_simtab_input (tab,file,NULL);
    (void) fclose (file);
    in = sf_simtab_getstring (tab,"in");
    if (NULL == in) {
	if (force) {
	    sf_warning ("%s:  File %s has no in=",__FILE__,filename);
	} else {
	    sf_error ("%s:  File %s has no in=",__FILE__,filename);
	}
    }
    if (0 != remove(filename)) 
	sf_error ("%s: Trouble removing header file %s:",__FILE__,filename);

    if (NULL != in && 0 != strcmp(in,"stdin")) {
	if (verb) sf_warning("sf_rm: Removing data %s",in);
	if (!force) {
	    if (0 != stat(in,&buf)) 
		sf_error ("%s: Trouble with file %s:",__FILE__,in);
	    /* (owner and can write) or (not-owner and others can write) */
	    mod = (buf.st_uid == getuid())? S_IWUSR:S_IWOTH;
	    if (0 == (buf.st_mode & mod)) {		    
		fprintf (query,"sf_rm: Remove protected file '%s'? ",in);
		c2 = c = getc (query);
		while (c2 != EOF && (cc = (char) c2) != '\n' && cc != '\0') 
		    c2 = getc(query);
		cc = (char) c;
		if ('y' != cc && 'Y' != cc) return;
	    }
	}
	if (0 != remove(in)) {
	    if (force) {
		sf_warning ("%s: Trouble removing data file %s:",__FILE__,in);
	    } else {
		sf_error ("%s: Trouble removing data file %s:",__FILE__,in);
	    }
	}
    }
    if (NULL != query)
        (void) fclose (query);
    sf_simtab_close (tab);
}

off_t sf_shiftdim(sf_file in, sf_file out, int axis) 
/*< shift grid after axis by one dimension forward >*/
{
    int j, ni;
    float f;
    off_t n3;
    char key1[7], key2[7], *val;

    n3 = 1;
    for (j=axis; j < SF_MAX_DIM; j++) {
	sprintf(key2,"n%d",j+1);
	sprintf(key1,"n%d",j);
	if (!sf_histint(in,key1,&ni)) {
	    sf_putint(out,key2,1);
	    break;
	}
	sf_putint(out,key2,ni);
	n3 *= ni;
	
	sprintf(key2,"o%d",j+1);
	sprintf(key1,"o%d",j);
	if (sf_histfloat(in,key1,&f)) sf_putfloat(out,key2,f);

	sprintf(key2,"d%d",j+1);
	sprintf(key1,"d%d",j);
	if (sf_histfloat(in,key1,&f)) sf_putfloat(out,key2,f);

	sprintf(key2,"label%d",j+1);
	sprintf(key1,"label%d",j);
	if (NULL != (val = sf_histstring(in,key1))) 
	    sf_putstring(out,key2,val);

	sprintf(key2,"unit%d",j+1);
	sprintf(key1,"unit%d",j);
	if (NULL != (val = sf_histstring(in,key1))) 
	    sf_putstring(out,key2,val);
    }

    return n3;
}

off_t sf_shiftdim2(sf_file in, sf_file out, int axis) 
/*< shift grid after axis by two dimension forward >*/
{
    int j, ni;
    float f;
    off_t n3;
    char key1[7], key2[7], *val;

    n3 = 1;
    for (j=axis; j < SF_MAX_DIM; j++) {
	sprintf(key2,"n%d",j+2);
	sprintf(key1,"n%d",j);
	if (!sf_histint(in,key1,&ni)) {
	     sf_putint(out,key2,1);
	     break;
	}
	sf_putint(out,key2,ni);
	n3 *= ni;
	
	sprintf(key2,"o%d",j+2);
	sprintf(key1,"o%d",j);
	if (sf_histfloat(in,key1,&f)) sf_putfloat(out,key2,f);

	sprintf(key2,"d%d",j+2);
	sprintf(key1,"d%d",j);
	if (sf_histfloat(in,key1,&f)) sf_putfloat(out,key2,f);

	sprintf(key2,"label%d",j+2);
	sprintf(key1,"label%d",j);
	if (NULL != (val = sf_histstring(in,key1))) 
	    sf_putstring(out,key2,val);

	sprintf(key2,"unit%d",j+2);
	sprintf(key1,"unit%d",j);
	if (NULL != (val = sf_histstring(in,key1))) 
	    sf_putstring(out,key2,val);
    }

    return n3;
}

off_t sf_shiftdimn(sf_file in, sf_file out, int axis, int n) 
/*< shift grid after axis by n dimension forward >*/
{
    int j, ni;
    float f;
    off_t n3;
    char key1[7], key2[7], *val;

    n3 = 1;
    for (j=axis; j < SF_MAX_DIM; j++) {
	if ((j+n) >= SF_MAX_DIM) sf_error ("Dimension shift is out of bounds");
	sprintf(key2,"n%d",j+n);
	sprintf(key1,"n%d",j);
	if (!sf_histint(in,key1,&ni)) {
	     sf_putint(out,key2,1);
	     break;
	}
	sf_putint(out,key2,ni);
	n3 *= ni;
	
	sprintf(key2,"o%d",j+n);
	sprintf(key1,"o%d",j);
	if (sf_histfloat(in,key1,&f)) sf_putfloat(out,key2,f);

	sprintf(key2,"d%d",j+n);
	sprintf(key1,"d%d",j);
	if (sf_histfloat(in,key1,&f)) sf_putfloat(out,key2,f);

	sprintf(key2,"label%d",j+n);
	sprintf(key1,"label%d",j);
	if (NULL != (val = sf_histstring(in,key1))) 
	    sf_putstring(out,key2,val);

	sprintf(key2,"unit%d",j+n);
	sprintf(key1,"unit%d",j);
	if (NULL != (val = sf_histstring(in,key1))) 
	    sf_putstring(out,key2,val);
    }

    return n3;
}

off_t sf_unshiftdim(sf_file in, sf_file out, int axis) 
/*< shift grid after axis by one dimension backward >*/
{
    int j, ni;
    off_t n3;
    float f;
    char key1[7], key2[7], *val;

    n3 = 1;
    for (j=axis; j < SF_MAX_DIM; j++) {
	sprintf(key2,"n%d",j);
	sprintf(key1,"n%d",j+1);
	if (!sf_histint(in,key1,&ni)) {
	    sf_putint(out,key2,1);
	    break;
	}
	sf_putint(out,key2,ni);
	n3 *= ni;
	
	sprintf(key2,"o%d",j);
	sprintf(key1,"o%d",j+1);
	if (sf_histfloat(in,key1,&f)) sf_putfloat(out,key2,f);

	sprintf(key2,"d%d",j);
	sprintf(key1,"d%d",j+1);
	if (sf_histfloat(in,key1,&f)) sf_putfloat(out,key2,f);

	sprintf(key2,"label%d",j);
	sprintf(key1,"label%d",j+1);
	if (NULL != (val = sf_histstring(in,key1))) 
	    sf_putstring(out,key2,val);

	sprintf(key2,"unit%d",j);
	sprintf(key1,"unit%d",j+1);
	if (NULL != (val = sf_histstring(in,key1))) 
	    sf_putstring(out,key2,val);
    }

    return n3;
}

off_t sf_unshiftdim2(sf_file in, sf_file out, int axis) 
/*< shift grid after axis by two dimension backward >*/
{
    int j, ni;
    off_t n3;
    float f;
    char key1[7], key2[7], *val;

    n3 = 1;
    for (j=axis; j < SF_MAX_DIM; j++) {
	sprintf(key2,"n%d",j);
	sprintf(key1,"n%d",j+2);
	if (!sf_histint(in,key1,&ni)) {
	    sf_putint(out,key2,1);
	    sprintf(key2,"n%d",j+1);
	    sf_putint(out,key2,1);
	    break;
	}
	sf_putint(out,key2,ni);
	n3 *= ni;
	
	sprintf(key2,"o%d",j);
	sprintf(key1,"o%d",j+2);
	if (sf_histfloat(in,key1,&f)) sf_putfloat(out,key2,f);

	sprintf(key2,"d%d",j);
	sprintf(key1,"d%d",j+2);
	if (sf_histfloat(in,key1,&f)) sf_putfloat(out,key2,f);

	sprintf(key2,"label%d",j);
	sprintf(key1,"label%d",j+2);
	if (NULL != (val = sf_histstring(in,key1))) 
	    sf_putstring(out,key2,val);

	sprintf(key2,"unit%d",j);
	sprintf(key1,"unit%d",j+2);
	if (NULL != (val = sf_histstring(in,key1))) 
	    sf_putstring(out,key2,val);
    }

    return n3;
}

bool sf_endian (void)
/*< Endianness test, returns true for little-endian machines >*/
{
    bool little_endian;

    union {
	unsigned char c[4];
	int i;
    } test;

    test.i=0;
    test.c[0] = (unsigned char) 1;
    
    assert (2 == sizeof(short) && 4 == sizeof(int)); /* fix this later */
    little_endian = (bool) (0 != (test.i << 8));
    
    return little_endian;
}

/* 	$Id: files.c 8223 2012-03-15 04:53:11Z vovizmus $	 */

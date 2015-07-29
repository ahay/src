/* Sort a dataset according to a header key. */
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
#include <stdlib.h>

#include <rsf.h>

struct skey {
    float key;
    off_t pos;
};

static int key_compare (const void *k1, const void *k2)
{
    float f1 = ((struct skey*) k1)->key;
    float f2 = ((struct skey*) k2)->key;
    return (f1 < f2)? -1: (f1 > f2)? 1: 0;
}

int main(int argc, char* argv[])
{
    int n1, n2, i2, esize;
    off_t pos;
    struct skey *sorted;
    float *unsorted;
    char *trace, *header;
    sf_file in, head, out;

    sf_init (argc,argv);
    in = sf_input ("in");
    out = sf_output ("out");
 
    header = sf_getstring("head");
    /* header file */
    if (NULL == header) { 
	header = sf_histstring(in,"head");
	if (NULL == header) sf_error("Need head=");
    }

    head = sf_input(header);
    if (SF_FLOAT != sf_gettype(head))
	sf_error("Need float header");
    n2 = sf_filesize(head);
 
    unsorted = sf_floatalloc(n2);
    sorted = (struct skey*) sf_alloc(n2,sizeof(struct skey));
    
    sf_floatread(unsorted,n2,head);
    for (i2 = 0; i2 < n2; i2++) {
	sorted[i2].key = unsorted[i2];
	sorted[i2].pos = i2;
    }
    free (unsorted);
    sf_fileclose(head);

    qsort(sorted,n2,sizeof(struct skey),key_compare);
 
    if (!sf_histint(in,"n1",&n1)) n1=1;
    esize = sf_esize(in);
    n1 *= esize;

    trace = sf_charalloc(n1);

    sf_unpipe(in,((off_t) n1)*((off_t) n2));
    sf_fileflush(out,in);
    sf_setform(in,SF_NATIVE);
    sf_setform(out,SF_NATIVE);

    pos = sf_tell(in);
    for (i2=0; i2<n2; i2++) {
	sf_seek(in,pos+(sorted[i2].pos)*n1,SEEK_SET);
	sf_charread(trace,n1,in);
	sf_charwrite(trace,n1,out);
    }


    exit(0);
}
    
/* 	$Id: headersort.c 13857 2015-02-22 19:18:06Z sfomel $	 */

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
    int ikey;
    float fkey;
    off_t pos;
};

static int int_key_compare (const void *k1, const void *k2)
{
    int f1 = ((struct skey*) k1)->ikey;
    int f2 = ((struct skey*) k2)->ikey;
    return (f1 < f2)? -1: (f1 > f2)? 1: 0;
}


static int float_key_compare (const void *k1, const void *k2)
{
    float f1 = ((struct skey*) k1)->fkey;
    float f2 = ((struct skey*) k2)->fkey;
    return (f1 < f2)? -1: (f1 > f2)? 1: 0;
}

int main(int argc, char* argv[])
{
    int n1, n2, i2, esize;
    off_t pos;
    sf_datatype type;
    struct skey *sorted;
    int *input=NULL;
    float *finput=NULL;
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
    n2 = sf_filesize(head);
    type = sf_gettype(head);
    
    if (SF_FLOAT == type) {
	finput = sf_floatalloc(n2);
    } else if (SF_INT == type) {
	input = sf_intalloc(n2);
    } else {	
	sf_error("Need int or float header");
    }

    sorted = (struct skey*) sf_alloc(n2,sizeof(struct skey));

    if (SF_FLOAT == type) {
	sf_floatread(finput,n2,head);
	for (i2 = 0; i2 < n2; i2++) {
	    sorted[i2].fkey = finput[i2];
	    sorted[i2].pos = i2;
	}
	free (finput);
    } else {
	sf_intread(input,n2,head);
	for (i2 = 0; i2 < n2; i2++) {
	    sorted[i2].ikey = input[i2];
	    sorted[i2].pos = i2;
	}
	free (input);
    }
    sf_fileclose(head);
	
    qsort(sorted,n2,sizeof(struct skey),
	  (SF_FLOAT==type)? float_key_compare: int_key_compare);

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


/* Window a dataset based on a header mask.

The input data is a collection of traces n1xn2,
mask is an integer array os size n2, windowed is n1xm2,
where m2 is the number of nonzero elements in mask.
*/
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

#include <rsf.h>

int main(int argc, char* argv[])
{
    int n1, n2, j2, i2, j, esize, *mask;
    off_t pos;
    char *trace, key[3];
    sf_file in, head, out;

    sf_init (argc,argv);
 
    head = sf_input("mask");
    if (SF_INT != sf_gettype(head))
	sf_error("Need integer mask");
    n2 = sf_filesize(head);
 
    mask = sf_intalloc(n2);
    
    sf_intread(mask,n2,head);
    sf_fileclose(head);

    in = sf_input ("in");
    out = sf_output ("out");
 
    if (!sf_histint(in,"n1",&n1)) n1=1;

    j2 = sf_leftsize(in,1);
    if (j2 != n2) sf_error("Wrong input dimensions, need %d by %d",n1,n2);

    esize = sf_esize(in);
    n1 *= esize;
    trace = sf_charalloc(n1);
    
    for (j2=i2=0; i2 < n2; i2++) {
	if (mask[i2]) j2++;
    }
    sf_putint(out,"n2",j2);
    for (j=2; j < SF_MAX_DIM; j++) {
	(void) snprintf(key,3,"n%d",j+1);
	if (!sf_histint(in,key,&j2)) break;
	sf_putint(out,key,1);
    }

    sf_unpipe(in,(off_t) n1*n2);
    sf_fileflush(out,in);
    sf_setform(in,SF_NATIVE);
    sf_setform(out,SF_NATIVE);

    pos = sf_tell(in);
    for (i2=0; i2<n2; i2++) {
	if (mask[i2]) {
	    sf_seek(in,pos+ (off_t) i2*n1,SEEK_SET);
	    sf_charread(trace,n1,in);
	    sf_charwrite(trace,n1,out);
	}
    }


    exit(0);
}
    
/* 	$Id: headerwindow.c 7107 2011-04-10 02:04:14Z ivlad $	 */

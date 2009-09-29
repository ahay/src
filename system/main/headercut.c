/* Zero a portion of a dataset based on a header mask.

The input data is a collection of traces n1xn2,
mask is an integer array of size n2.
*/
/*
  Copyright (C) 2005 University of Texas at Austin
  Copyright (C) 2005 University of British Columbia
  
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
/* Modified from Gilles Hennenfent's applymask */

#include <stdio.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int n1, n2, j2, i2, esize, *mask;
    off_t pos;
    char *trace, *zero;
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
    zero = sf_charalloc(n1);
    memset(zero,0,n1);

    sf_unpipe(in,n1*n2);
    sf_fileflush(out,in);
    sf_setform(in,SF_NATIVE);
    sf_setform(out,SF_NATIVE);
   
    pos = sf_tell(in);
    for (i2=0; i2<n2; i2++) {
	if (mask[i2]) {
	    sf_seek(in,pos+i2*n1,SEEK_SET);
	    sf_charread(trace,n1,in);
	    sf_charwrite(trace,n1,out);
	} else {
	    sf_charwrite(zero,n1,out);
	}
    }


    exit(0);
}
    
/* 	$Id: headerwindow.c 1303 2005-08-17 02:08:33Z fomels $	 */

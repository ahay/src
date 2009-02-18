/* Denoising using stationary wavelet transform. */
/*
  Copyright (C) 2006 University of Texas at Austin
   
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
#include <rsf.h>

#include "swt.h"

int main(int argc, char* argv[])
{   
    int i,k,len_filter,n_layer,j;

    float *s,*ainverse;
    float *ch,*cl,*ch2,*cl2;
    float *s_orgi;
    
    int n_trace,n_number;
    float ratio;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    /* get number of samples in a trace */
    if (!sf_histint(in,"n1",&n_number)) sf_error("No n1= in input");
    /* get number of traces */
    n_trace = sf_leftsize(in,1); 

    if (!sf_getfloat("ratio",&ratio)) ratio=1.;
    /* ratio for denoising */

    if (!sf_getint("len_filter",&len_filter)) len_filter=2;
    /* filter length */

    if (!sf_getint("n_layer",&n_layer)) n_layer=2;
    /* number of wavelet transform layers */

    ch=sf_floatalloc(n_layer*n_number);
    cl=sf_floatalloc(n_number);
    ch2=sf_floatalloc(n_layer*n_number);
    cl2=sf_floatalloc(n_number);
    
    s=sf_floatalloc(n_number);
    s_orgi=sf_floatalloc(n_number);
    ainverse=sf_floatalloc(n_number);
    
    /* loop over traces */
    for(k=0;k<n_trace;k++)
    {
        sf_floatread(s,n_number,in);        /* read data from file */
		
	for(j=0;j<n_number;j++)
	    s_orgi[j]=s[j];

        multi_fwavelet(s,n_number,n_layer,cl,ch);

	multi_fwavelet(ch,n_number,n_layer,cl2,ch2);

	for(i=0;i<n_number;i++)
	    ch2[i]=ratio*ch2[i];

	multi_iwavelet(cl2,ch2,n_number,n_layer,ch);
	multi_iwavelet(cl,ch,n_number,n_layer,ainverse);
	
	sf_floatwrite(ainverse,n_number,out);
    }

    exit(0);
}


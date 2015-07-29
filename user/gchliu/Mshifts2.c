/* Generate shifts for 2-D regularized autoregression in complex domain. From (x,y,f) to (x,y,s,f) */
/*
  Copyright (C) 2012 University of Texas at Austin
  
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

int main(int argc, char* argv[])
{
    int n1, n2, nsx, nsy, n12, isx,isy;
    off_t i1, i2, i3, n3;
    sf_complex **trace, **trace1;
    sf_file in, shifts;

    sf_init(argc,argv);
    in = sf_input("in");
    shifts = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    
    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");
    if (!sf_getint("ns1",&nsx)) sf_error("Need ns1=");
    /* number of shifts in first dim */

    if (!sf_getint("ns2",&nsy)) sf_error("Need ns2=");
    /* number of shifts in second dim */

    if (nsx >= n1) sf_error("nsx=%d is too large",nsx);
    if (nsy >= n2) sf_error("nsy=%d is too large",nsy);

    sf_putint(shifts,"n3",(2*nsx+1)*(2*nsy+1)-1);
    sf_shiftdim(in, shifts, 3);

    sf_fileflush(shifts,in);
    sf_setform(in,SF_NATIVE);
    sf_setform(shifts,SF_NATIVE);
    

    n12 = n1*n2;

    trace  = sf_complexalloc2(n1,n2);
    trace1 = sf_complexalloc2(n1,n2);

/*(x,y,f)----(x,y,s,f)*/

    for (i3=0; i3 < n3; i3++) { 
        sf_complexread(trace[0],n12,in);
        for (isy = -nsy; isy < nsy+1; isy++) { 
        for (isx = -nsx; isx < nsx+1; isx++) { 
        if (!(isx==0 && isy==0)) { 
                   for (i2=0; i2 < n2; i2++) { 
                   for (i1=0; i1 < n1; i1++) { 
                         trace1[i2][i1] = sf_cmplx(0,0);
                   } 
                   } 

                   for (i2=0; i2 < n2; i2++) { 
                   for (i1=0; i1 < n1; i1++) { 
                   if (i2-isy < n2 && i1-isx < n1 && i2-isy >= 0 && i1-isx >= 0) {
                         trace1[i2][i1] = trace[i2-isy][i1-isx];
                   } 
                   } 
                   } 
	          sf_complexwrite(trace1[0],n12,shifts);    
        }         
        } 
        } 
         
        
    } 
    
    exit(0);
}

/* Stack a dataset over the second dimensions by SNR weighted method. */

/*
  Copyright (C) 2009 China university of Petroleum, Beijing, China
                     and University of Texas at Austin
       
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
#include <math.h>



int main(int argc, char* argv[])
{
    int nt, it, iw, w, shift, n2, i;
    float *indata, *outdata, *otherdata, win1, win2, sum1, sum2, eps, dt, t0, mul;
    

    sf_file in, out, other; 
    
    /* initialization */
    sf_init(argc,argv);
    in = sf_input("in");
    other = sf_input("other");
    out = sf_output("out");
    
    
    /* get rsf_file parameters */
    if (!sf_histint(in,"n1",&nt)) sf_error("Need n1=");
    if (!sf_histfloat(in,"d1",&dt)) dt=1.;
    if (!sf_histfloat(in,"o1",&t0)) t0=0.;
    
    n2 = sf_leftsize(in,1);


    if (!sf_getint("w",&w)) w = 50;  
    /*size of window*/
    if (!sf_getfloat("eps",&eps)) eps=0; 
    /* stable parameter*/      
    
    sf_putint(out,"n1",nt);

    sf_putfloat(out,"d1",dt);
    sf_putfloat(out,"o1",t0);
                  
    indata = sf_floatalloc(nt);
    otherdata = sf_floatalloc(nt);
    outdata = sf_floatalloc(nt);

 
    for (i=0; i< n2; i++) {
         
         sf_floatread(indata,nt,in);
         sf_floatread(otherdata,nt,other);
         for (it=0; it < nt; it++) {
              sum1 = 0.;
              sum2 = 0.;
              mul  = 0.;
              shift = SF_MAX(0,SF_MIN(nt-w, it-w/2-1));
              for (iw=0; iw <w; iw++){
                        win1  = indata[iw+shift];
                        win2  = otherdata[iw+shift];
                        mul  += win1*win2;
                        sum1 += win1*win1;
                        sum2 += win2*win2; 
              }
              outdata[it] = mul/(sqrt(sum1*sum2)+eps);
         }/*it*/
         sf_floatwrite(outdata,nt,out);
    } /* n2*/

    exit(0);
}



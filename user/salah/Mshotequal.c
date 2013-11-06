/* sfshotequal projects amplitudes of each shot to Z-score distribution*/
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

int main (int argc, char* argv[])
{
    int n1,n2,i3,i2,n4,j,n12,fold,ndim, mdim;
    int n[SF_MAX_DIM], m[SF_MAX_DIM], *mk;
    float *d, *sht, *x, mean, std, sum, amp=0.0;
    bool verb;
    sf_file in, out, msk, scl;

    sf_init (argc,argv);
    in  = sf_input("in");
    out = sf_output("out");
    
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float");
    
    if (!sf_histint(in,"n1",&n1)) sf_error("Need n1=");
    if (!sf_histint(in,"n2",&n2)) sf_error("Need n2=");
    n4 = sf_leftsize(in,2);
    n12=n1*n2;

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    /* input file*/
    if (! (sf_getstring("mask") || 
               sf_histstring(in,"mask")))
         sf_error("Need mask=");

    msk = sf_input("mask");
    if (SF_INT != sf_gettype(msk)) sf_error("Need integer mask");

    sf_getfloat("amp",&amp);
    /*Exclude amplitudes greater than amp && less than -amp for statistics computations*/

    ndim = sf_filedims(in,n);
    mdim = sf_filedims(msk,m);

    if (mdim != ndim) 
	sf_error("Wrong dimensions: %d != %d",mdim,ndim);
    
    for (j=0; j < ndim; j++) {
	if (m[j] != n[j])
            sf_error("Size mismatch [n%d]: %d != %d",
		     j+1,m[j],n[j]);
    }

    /* output file*/
    if (NULL != sf_getstring("scaler")){
        scl = sf_output("scaler");
        sf_unshiftdim(in, scl, 2);
        sf_putint(scl,"n1",2);
    } else {
	scl = NULL;
    }
    
    d   = sf_floatalloc(n12*n4);
    x   = sf_floatalloc(2*n4);
    sht = sf_floatalloc(n12);
    mk  = sf_intalloc(n12);
    
    /* loop through time samples */
    for (i3=0; i3 < n4; i3++) {
        mean=0.0;
        fold=0;
        std =0.0;
        sum=0.0;
	sf_floatread(sht,n12,in);
        sf_intread(mk,n12,msk);
        /* compute mean */
	for (i2=0; i2 < n12; i2++) {
           if (sht[i2]!=0.0 && mk[i2]!=0 ) {
                if ( amp==0.0 || fabsf(sht[i2]) < fabsf(amp) ){        
                    sum  +=sht[i2];
                    fold +=1;
                }
           }     
        }
        if (fold > 0) mean=sum/fold;
        
        /* compute standard deviation */
        for (i2=0; i2 < n12; i2++) {
           if (sht[i2]!=0.0 && mk[i2]!=0 ) {
               //if (!(amp && fabsf(sht[i2]) < fabsf(amp)))
               //     continue;
               if ( amp==0.0 || fabsf(sht[i2]) < fabsf(amp) )
                  std += (sht[i2]-mean)*(sht[i2]-mean);
           }     
        }
        if (fold > 0) std =sqrtf(std/(fold-1));
        
        /* scale time samples*/
        for (i2=0; i2 < n12; i2++) 
            if (sht[i2]!=0.0 && mean!=0.0 && std!=0.0)
                d[i2+i3*n12]=(sht[i2]-mean)/std;
            else
                d[i2+i3*n12]=sht[i2];
        x[0+i3*2]=mean;
        x[1+i3*2]=std;

        if (verb) sf_warning("shot %8d-> mean:%8g std:%8g, fold:%8d\n;",i3+1,mean,std,fold);
    }
    sf_floatwrite (d,n4*n12,out);
    if (scl)
       sf_floatwrite (x,n4*2,scl);
    exit(0);
}

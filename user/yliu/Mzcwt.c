/*  Improve signal resolution using zero-crossing of wavelet transform. */
/*
  Copyright (C) 2008 Jilin University
  
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
#include <stdio.h>
#include <math.h>

int main (int argc, char* argv[]) 
{
    int n1, n2, n3, n12; /*n1 is trace length, n2 is the number of traces, n3 is the number of 3th axis*/
    int maxscale, outscale;
    int i, ii, jj, kk, mm, number, nn3, kc, num;
    int ncoef, temp;
    float *h,*g,*x,*result;
    float *a, *d;
    int adn1, adn2;
    sf_file in, out;
    
    sf_init (argc, argv); 
    in = sf_input("in");
    out = sf_output("out");
    
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);
    n12=n1*n2;
    /* get the trace length (n1) and the number of traces (n2) and n3*/
    
    if (!sf_getint("maxscale",&maxscale)) maxscale=5;
    /* The maximum decomposition scale (default=5)*/

    if (!sf_getint("outscale",&outscale)) outscale=1;
    /* The output scale (default=1)*/
 
    if (outscale > maxscale)  sf_error("output scales need to be less than maximum decomposition scales"); 

    ncoef=5;
    temp=ncoef*(pow(2,(maxscale-1)));
    adn2=maxscale+1;
    adn1=n1+2*temp+1;

    h = sf_floatalloc(2*ncoef+1);
    g = sf_floatalloc(2*ncoef+1);
    x = sf_floatalloc(n12);
    result = sf_floatalloc(n12);

    a = sf_floatalloc(adn1*adn2);
    d = sf_floatalloc(adn1*adn2);

    h[0]=0.0032;
    h[1]=-0.0132;
    h[2]=-0.0393;
    h[3]=0.0450;
    h[4]=0.2864;
    h[5]=0.4347;
    
    g[0]=0.0039;
    g[1]=0.0062;
    g[2]=-0.0226;
    g[3]=-0.1120;
    g[4]=-0.2309;
    g[5]=0.7118;

    for(i = 1; i <= 5; i++) {
	h[2*ncoef+1-i]=h[i-1];
	g[2*ncoef+1-i]=g[i-1];
    }
   
    for(nn3 = 0; nn3 < n3; nn3++){
	sf_warning("slice %d of %d;",nn3+1,n3);
	sf_floatread(x,n12,in);
	for(number = 0; number < n2; number++) {
	    for(ii = 1; ii<= maxscale; ii++) {
		for(jj = 0; jj < n1; jj++)	{
		    a[0*adn1+temp+jj]=x[number*n1+jj];
		} 
		kc=1;
		for(jj = 1; jj <= ii; jj++) {
		    num=ncoef*(pow(2,(jj-1)));
		    for(kk = 0; kk < num; kk++) {
			a[(jj-1)*adn1+n1+temp+kk]=a[(jj-1)*adn1+temp+n1-2-kk];
		    }
		    for(kk = -num; kk < 0; kk++) {
			a[(jj-1)*adn1+temp+kk]=a[(jj-1)*adn1+temp-kk+2];
		    }
		    for(kk = 0; kk < n1; kk++) {
			a[jj*adn1+temp+kk]=0.;
			d[jj*adn1+temp+kk]=0.;
		    }
		    for(kk = 0; kk < n1; kk++) {
			for(mm = -ncoef; mm <= ncoef; mm++) {
			    a[jj*adn1+temp+kk]=a[jj*adn1+temp+kk]+a[(jj-1)*adn1+temp+kk-kc*mm]*h[mm+ncoef];
			    d[jj*adn1+temp+kk]=d[jj*adn1+temp+kk]+a[(jj-1)*adn1+temp+kk-kc*mm]*g[mm+ncoef];
			}
		    }
		    kc=kc*2;
		}
		for( jj = 0; jj < n1; jj++) {
		    result[number*n1+jj]=d[outscale*adn1+temp+jj];
		}
	    }
	} 
	sf_floatwrite(result,n12,out);
    }
    sf_warning(".");
    exit (0);
}

/* 	$Id$	 */


/* Signal-to-noise ratio (decibel)
*/

/*
  Copyright (C) 2014 Xi'an Jiaotong University, UT Austin (Pengliang Yang)
   
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

int main(int argc, char *argv[])
{
    int i,n1,n2,n;
    float a, sum1, sum2;
    float *ss, *rr;
    sf_file in, ref, out;
	
    sf_init(argc,argv);
    in = sf_input("in");/* input signal */
    ref= sf_input("ref");/* reference signal */
    out= sf_output("out");/* computed SNR */

    sf_histint(in,"n1",&n1);
    n2=sf_leftsize(in,1);
    n=n1*n2;

    ss=sf_floatalloc(n);
    rr=sf_floatalloc(n);

    sf_floatread(ss,n,in);
    sf_floatread(rr,n,ref);
    sum1=sum2=0.;
    for(i=0; i<n; i++){
	sum1+=ss[i]*ss[i];
	a=ss[i]-rr[i];
	sum2+=a*a;
    }

    a=10.*log10f(sum1/sum2);
    sf_warning("SNR=%g dB",a);
    sf_floatwrite(&a, 1, out);
}

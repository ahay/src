/* Pick values at equal division points of 1-D array  

for example: There is an array a[5]={2.,1.,8.,5.,-4.}.
Value at bisection point is 2.0; Values at trisection is 1.0 and 5.0.
*/
/*
  Copyright (C) 2013 University of Texas at Austin
  
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

int main(int argc, char* argv[])
{
    bool abs;
    int i, n, m;
    float *dat,*val;
    sf_file in;
    
    sf_init(argc, argv);
    in=sf_input("in");
    
    n=sf_filesize(in);
    dat=sf_floatalloc(n);
    sf_floatread(dat, n, in);
    
    if(!sf_getint("divnum", &m)) m=2;
    /* number of equal division points */
    if(!sf_getbool("abs", &abs)) abs=false;
    /* sorting based on absolute value */
    
    if(abs) for(i=0; i<n; i++) dat[i]=fabsf(dat[i]);
    
    val=sf_floatalloc(m-1);
    for(i=1; i<=m-1; i++){
       val[i-1]=sf_quantile((int)(n/m*i),n,dat);
       printf("value at %dth of %d points: %f\n", i, m-1, val[i-1]);
    }
    
    exit(0);
}
    

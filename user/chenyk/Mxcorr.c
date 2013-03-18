/* Cross-correlation function  
  C=XCORR(X,Y,L), computes the (auto/cross) correlation over the range of lags:
  -L to L, i.e., 2*L+1 lags. If L is left out, default is L=n1-1, 
  where n1 is the length of Y.
*/

/* 
  Copyright (C) 2013 the University of Texas at Austin 

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
  int  i, ix, n1x, n1y, maxlag;
  float  *x, *y, *corr, *lag;
  sf_file in, match, out, lagfile;
    
  sf_init(argc,argv);
  in = sf_input("in");     /*input vector x*/
  match=sf_input("match"); /*input vector y*/
  out = sf_output("out");  /*output auto/cross correlation function*/
  if(!sf_histint(in,"n1",&n1x)) sf_error("No n1= in input");
  if(!sf_histint(match,"n1",&n1y)) sf_error("No n1= in match");
  if(!sf_getint("l",&maxlag)) maxlag=n1y-1;    /*maxlag of auto/cross correlation function*/
  
  sf_putfloat(out,"o1",-maxlag);
  sf_putfloat(out,"d1",1);
  sf_putint(out,"n1",2*maxlag+1);
  
  x = sf_floatalloc(n1x);
  y = sf_floatalloc(n1y);
  corr=sf_floatalloc(2*maxlag+1);
  sf_floatread(x,n1x,in);
  sf_floatread(y,n1y,match);
  
  for (i=0;i<2*maxlag+1;i++)
    {   corr[i]=0;	
      for (ix=0; ix < n1x; ix++)
	{	
	  if(ix-(i-maxlag) < 0 || ix -(i-maxlag) >n1y-1)
	    corr[i] += 0;
	  else 
	    corr[i] += x[ix]*y[ix-(i-maxlag)];   		
	}
    }
  
  sf_floatwrite(corr,2*maxlag+1,out);
  
  /* write to the lagfile */
  if(NULL!=sf_getstring("lagfile"))
    { lagfile=sf_output("lagfile"); /*output auto/cross correlation lag file*/ 
      sf_putfloat(lagfile,"o1",1);
      sf_putfloat(lagfile,"d1",1);
      sf_putint(lagfile,"n1",2*maxlag+1);
      lag=sf_floatalloc(2*maxlag+1);
      for(i=0;i<2*maxlag+1;i++)
	lag[i]=-maxlag+i;
      sf_floatwrite(lag,2*maxlag+1,lagfile);}	
  
  
  exit(0);
}

     

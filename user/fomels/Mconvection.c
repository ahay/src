/* Convection filter. */
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
#include <math.h>

#include <rsf.h>

int main (int argc, char* argv[])
{
  bool verb;
  int niter, n1, n2, n12, i1, i2, rect[2], m[2];
  float **data, **rhs, ***conv, ***lhs, dif, mean;
  sf_file in, flt;

  sf_init(argc, argv);
  in = sf_input("in");
  flt = sf_input("out");

  if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
  if (!sf_histint(in,"n1",&n1)) sf_error("Need n1= in input");
  if (!sf_histint(in,"n2",&n2)) sf_error("Need n2= in input");

  sf_putint(flt, "n3", 2);

  if (!sf_getint("rect1",&rect[0])) rect[0]=1;
  if (!sf_getint("rect2",&rect[1])) rect[1]=1;
  /* smoothing */

  if (!sf_getint("niter",&niter)) niter=100;
  /* number of iterations */

  if (!sf_getbool("verb",&verb)) verb=true;
  /* verbosity flag */
  
  data = sf_floatalloc2(n1,n2);
  rhs = sf_floatalloc2(n1,n2);
  
  conv = sf_floatalloc3(n1,n2,2);
  lhs = sf_floatalloc3(n1,n2,2);

  n12 = n1*n2;
  m[0] = n1;
  m[1] = n2;

  sf_multidivn_init(2, 2, n12, m, rect, **lhs, NULL, verb); 

  sf_floatread(data[0],n12,in);

  mean = 0.;
  for (i2=0; i2 < n2-1; i2++) {
    rhs[i2][0] = 0.0f;
    lhs[0][i2][0] = 0.0f;
    lhs[1][i2][0] = 0.0f;
    for (i1=1; i1 < n1-1; i1++) {
      dif = data[i2+1][i1] - data[i2][i1];
      rhs[i2][i1] = dif;
      lhs[0][i2][i1] = dif - data[i2+1][i1+1] + data[i2][i1-1];
      lhs[1][i2][i1] = dif - data[i2+1][i1-1] + data[i2][i1+1];
      mean += lhs[0][i2][i1]*lhs[0][i2][i1] + lhs[1][i2][i1]*lhs[1][i2][i1];
    }
    rhs[i2][n1-1] = 0.0f;
    lhs[0][i2][n1-1] = 0.0f;
    lhs[1][i2][n1-1] = 0.0f;
  }
  for (i1=0; i1 < n1; i1++) {
    rhs[n2-1][i1] = 0.0f;
    lhs[0][n2-1][i1] = 0.0f;
    lhs[1][n2-1][i1] = 0.0f;
  }
  mean = sqrtf (mean/(n12*2));
  for (i2=0; i2 < n2-1; i2++) {
    for (i1=1; i1 < n1-1; i1++) {
      rhs[i2][i1] /= mean;
      lhs[0][i2][i1] /= mean;
      lhs[1][i2][i1] /= mean;
    }
  }

  sf_multidivn (*rhs,**conv,niter);
  sf_floatwrite(**conv,n12*2,flt);
  
  exit(0);
}
  

  

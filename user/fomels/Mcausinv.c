/* Smooth derivative by regularized causint inversion */
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
#include <rsf.h>

int main(int argc, char* argv[])
{
  int n1, niter;
  float eps;
  bool verb;
  float *x, *y, **xmov;
  sf_file in, out, mov;

  sf_init(argc,argv);
  in = sf_input("in");
  out = sf_output("out");

  mov = sf_output("mov");

  if (!sf_histint(in, "n1", &n1)) sf_error("No n1= in input");
  if (!sf_getint("niter", &niter)) niter=1;
  if (!sf_getfloat("eps", &eps)) eps=0.0;
  if (!sf_getbool("verb", &verb)) verb=false;

  x = sf_floatalloc(n1);
  y = sf_floatalloc(n1);
  xmov = sf_floatalloc2(n1,niter);

  sf_floatread(y, n1, in);
  
  sf_solver_reg (sf_causint_lop, sf_cgstep, sf_igrad1_lop, n1, n1, n1, x, y, niter, eps, "verb", verb, "xmov", xmov, "end");

  sf_floatwrite(x, n1, out);
  sf_floatwrite(xmov[0], n1*niter, mov);

  exit(0);
}

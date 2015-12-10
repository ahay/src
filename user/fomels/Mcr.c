/* Column-row matrix decomposition */
 /*
  Copyright (C) 2015 University of Texas at Austin
  
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

#include "pcg.h"
#include "cr.h"

int main(int argc, char* argv[])
{
    bool prec;
    int nr, nc, nx, niter;
    float *x, tol;
    sf_file row_in, col_in, row_out, col_out;

    sf_init(argc,argv);
    row_in = sf_input("in");
    col_in = sf_input("col_in");

    row_out = sf_output("out");
    col_out = sf_output("col_out");

    if (!sf_histint(row_in,"n1",&nr)) sf_error("No n1= in input");
    if (!sf_histint(col_in,"n1",&nc)) sf_error("No n1= in input");

    sf_putint(col_out,"n1",nc);

    nx = nr+nc;
    x = sf_floatalloc(nx);

    /* read input - B'd */
    sf_floatread(x,nr,row_in);
    sf_floatread(x+nr,nc,col_in);

    if (!sf_getint("niter",&niter)) niter=10; /* number of iterations */
    if (!sf_getfloat("tol",&tol)) tol=0.0f;   /* CG tolerance */
    if (!sf_getbool("prec",&prec)) prec=true; /* If apply preconditioning */

    /* Run PCG */
    cr_init(nr,nc);
    conjgrad(cr_apply,prec? cr_prec:NULL,nx,x,x,x,niter,tol);
    
    /* write output */
    sf_floatwrite(x,nr,row_out);
    sf_floatwrite(x+nr,nc,col_out);

    exit(0);
}
    
    

/* Signal component separation using plane-wave shaping. 

The program works with 2-D data.
*/
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
    
#include "copyk.h"
#include "trislk.h"

int main(int argc, char* argv[])
{
    int i, n1, n2, n12, n3, n123, rect1, rect2, niter;
    float eps, *d, *s, *p, ***pp;
    bool verb;
    sf_file in, out, dips;

    sf_init (argc,argv);
    in = sf_input("in");
    dips = sf_input("dips");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");

    if (!sf_histint(dips,"n1",&i) || i != n1) sf_error("Wrong n1= in dips");
    if (!sf_histint(dips,"n2",&i) || i != n2) sf_error("Wrong n2= in dips");
    if (!sf_histint(dips,"n3",&n3)) sf_error("No n3= in dips");

    n12 = n1*n2;
    n123 = n12*n3;
    sf_putint (out,"n3",n3);

    if (!sf_getint ("niter",&niter)) niter=50;
    /* maximum number of iterations */

    if (!sf_getfloat ("eps",&eps)) eps=1.;
    /* regularization parameter */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    if (!sf_getint("rect1",&rect1)) rect1=1;
    /* vertical smoothing radius */
    if (!sf_getint("rect2",&rect2)) rect2=1;
    /* lateral smoothing radius */

    s = sf_floatalloc(n123);
    p = sf_floatalloc(n123);
    d = sf_floatalloc(n12);
    pp = sf_floatalloc3(n1,n2,n3);

    sf_floatread (d,n12,in);
    sf_floatread (pp[0][0],n123,dips);

    trislk_init(n3,n1,n2,rect1,rect2,pp);
    copyk_init(n3,n12);

    sf_conjgrad_init(n123,n123,n12,n12,eps,1.e-6,verb,false);
    sf_conjgrad(NULL,copyk_lop,trislk_lop,p,s,d,niter);

    sf_floatwrite(s,n123,out);

    exit(0);
}

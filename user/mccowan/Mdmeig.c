/* Find eigenvalues and eigenvectors of an spd matrix. */
/*
  Copyright (C) 2007 The University of Texas at Austin

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

#include "dmeig.h"

int main(int argc, char* argv[])
{
    int n, n2, i3, n3;
    float *a, *eval, *evec;
    sf_file mat=NULL, evals=NULL, evecs=NULL;

    sf_init(argc,argv);
    mat = sf_input("in");
    evals = sf_output("eval");
    evecs = sf_output("out");

    if (SF_FLOAT != sf_gettype(mat)) sf_error("Need float input");
    if (!sf_histint(mat,"n1",&n)) sf_error("No n1= in input");
    if (!sf_histint(mat,"n2",&n2) || n2 != n) sf_error("Need n1=n2 in input");
    n3 = sf_leftsize(mat,2);
    n2 = n*n;

    sf_putint(evals,"n2",1);

    a = sf_floatalloc(n2);
    eval = sf_floatalloc(n);
    evec = sf_floatalloc(n2);

    for (i3=0; i3 < n3; i3++) {
	sf_floatread(a,n2,mat);
	
	dmeig(n,a,eval,evec);

	sf_floatwrite(eval,n, evals);
	sf_floatwrite(evec,n2,evecs);
    }


    exit(0);
}

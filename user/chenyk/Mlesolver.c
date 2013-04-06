/* Linear equations solver using Gauss Elimination */
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
#include "solver.h"
int main(int argc, char *argv[])
{
    int i3, i, j, n1, n2, n3; /*n1 is trace length, n2 is the number of traces, n3 is the number of 3th axis*/
    sf_file in, d, out;
    float **a, *rhs, *m;
    sf_init(argc,argv);

    in  = sf_input("in");
    d   = sf_input("d");
    out = sf_output("out");

    if(!sf_histint(in,"n1",&n1))sf_error("No n1 in input!");
    if(!sf_histint(in,"n2",&n2))sf_error("No n2 in input!");
    if(n1!=n2)sf_error("Not square matrix!");
    n3=sf_leftsize(in,2);

    sf_putint(out,"n1",n1);
    sf_putint(out,"n2",n3);

    gaussel_init(n1);
    rhs = sf_floatalloc(n1);
    m   = sf_floatalloc(n1);
    a   = sf_floatalloc2(n1,n1);

    for(i3=0;i3<n3;i3++)
	{
		sf_floatread(a[0],n1*n1,in);
		for(i=0;i<n1;i++)
		for(j=0;j<n1;j++)
		sf_warning("a[%d][%d]=%f",i+1,j+1,a[i][j]);

		sf_floatread(rhs,n1,d);
		gaussel_solve(a,rhs,m);
   		sf_floatwrite(m,n1,out);

		for(j=0;j<n1;j++)
		sf_warning("m[%d]=%f",j+1,m[j]);	
	}
    gaussel_close();
    exit(0);
}


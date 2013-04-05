/* Matrix algebra operation
Implement C=Oper(A,B).
Specially, when "Oper" stands for multiplication, C=AB,  where: 
C is a m*k matrix, m corresponds to dimension 2 while k corrsponds to dimension 1.
A is a m*n matrix, m corresponds to dimension 2 while n corrsponds to dimension 1.
B is a n*k matrix, n corresponds to dimension 2 while k corrsponds to dimension 1.
Like matlab matrix operation.*/
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
#include "matrix.h"

int main(int argc, char *argv[])
{
    int i, m, n, k, n1a, n1b, n2a, n2b, n3a, n3b; /*n1 is trace length, n2 is the number of traces, n3 is the number of 3th axis*/
    sf_file in, mat, out;
    float *a, *b, *c;
    char *type;
    sf_init(argc,argv);

    in =  sf_input("in");    
    /*  in is the input matrix A, with m rows and n columns */
    mat= sf_input("mat");
    /* mat is the input matrix B, with n rows and k columns */
    out = sf_output("out");
    /* out is the output matrix C, with m rows and k columns */

    if (NULL == (type=sf_getstring("type"))) type="mul";
    /* [mul, add, sub, dotmul] operation type, the default is mul  */


    if (!sf_histint(in,"n1",&n1a)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2a)) sf_error("No n2= in input");
    n3a = sf_leftsize(in,2);

    if (!sf_histint(mat,"n1",&n1b)) sf_error("No n1= in input");
    if (!sf_histint(mat,"n2",&n2b)) sf_error("No n2= in input");
    n3b = sf_leftsize(mat,2);
    
    if (type[0]=='m')
	{
	if (n1a!=n2b) sf_error("n1 of A must be equal to n2 of B for multiplication");
	m=n2a;n=n1a;k=n1b;
	}

    if (type[0]!='m')
	{
	if (n1a!=n1b || n2a!=n2b) sf_error(" n1 of A must be equal to n1 of B and n2 of A must be equal to n2 of B");
	m=n2a;n=n1a;k=n;	
	}
    if (n3a != n3b) sf_error("n3 must be equal");
    
    sf_putint(out, "n1", k);   
    sf_putint(out, "n2", m);
    sf_putint(out, "n3", n3a);

    a = sf_floatalloc(m*n);
    b = sf_floatalloc(n*k);
    c = sf_floatalloc(m*k);

    for(i=0;i<n3a;i++)  {
	sf_floatread(a, m*n, in);

	if(type[0]!='m')
	sf_floatread(b, m*n, mat);
	else
	sf_floatread(b,n*k,mat);

	switch(type[0])
	{
		case 'm': mul(a, b, m, n, k, c); break;
		case 'a': add(a, b, m, n, c); break;
		case 's': sub(a, b, m, n, c); break;
		case 'd': dotmul(a, b, m, n, c);
	}
	sf_floatwrite(c, m*k, out);
    }

    exit(0);
}


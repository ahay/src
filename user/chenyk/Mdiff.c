/* Compare the difference of two rsf data sets with the same size. */
/*
  Copyright (C) 2012 University of Texas at Austin
  
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
#include <stdio.h>
#include <math.h>

#include <rsf.h>


int main(int argc, char* argv[])
{
    int  i1, i2, n1,n2,n3,n12;
    float  s=0, *pp1, *pp2;
    sf_file inp1, inp2,  dif;

    sf_init(argc,argv);
    inp1 = sf_input("in");
    inp2 = sf_input("match");
    dif = sf_output("out");
    
    sf_putfloat(dif,"o2",0);
    sf_putfloat(dif,"o1",0);
    sf_putint(dif,"n1",1);  
    sf_putint(dif,"n2",1);    
	
    if (!sf_histint(inp1,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(inp1,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(inp1,"n3",&n3)) n3=1;

    n2 = sf_leftsize(inp1,1); /* left dimensions after the first one */
    n12 = n1*n2;
    
    if (n3!=1) {sf_putint(dif,"n3",1); sf_putint(dif,"d3",1); }
    
    pp1 = sf_floatalloc(n12);
    pp2 = sf_floatalloc(n12);
    sf_floatread(pp1,n12,inp1);
    sf_floatread(pp2,n12,inp2);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    s = s + pow((pp1[i2*n1+i1]-pp2[i2*n1+i1]),2);
	}
    }
    
    sf_warning("The difference is %f", s );
    sf_floatwrite(&s,1,dif);
    exit(0);
}




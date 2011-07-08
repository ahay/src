/* 1-D trapezoidal integration (2-D dataset) */
/*
  Copyright (C) 2011 KAUST
  
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
    sf_file Fi,Fo;
    sf_axis a1,a2;
    int n1,n2,i1,i2,i;
    float d1,*f,*g,s;

    sf_init(argc,argv);
    Fi = sf_input("in");
    Fo = sf_output("out");

    a1 = sf_iaxa(Fi,1);
    a2 = sf_iaxa(Fi,2);
    n1 = sf_n(a1);
    n2 = sf_n(a2);
    d1 = sf_d(a1);

    sf_oaxa(Fo,a1,1);
    sf_oaxa(Fo,a2,2);

    f = sf_floatalloc(n1*n2);
    g = sf_floatalloc(n1*n2);

    sf_floatread(f,n1*n2,Fi);

    for (i2=0; i2 < n2; i2++) {
	s = 0;
	g[i2*n1] = 0;
	for (i1=1; i1 < n1; i1++) {
	    i = i2*n1+i1;
	    s += 0.5f*d1*(f[i] + f[i-1]);
	    g[i] = s;
	}
    }
    
    sf_floatwrite(g,n1*n2,Fo);

    exit(0);
}

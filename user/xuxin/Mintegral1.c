/* integration */

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

void trape(/*@out@*/ float *g,/* g(t) */
	   const float *f,    /* f(t) = dg/dt */
	   int n,             /* nt   */
	   float d,           /* dt   */
	   float g0           /* g(0) */)
{
    int i = 1;
    float w = .5*d;

    g[0] = g0;
    while (i < n) {
	g[i] = g[i-1] + w*(f[i-1] + f[i]);
	i++;
    }
}

void simpson(/*@out@*/ float *g,/* g(t) */
	     const float *f,    /* f(t) = dg/dt */
	     int n,			    /* nt   */
	     float d,           /* dt   */
	     float g0           /* g(0) */)
{
    int i;
    float we,wo;

    we = 4./3 * d;
    wo = 1./3 * d;
    g[0] = g0;
    /* 3/8 rule for first interval */
    g[1] = g0 + (3./8*f[0]
		 + (9./8-1./3)*f[1]
		 + (9./8-4./3)*f[2]
		 + (3./8-1./3)*f[3])*d;
	
    i = 2;
    while (i < n) {
	g[i] = g[i-2] + we*f[i-1] + wo*(f[i] + f[i-2]);
	i++;
    }
}

int main(int argc, char* argv[])
{
    int n1,n2,i2;
    float d1,*f,*g;
    char *rule;
    sf_file Fin,Fout;

    sf_init(argc,argv);

    Fin = sf_input("in");
    if (SF_FLOAT != sf_gettype(Fin)) sf_error("Need type=float in input");
    if (!sf_histint(Fin,"n1",&n1))   sf_error("Need n1= in input");
    if (!sf_histfloat(Fin,"d1",&d1)) sf_error("Need d1= in input");

    n2 = sf_leftsize(Fin,1);
	
    rule = sf_getstring("rule");
    /* t, s : quadrature rules */

    Fout= sf_output("out");
    sf_putint(Fout,"n1",n1);
    sf_putint(Fout,"n2",n2);
    sf_putfloat(Fout,"d1",d1);

    f = sf_floatalloc(n1);
    g = sf_floatalloc(n1);
    for (i2=0; i2 < n2; i2++) {
	sf_floatread(f,n1,Fin);
	switch (rule[0]) {
	    case 's':
		simpson(g,f,n1,d1,0.); break;
	    case 't':
	    default:
		trape(g,f,n1,d1,0.); break;
	}
	sf_floatwrite(g,n1,Fout);
    }

    return 0;
}

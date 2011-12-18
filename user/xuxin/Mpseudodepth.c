/* pseudo-depth mapping  */
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
#include "spline.h"

int main(int argc,char *argv[])
{
	int n1,n2,n,i;
	float o,d,*z,*zz,*u,*uu;
	sf_file Fin,Fz,Fout;
	Spl *spl;

	sf_init(argc,argv);

	Fin = sf_input("in");    /* f(z)  */
	Fz  = sf_input("depth"); /* zz(z) */
	Fout= sf_output("out");	 /* f(zz) */

	n2 = sf_leftsize(Fin,1);
	if (!sf_histint(Fin,"n1",&n1) || n1 < 3  ||
		!sf_histint(Fz, "n1",&n ) || n != n1 ||
		sf_leftsize(Fz,1) != n2)
		sf_error("Need same n1=(n1>2) and n2= in input and x.");
	if (!sf_getint("n1",&n))   sf_error("Need n1="); /* new depth n */
	if (!sf_getfloat("o1",&o)) o = 0.;               /* new depth o */
	if (!sf_getfloat("d1",&d)) sf_error("Need d1="); /* new depth d */

	zz= sf_floatalloc(n);
	for (zz[0]=o, i=1; i < n; i++) zz[i] = zz[i-1] + d;
	
	sf_putint(Fout,"n1",n);	sf_putfloat(Fout,"o1",o);
	sf_putint(Fout,"n2",n2);sf_putfloat(Fout,"d1",d); 

	u = sf_floatalloc(n1);
	z = sf_floatalloc(n1);
	uu= sf_floatalloc(n);
	for (i=0; i < n2; i++) {
		sf_floatread(u,n1,Fin);
		sf_floatread(z,n1,Fz );
		spl = spline_init(u,z,n1);
		spline_eval(uu,zz,n,spl);
		spline_free(spl);
		sf_floatwrite(uu,n,Fout);
	}

	sf_fileclose(Fin); sf_fileclose(Fz); sf_fileclose(Fout);
	return 0;
}

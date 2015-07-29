/* depth to vertical-time interpolation
   z to tau : pseudodepth < FZ.rsf inv=n tau=tau.rsf n=ntau o=otau d=dtau > FT.rsf
   tau to z : pseudodepth < FT.rsf inv=y tau=tau.rsf > FZ.rsf */
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
    int nt,nz,nx,i;
    float ot,dt,oz,dz,*t,*T,*u,*U;
    bool inv,linear;
    sf_file Fin,Ftau,Fout;
    Spl *spl;

    sf_init(argc,argv);

    Fin  = sf_input ("in");
    Ftau = sf_input ("tau");
    Fout = sf_output("out");

    if (!sf_getbool("inv",&inv))       inv = false;
    /* if y, tau to z; if n, tau to z */
    if (!sf_getbool("linear",&linear)) linear = true;
    /* if y, linear spline; if n, cubic spline (buggy) */

    if (!sf_histint  (Ftau,"n1",&nz) || nz < 3 ||
	!sf_histfloat(Ftau,"o1",&oz) ||
	!sf_histfloat(Ftau,"d1",&dz))
	sf_error("Need n1=(n1>2) o1= d1= in tau");
    nx = sf_leftsize(Ftau,1);
    // sf_putint(Fout,"n2",nx);

    if (inv) {
	if (!sf_histint  (Fin,"n1",&nt) || nt < 3) sf_error("Need n1>2 in input");
	if (!sf_histfloat(Fin,"o1",&ot)) sf_error("Need o1= in input");
	if (!sf_histfloat(Fin,"d1",&dt)) sf_error("Need d1= in input");
	if (sf_leftsize(Fin,1) != nx) sf_error("Need n2=%d in input",nx);

    sf_putint  (Fout,"n1",nz); sf_putstring(Fout,"label1","z");
    sf_putfloat(Fout,"o1",oz); sf_putstring(Fout,"unit1","m");
    sf_putfloat(Fout,"d1",dz);
} else {
    if (!sf_histint(Fin,"n1",&i) || i != nz ||
	sf_leftsize(Fin,1) != nx)
	sf_error("Need n1=%d and n2=%d in input",nz,nx);

    if (!sf_getint  ("n",&nt)) sf_error("Need n="); /* tau n */
    if (!sf_getfloat("o",&ot)) ot = 0.;             /* tau o */
    if (!sf_getfloat("d",&dt)) sf_error("Need d="); /* tau d (>0) */

    sf_putint  (Fout,"n1",nt); sf_putstring(Fout,"label1","tau");
    sf_putfloat(Fout,"o1",ot); sf_putstring(Fout,"unit1","sec");
    sf_putfloat(Fout,"d1",dt);
}

u = sf_floatalloc(nz); /* u(z) */
t = sf_floatalloc(nz); /* t(z) */
T = sf_floatalloc(nt); /* T equally spaced t(z) */
U = sf_floatalloc(nt); /* u(T) */

for (T[0]=ot, i=1; i < nt; i++) T[i] = T[i-1] + dt;

for (i=0; i < nx; i++) {
    sf_floatread(t,nz,Ftau);

    if (inv) {
	sf_floatread(U,nt,Fin);
	spl = spline_init(U,T,nt,linear);
	spline_eval(u,t,nz,spl);
	spline_free(spl);
	sf_floatwrite(u,nz,Fout);
    } else {
	sf_floatread(u,nz,Fin);
	spl = spline_init(u,t,nz,linear);
	spline_eval(U,T,nt,spl);
	spline_free(spl);
	sf_floatwrite(U,nt,Fout);
    }
}

sf_fileclose(Fin);
sf_fileclose(Ftau);
sf_fileclose(Fout);
return 0;
}

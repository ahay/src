/* Compute a variogram of data values. */
/*
  Copyright (C) 2011 Jilin University

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
    int i, n, i2, n2, i1, n1, dh, nh, oh;
    float o1, d1, *fbuf, *vari;
    bool semi;
    sf_file in=NULL, out=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&n1)) sf_error("Need n1=");
    if (!sf_histfloat(in,"o1",&o1)) sf_error("Need o1=");
    if (!sf_histfloat(in,"d1",&d1)) sf_error("Need d1=");
    n2 = sf_leftsize(in,1);

    if (!sf_getint("dh",&dh)) dh=1;
    /* interval (number) of variogram lag */
    if (!sf_getint("oh",&oh)) oh=0;
    /* origin (number) of variogram lag */
    if (oh < 0 || oh >= n1/dh) oh=0;
    if (!sf_getint("nh",&nh)) nh=n1/dh-oh;
    /* number of variogram lag */

    if (!sf_getbool("semi",&semi)) semi=true;
    /* if y, output semivariogram */

    sf_settype(out,SF_FLOAT);
    sf_putint(out,"n1",nh);
    sf_putfloat(out,"o1",oh*dh*d1);
    sf_putfloat(out,"d1",dh*d1);
    sf_putstring(out,"label1","Lag");

    fbuf = sf_floatalloc(n1);
    vari = sf_floatalloc(nh);

    for (i2 = 0; i2 < n2; i2++) {
	sf_floatread(fbuf,n1,in);
	for (i1=0; i1 < nh; i1++) {
	    vari[i1]=0.;
	}	

      	for (i1=0; i1 < nh; i1++) {
	    n = 0;
	    for (i=0; i < n1; i++) {
		if (i+(i1+oh)*dh < n1) {
		    vari[i1] += (fbuf[i+(i1+oh)*dh] - fbuf[i])
			*(fbuf[i+(i1+oh)*dh] - fbuf[i]);
		    n++;
		}
	    }
	    if (semi) {
		vari[i1] /= (n*2.+FLT_EPSILON);
	    } else {
		vari[i1] /= (n*1.+FLT_EPSILON);
	    }
	}
	sf_floatwrite(vari,nh,out);
    }
    
    exit(0);
}

/* 	$Id$	 */

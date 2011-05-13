/* Compute a horizontal variogram of data slice. */
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
    int i, j, n, i1, n1, i2, n2, i3, n3, dh1, nh1, dh2, nh2;
    float o1, d1, o2, d2, *fbuf, *vari;
    bool semi, verb;
    sf_file in=NULL, out=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&n1)) sf_error("Need n1=");
    if (!sf_histfloat(in,"o1",&o1)) sf_error("Need o1=");
    if (!sf_histfloat(in,"d1",&d1)) sf_error("Need d1=");
    if (!sf_histint(in,"n2",&n2)) sf_error("Need n2=");
    if (!sf_histfloat(in,"o2",&o2)) sf_error("Need o2=");
    if (!sf_histfloat(in,"d2",&d2)) sf_error("Need d2=");
    n3 = sf_leftsize(in,2);

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity */
    if (!sf_getint("dh1",&dh1)) dh1=1;
    /* interval (jump) of variogram lag in first axis */
    if (!sf_getint("nh1",&nh1)) nh1=n1/dh1;
    /* number of variogram lag in first axis */

    if (!sf_getint("dh2",&dh2)) dh2=1;
    /* interval (jump) of variogram lag in second axis */
    if (!sf_getint("nh2",&nh2)) nh2=n2/dh2;
    /* number of variogram lag in second axis */

    if (!sf_getbool("semi",&semi)) semi=true;
    /* if y, output semivariogram */

    sf_settype(out,SF_FLOAT);
    sf_putint(out,"n1",2*nh1-1);
    sf_putfloat(out,"o1",(-nh1+1)*dh1*d1);
    sf_putfloat(out,"d1",dh1*d1);
    sf_putstring(out,"label1","Inline Lag");

    sf_putint(out,"n2",2*nh2-1);
    sf_putfloat(out,"o2",(-nh2+1)*dh2*d2);
    sf_putfloat(out,"d2",dh2*d2);
    sf_putstring(out,"label2","Crossline Lag");

    fbuf = sf_floatalloc(n1*n2);
    vari = sf_floatalloc((2*nh1-1)*(2*nh2-1));

    for (i3 = 0; i3 < n3; i3++) {
	sf_floatread(fbuf,n1*n2,in);
	for (i1=0; i1 < (2*nh1-1)*(2*nh2-1); i1++) {
	    vari[i1]=0.;
	}	

      	for (i2=(-nh2+1); i2 < nh2; i2++) {
	    if (verb)
		sf_warning("Slice %d of %d;",i2+nh2,2*nh2);
	    for (i1=(-nh1+1); i1 < nh1; i1++) {
		n = 0;
		for (j=0; j < n2; j++) {
		    for (i=0; i < n1; i++) {

			if (i+i1*dh1 < n1 && i+i1*dh1 >=0 && 
			    j+i2*dh2 < n2 && j+i2*dh2 >=0) {
			    vari[(i2+nh2-1)*(2*nh1-1)+i1+nh1-1] += 
				(fbuf[(j+i2*dh2)*n1+(i+i1*dh1)] - fbuf[j*n1+i])
				*(fbuf[(j+i2*dh2)*n1+(i+i1*dh1)] - fbuf[j*n1+i]);
			    n++;
			}
		    }
		}
		if (semi) {
		    vari[(i2+nh2-1)*(2*nh1-1)+i1+nh1-1] /= (n*2.+FLT_EPSILON);
		} else {
		    vari[(i2+nh2-1)*(2*nh1-1)+i1+nh1-1] /= (n*1.+FLT_EPSILON);
		}
	    }
	}
	sf_floatwrite(vari,(2*nh1-1)*(2*nh2-1),out);
    }
    
    exit(0);
}

/* 	$Id$	 */

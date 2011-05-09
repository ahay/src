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
    int i, j, n, i1, n1, i2, n2, i3, n3, dh1, nh1, oh1, dh2, nh2, oh2;
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
    /* interval (number) of variogram lag in first axis */
    if (!sf_getint("oh1",&oh1)) oh1=0;
    /* origin (number) of variogram lag in first axis */
    if (oh1 < 0 || oh1 >= n1/dh1) oh1=0;
    if (!sf_getint("nh1",&nh1)) nh1=n1/dh1-oh1;
    /* number of variogram lag in first axis */

    if (!sf_getint("dh2",&dh2)) dh2=1;
    /* interval (number) of variogram lag in second axis */
    if (!sf_getint("oh2",&oh2)) oh2=0;
    /* origin (number) of variogram lag in second axis */
    if (oh2 < 0 || oh2 >= n2/dh2) oh2=0;
    if (!sf_getint("nh2",&nh2)) nh2=n2/dh2-oh2;
    /* number of variogram lag in second axis */

    if (!sf_getbool("semi",&semi)) semi=true;
    /* if y, output semivariogram */

    sf_settype(out,SF_FLOAT);
    sf_putint(out,"n1",nh1);
    sf_putfloat(out,"o1",oh1*dh1*d1);
    sf_putfloat(out,"d1",dh1*d1);
    sf_putstring(out,"label1","X Lag");

    sf_putint(out,"n2",nh2);
    sf_putfloat(out,"o2",oh2*dh2*d2);
    sf_putfloat(out,"d2",dh2*d2);
    sf_putstring(out,"label2","Y Lag");

    fbuf = sf_floatalloc(n1*n2);
    vari = sf_floatalloc(nh1*nh2);

    for (i3 = 0; i3 < n3; i3++) {
	sf_floatread(fbuf,n1*n2,in);
	for (i1=0; i1 < nh1*nh2; i1++) {
	    vari[i1]=0.;
	}	

      	for (i2=0; i2 < nh2; i2++) {
	    if (verb)
		sf_warning("Slice %d of %d;",i2+1,nh2);
	    for (i1=0; i1 < nh1; i1++) {
		n = 0;
		for (j=0; j < n2; j++) {
		    for (i=0; i < n1; i++) {

			if (i+(i1+oh1)*dh1 < n1 && j+(i2+oh2)*dh2 < n2) {
			    vari[i2*nh1+i1] += (fbuf[(j+(i2+oh2)*dh2)
*n1+(i+(i1+oh1)*dh1)] - fbuf[j*n1+i])
				*(fbuf[(j+(i2+oh2)*dh2)*n1+(i+(i1+oh1)*dh1)] - fbuf[j*n1+i]);
			    n++;
			}
		    }
		}
		if (semi) {
		    vari[i2*nh1+i1] /= (n*2.+FLT_EPSILON);
		} else {
		    vari[i2*nh1+i1] /= (n*1.+FLT_EPSILON);
		}
	    }
	}
	sf_floatwrite(vari,nh1*nh2,out);
    }
    
    exit(0);
}

/* 	$Id$	 */

/* Apply power gain. 

March 2013 program of the month:
http://www.ahay.org/rsflog/index.php?/archives/327-Program-of-the-month-sfpow.html
*/
/*
  Copyright (C) 2004 University of Texas at Austin

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

#include <math.h>
#include <rsf.h>

int main(int argc, char *argv[])
{
    int n[SF_MAX_DIM], ii[SF_MAX_DIM], j, nd, id, ix, nx, nbuf;
    off_t i, nsiz;
    float d, o, p, *gain[SF_MAX_DIM], *buf=NULL;
    char key[6], *prog=NULL;
    sf_file in=NULL, out=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    nd = sf_filedims (in, n);
    nsiz = 1;

    nbuf = BUFSIZ/sizeof(float);
    buf = sf_floatalloc(nbuf);

    for (id=0; id < nd; id++) {
	nsiz *= n[id];

	(void) snprintf(key,5,"pow%d",id+1);
	if (!sf_getfloat(key,&p)) p = 0.;
	/*( pow#=(0,0,...) power on #-th axis )*/

	if (p != 0.) {
	    nx = n[id];
	    gain[id] = sf_floatalloc(nx);

	    (void) snprintf(key,3,"d%d",id+1);
	    if (!sf_histfloat(in,key,&d)) d = 1.;

	    (void) snprintf(key,3,"o%d",id+1);
	    if (!sf_histfloat(in,key,&o)) o = 0.;

	    for (ix=0; ix < nx; ix++) {
		gain[id][ix] = powf(o+(ix+1)*d,p);
	    }
	} else {
	    gain[id] = NULL;
	}
    }

    prog = sf_getprog();
    if (NULL != strstr(prog,"tpow")) {
	if (!sf_getfloat("tpow",&p)) p = 0.;
	/* power on the first axis */

	if (p != 0. && NULL == gain[0]) {
	    nx = n[0];
	    gain[0] = sf_floatalloc(nx);

	    if (!sf_histfloat(in,"d1",&d)) d = 1.;
	    if (!sf_histfloat(in,"o1",&o)) o = 0.;

	    for (ix=0; ix < nx; ix++) {
		gain[0][ix] = powf(o+(ix+1)*d,p);
	    }
	} else {
	    gain[0] = NULL;
	}
    }

    for (i=0; nsiz > 0; nsiz -= nbuf) {
	if (nbuf > nsiz) nbuf = nsiz;
	sf_floatread(buf,nbuf,in);

	for (j=0; j < nbuf; j++, i++) {
	    sf_line2cart(nd,n,i,ii);
	    for (id=0; id < nd; id++) {
		if (NULL != gain[id]) buf[j] *= gain[id][ii[id]]; 
	    }
	}

	sf_floatwrite(buf,nbuf,out);
    }

    exit(0);
}

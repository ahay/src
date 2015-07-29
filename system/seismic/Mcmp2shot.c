/* Convert CMPs to shots for regular 2-D geometry. */
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

#include <rsf.h>

int main(int argc, char* argv[])
{
    int nt,ns, ny,nh, iy,ih,is,it, type;
    off_t pos;
    bool sign;
    float ds, dy,dh, os, oy,oh;
    char *trace, *zero;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nh)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&ny)) sf_error("No n3= in input");

    if (!sf_histfloat(in,"d2",&dh)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"d3",&dy)) sf_error("No d3= in input");
    if (!sf_histfloat(in,"o2",&oh)) sf_error("No o2= in input");
    if (!sf_histfloat(in,"o3",&oy)) sf_error("No o3= in input");

    if (!sf_getbool("positive",&sign)) sign=true;
    /* initial offset orientation */

    type = 0.5 + dh/dy;

    ds = dh;
    os = sign? oy - oh - (nh-1)*dh: oy + oh;
    ns = (ny-1)/type + nh;

    sf_putint(out,"n2",type*nh);
    sf_putfloat(out,"d2",dh/type);

    sf_putint(out,"n3",ns);
    sf_putfloat(out,"d3",ds);
    sf_putfloat(out,"o3",os);

    sf_putstring(out,"label3","Shot");

    nt *= sf_esize(in);

    trace = sf_charalloc(nt);
    zero = sf_charalloc(nt);
    memset(zero,0,nt);

    sf_fileflush(out,in);
    sf_setform(in,SF_NATIVE);
    sf_setform(out,SF_NATIVE);

    sf_unpipe(in,(off_t) ny*nh*nt);
    pos = sf_tell(in);

    for (is=0; is < ns; is++) {
	for (ih=0; ih < nh; ih++) {
	    for (it=0; it < type; it++) {
		iy = sign? it + type*(is + ih - nh + 1): type*(is - ih) - it;
		if (iy >= 0 && iy < ny) {
		    sf_seek(in,pos+(off_t)(iy*nh+ih)*nt,SEEK_SET);
		    sf_charread(trace,nt,in);
		    sf_charwrite(trace,nt,out);
		} else {
		    sf_charwrite(zero,nt,out);
		}
	    }
	}
    }


    exit(0);
}

/* 	$Id: Mcmp2shot.c 7107 2011-04-10 02:04:14Z ivlad $	 */

/* Find the position of first non-zero value along time axis. */
/*
  Copyright (C) 2022 Jilin University
  
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
#include <math.h>
int main(int argc, char *argv[])
{
    bool  verb;

    int   n1, n2, i1, i2, *pos;                    
    float d1, o1, s1;
    float *sig;
    sf_file inp, out;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    /* verbosity */
    if (!sf_getbool("verb",&verb)) verb = false;

    if (!sf_histint(inp, "n1", &n1)) sf_error("No n1=in input");
    if (!sf_histfloat(inp, "d1", &d1)) d1 = 1. ;
    if (!sf_histfloat(inp, "o1", &o1)) o1 = 0. ;

    n2=sf_leftsize(inp,1);
    
    sf_unshiftdim(inp, out, 1);
    
    sf_putint(out,"n1",n2);
    sf_putfloat(out,"d1",1);
    sf_putfloat(out,"o1",o1);

    sf_settype(out,SF_INT);
    
    pos=sf_intalloc(n2);
    sig=sf_floatalloc(n1);
    
    for (i2=0; i2<n2; i2++) {
	
	sf_floatread(sig,n1,inp);
	
	sf_warning("slice %d of %d;",i2+1,n2);
	
	for (i1=0; i1<n1; i1++) {
	    s1 = sig[i1];
	    if (s1 == 0) {
		pos[i2] = 0;
		pos[i2] = i1;
	    } else break;
	}
    }
    sf_intwrite(pos,n2,out);

    exit(0);
}

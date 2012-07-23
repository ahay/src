/* S transform */
/*
  Copyright (C) 2010 University of Texas at Austin
   
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

#include "st.h"

int main(int argc, char *argv[])
{
    int i, n1, n2, nw, nflo, nfhi, nf;
    float *inp, d1, d2, o2, flo, fhi;
    bool inv;
    sf_complex *outp;
    sf_file in, out;

    sf_init(argc,argv);

    in  = sf_input("in");
    out = sf_output("out");
    
    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, do inverse transform */

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&d1)) d1=0.004;

    if (!inv) {
		n2 = sf_leftsize(in,1);
	
		if (!sf_getfloat("flo",&flo)) {
			/* Low frequency in band, default is 0 */
			flo=0.;
		} else if (0. > flo) {
			sf_error("Negative flo=%g",flo);
		} else {
			flo *= d1;
		}
	
		if (!sf_getfloat("fhi",&fhi)) {
			/* High frequency in band, default is Nyquist */	
			fhi=0.5;
			if (flo/d1 > fhi/d1) 
				sf_error("Need flo < fhi, "
						 "got flo=%g, fhi=%g(Nyquist)",flo/d1,fhi/d1);
		} else {
			fhi *= d1;	
			if (flo > fhi) 
				sf_error("Need flo < fhi, "
						 "got flo=%g, fhi=%g",flo/d1,fhi/d1);
			if (0.5 < fhi)
					sf_error("Need fhi < Nyquist, "
							 "got fhi=%g, Niquist=%g",fhi/d1,0.5/d1);
		}
	
		nflo = (int) (flo*n1+0.5);
		nfhi = (int) (fhi*n1+0.5);
		nf = nfhi-nflo+1;
		sf_shiftdim(in, out, 1);
		sf_putint(out,"n2",nf);
		sf_putfloat(out,"d2",1/(d1*n1));
		sf_putfloat(out,"o2",flo/d1);
		sf_putstring(out,"label2","Frequency");
		sf_putstring(out,"unit2","Hz");
		sf_settype(out,SF_COMPLEX);
    } else {
		if (!sf_histint(in,"n2",&nw)) sf_error("No n2= in input");
		if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input");
		if (!sf_histfloat(in,"o2",&o2)) sf_error("No o2= in input");
		flo = o2*d1;
		fhi = (o2+(nw-1)*d2)*d1;
		n2 = sf_leftsize(in,2);
		nflo = (int) (flo*n1+0.5);
		nfhi = (int) (fhi*n1+0.5);
		nf = nfhi-nflo+1;
		if(nw!=nf) {
			sf_warning("n2!=nf, (n2=%d, nf=%d)",nw,nf);
			nf = nw;
		}
		sf_unshiftdim(in, out, 2);
		sf_settype(out,SF_FLOAT);
    }
    
    inp = sf_floatalloc(n1);
    outp = sf_complexalloc(n1*nf);
    
    for (i=0; i < n2; i++)  {
		sf_warning("slice %d of %d;",i+1,n2);
		if (!inv) {
			sf_floatread(inp,n1,in);
			st(n1,d1,nflo,nfhi,inp,outp);
			sf_complexwrite(outp,n1*nf,out);
		} else {
			sf_complexread(outp,n1*nf,in);
			ist(n1,d1,nflo,nfhi,inp,outp);
			sf_floatwrite(inp,n1,out);
		}
    }
    sf_warning(".");
    exit(0);
}

/* 	$Id$	 */

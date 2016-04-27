/* Pre-stack 2-D oriented velocity continuation. 

   Axes: (Omega,h,p,k) -> (Omega,v,p,k)

   Make sure you use half-offsets for h.

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
#include <rsf.h>

int main(int argc, char* argv[])
{
    bool verb;
    int nx,nv,np,nw, ix,iv,ip,iw,iw2, nh,ih;
    float v0,v1,v2,v,dv, dx,dp,dw, x0,p0,w0, x,p,w, dh,h,h0;
    sf_complex *ctrace, *ctrace2, shift, **cstack, *trace;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");
    if (!sf_histint(in,"n1",&nw)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nh)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&np)) sf_error("No n3= in input");
    if (!sf_histint(in,"n4",&nx)) sf_error("No n4= in input");
    
    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity flag */

    if (!sf_histfloat(in,"o1",&w0)) w0=0.;  
    if (!sf_histfloat(in,"d1",&dw)) sf_error("No d1= in input");

    if (!sf_getint("nv",&nv)) sf_error("Need nv=");
    /* velocity steps */
    if (!sf_getfloat("dv",&dv)) sf_error("Need dv=");
    /* velocity step size */
    if (!sf_getfloat("v0",&v0) && 
	!sf_histfloat(in,"v0",&v0)) sf_error("Need v0=");
    /*( v0 starting velocity )*/

    if(!sf_histfloat(in,"d2",&dh)) sf_error("No d2= in input");
    if(!sf_histfloat(in,"d3",&dp)) sf_error("No d3= in input");
    if(!sf_histfloat(in,"d4",&dx)) sf_error("No d4= in input");
    
    if(!sf_histfloat(in,"o2",&h0)) sf_error("No n2= in input");
    if(!sf_histfloat(in,"o3",&p0)) sf_error("No o3= in input");
    if(!sf_histfloat(in,"o4",&x0)) x0=0.;

    sf_putfloat(out,"o2",v0+dv);
    sf_putfloat(out,"d2",dv);
    sf_putint(out,"n2",nv);

    sf_putstring(out,"label2","Velocity");

    //sf_shiftdim(in, out, 2);

    //sf_putfloat(out,"o5",h0);
    //sf_putfloat(out,"d5",dh);
    //sf_putint(out,"n5",nh);

    dx *= 2.*SF_PI;
    x0 *= 2.*SF_PI;

    dw *= 2.*SF_PI;
    w0 *= 2.*SF_PI;

    ctrace  = sf_complexalloc(nw);
    ctrace2 = sf_complexalloc(nw);
    cstack = sf_complexalloc2(nw,nv);
    trace  = sf_complexalloc(nw);
    
for (ix=0; ix < nx; ix++) {
	x = x0+ix*dx;
	sf_warning("x %d of %d;\n", ix,nx); 
       
    	for (ip=0; ip < np; ip++) {
	//if (verb) sf_warning("slope %d of %d;", ip+1,np);
		p = p0+ip*dp;
		for (iv=0; iv < nv; iv++) {
			for (iw2=0; iw2 < nw; iw2++) {

				cstack[iv][iw2]= sf_cmplx(0.,0.);

		}}
		
	
		for (ih=0; ih < nh; ih++) {
	    		//if (verb) sf_warning("offset %d of %d;\n", ih+1,nh);
	    		h = h0 + ih*dh;
			sf_complexread(ctrace,nw,in);
	
			for (iv=0; iv < nv; iv++) {
				v = v0 + (iv+1)*dv;
				v1 = 4./(v*v) - 4./(v0*v0);
				v2 = 0.25*v*v - 0.25*v0*v0;  
                	//V2 equals to zero is unaccepatble!!!    
			
				for (iw=0; iw < nw; iw++) {
				    w = iw*dw;
				    w = (0.25*w*p*p + 0.5*x*p)*v2 + w*h*h*v1;
				    
				    //sf_warning("w %f\n", w);
                                    //sf_warning("offset %d of %d;\n", ih+1,nh);
				    //sf_warning("v is %f\n", v);				    
				    //sf_warning("v1 is %f\n", v1);
				    //sf_warning("v2 is %f\n", v2);
				    //sf_warning("w is %f\n", w);
                                    
				    shift = sf_cmplx(cosf(w),sinf(w));
				    //shift = sf_cmplx(1.,0.);
#ifdef SF_HAS_COMPLEX_H
				    ctrace2[iw] = ctrace[iw] * shift;
#else
				    ctrace2[iw] = sf_cmul(ctrace[iw],shift);
#endif
				    cstack[iv][iw] += ctrace2[iw];
				} /* w */
			} /* v */
		} /* h */
			for (iv=0; iv < nv; iv++){
                        trace[0] = sf_cmplx(0.,0.);
			for (iw=0; iw < nw; iw++) {
			trace[iw]=cstack[iv][iw];}
			sf_complexwrite(trace,nw,out);}
	} /* p */
} /* x */
    if (verb) sf_warning(".");

    exit (0);
}

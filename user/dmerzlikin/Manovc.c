/* Oriented anisotropy continuation: shifted hyperbola travel-time approximation. 

Axis order: t, p, x
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
//#include "warp3.h"

int main(int argc, char* argv[])
{
    bool lagrange, debug, plus, isotr, testwarp, full;
    int it, nt, ip, np, ix, nx, ntpx, is, ns;
    float eps, epsr, t, t0, dt, p, p0, dp, x, x0, dx, v0;//, v, dv, dv2, sq;
    float ***slice, ***slice0, ***tstr, ***pstr, ***xstr;
    float kappa1, kappa2, kappa3, alpha=-9999.9, beta=-9999.9; //coefficients
    float smax;
    float pwarp1, pwarp2, pwarp3;
    double root2d = -9999.9;
    float root2, droot2dp=-9999.9;
    float s0,s,ds;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    
    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&np)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&nx)) sf_error("No n3= in input");

    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&dp)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"d3",&dx)) sf_error("No d3= in input");

    if (!sf_histfloat(in,"o1",&t0)) t0=0.;
    if (!sf_histfloat(in,"o2",&p0)) p0=0.;
    if (!sf_histfloat(in,"o3",&x0)) x0=0.;

    ntpx = nt*np*nx;

    if (!sf_getfloat("eps",&eps)) eps=0.1;
    /* stretch regularization */

    //if (!sf_getint("nv",&nv)) nv=1; /* number of velocity steps */

    if (!sf_getbool("lagrange",&lagrange)) lagrange=false;
    /* Use Lagrangian method */

    if (!sf_getbool("plus",&plus)) plus=true;
    /* Plus or minus in coefficients: I have two versions */

    if (!sf_getbool("debug",&debug)) debug=true;
    /* Implement debugger: add it later */

    if (!sf_getbool("isotr",&isotr)) isotr=false;
    /* Implement debugger: add it later */

    if (!sf_getbool("testwarp",&testwarp)) testwarp=false;
    /* Implement debugger: add it later */

    if (!sf_getbool("full",&full)) full=false;
    /* full accuracy flag - considers all (s-1) terms in any power */

    if (!sf_getfloat("v0",&v0)){
		if (v0==0.0 || v0<0.0) sf_error("Need v0>0");
	} /* starting velocity */

    if (!sf_getint("ns",&ns)) ns=1;
    /* s steps */
    if (!sf_getfloat("ds",&ds)) sf_error("Need ds=");
    /* s step size */
    if (!sf_getfloat("s0",&s0)) {s0=1.0; sf_warning("We assume s0 is 1");}
    /*( s0 start )*/
    //if (!sf_getfloat("smax",&smax)) {sf_error("Need smax=");}

    if (!sf_getfloat("epsr",&epsr)) epsr=0.001;
    /* damper for root */

    //ds = (smax - s0)/ns;

    sf_putint(out,"n4",ns);
    sf_putfloat(out,"o4",s0+ds);
    sf_putfloat(out,"d4",ds);

    sf_putstring(out,"label4","S heterogeneity");
    
    tstr   = sf_floatalloc3(nt,np,nx);
    pstr   = sf_floatalloc3(nt,np,nx);
    xstr   = sf_floatalloc3(nt,np,nx); 

    warp3_init(nt, t0, dt,
	       np, p0, dp,
               nx, x0, dx,
	       nt, np, nx, eps); 

    /*if (lagrange) {*/
	//slice0 = sf_floatalloc3(nt,np,nx);
	//slice  = sf_floatalloc3(nt,np,nx);
    /*} else {*/
	slice  = sf_floatalloc3(nt,np,nx);	
	slice0 = slice;
    /*} */
		
    sf_floatread(slice0[0][0],ntpx,in);

    for (is=0; is < ns; is++) {
	sf_warning("step %d of %d;",is,ns);

	s = s0+is*ds;

	/*if (lagrange) {
	    dv2 = 0.25*(v*v-v0*v0);
	} else {
	    dv2 = 0.25*(2*v-dv)*dv;
	}*/
	    
	for (ix=0; ix < nx; ix++) {
	    
	    for (ip=0; ip < np; ip++) {

		for (it=0; it < nt; it++) {

		    x = x0+ix*dx;
		    p = p0+ip*dp;
		    t = t0+it*dt;    
                    
                    if(!full){    
                    //evaluating coefficients			
                    alpha = 3.0 - 2.0*s;

                    beta  = 2.0 - 3.0*s + s*s;

                    root2d = 1 + p*p*(alpha)*v0*v0 + p*p*p*p*(beta)*v0*v0*v0*v0;

		    if( (ix == (int)nx/2) && (ip == (int)np/2+30) && (it == (int)nt/2) && debug){
                    	sf_warning("root2d = root2^2 = %3.3f\n\n",root2d);
			}

                    root2d = fabs(root2d);

		    if( (ix == (int)nx/2) && (ip == (int)np/2+30) && (it == (int)nt/2) && debug){
                    	sf_warning("after abs root2d = root2^2 = %3.3f\n\n",root2d);
			}
				
				//if (debug) nantracker((float)root2d,1,nb,ip,p,ix,is,s,it);
				
                    root2 = (float)sqrt(root2d);

				//if (debug) nantracker(root2,2,nb,ip,p,ix,is,s,it);

                    // derivative over p
                    if (plus){

                                kappa1 = p*t*v0*v0*(alpha + (alpha + 2.0*beta*v0*v0*p*p)/(root2 + epsr));

                    }/*if*/ else {

				kappa1 = p*t*v0*v0*(alpha - (alpha + 2.0*beta*v0*v0*p*p)/(root2 + epsr));

                    }//else

				//if (debug) nantracker(kappa1,3,nb,ip,p,ix,is,s,it);

                     // function F

                     if (plus){			
	
				kappa3 = 0.5*t*( (2.0 + p*p*(alpha)*v0*v0) + 2.0*root2);

                     }/*if*/ else {
				
				kappa3 = 0.5*t*( (2.0 + p*p*(alpha)*v0*v0) - 2.0*root2);
					
                     }//else

				//if (debug) nantracker(kappa3,5,nb,ip,p,ix,is,s,it);

                     // derivative over t

                     if (plus){

				kappa2 = 0.5*( (2.0 + p*p*(alpha)*v0*v0) + 2.0*root2); 

                     }/*if*/ else {

				kappa2 = 0.5*( (2.0 + p*p*(alpha)*v0*v0) - 2.0*root2); 				
 
                     }//else

				//if (debug) nantracker(kappa2,4,nb,ip,p,ix,is,s,it);

                     }/*first order approximation */ else {
				

				root2 = sqrt(1.0 + (2.0-s)*p*p*v0*v0);

				droot2dp = p*(2.0-s)*v0*v0/(root2 + epsr);

				/* dF/dp */
				kappa1 = - (   (1.0 - root2) + epsr )/ (   ( 2.0 - s )*(1.0 - (s - 1.0)*root2)  + epsr  );

				kappa1 += (   (1.0 - root2)*(1.0 - root2)*(s - 1.0) + epsr   )/(   2.0*( 2.0 - s )*( 1.0 - (s-1.0)*root2 )*( 1.0 - (s-1.0)*root2 )  + epsr );
				
				kappa1 *= t*droot2dp;

				/* dF/dt */				
				kappa2 = (   (1.0 - root2)*(1.0 - root2) + epsr   )/(   epsr + 2.0*(2.0-s)*( 1.0 - (s-1.0)*root2 )   );
				
				/* F */
				kappa3 = t*(   (1.0 - root2)*(1.0 - root2) + epsr   )/(   epsr + 2.0*(2.0-s)*( 1.0 - (s-1.0)*root2 )   );
			
                     }	    

		     	tstr[ix][ip][it] = t + (kappa3 - kappa1*p)*ds;
		     	pstr[ix][ip][it] = p + (kappa2*p)*ds;
		     	xstr[ix][ip][it] = x + (-kappa1)*ds;	

                          if( (ix == (int)nx/2) && (ip == (int)np/2+30) && (it == (int)nt/2) && debug){

                            sf_warning("is = %d [ ix = %d ; ip = %d ; it = %d ]\n",is,ix,ip,it);
                            sf_warning("s = %3.3f [ x = %3.3f ; p = %3.3f ; t = %3.3f ]\n",s,x,p,t);
			    sf_warning("kappa3 = F = 0.5*t*( (2.0 + p*p*(alpha)*v0*v0) - 2.0*root2)\n");
			    sf_warning("kappa3 = F = 0.5*%3.3f*((2.0 + %3.3f*%3.3f*%3.3f*%3.3f*%3.3f\n",t,p,p,alpha,v0,v0);
			    sf_warning("+/- 2.0*%3.3f\n",root2);
                            sf_warning("root2d = root2^2 = %3.3f =  1 + %3.3f^2*(%3.3f)*%3.3f^2 + %3.3f^4*(%3.3f)*%3.3f^4\n",root2d,p,alpha,v0,p,beta,v0);
                            sf_warning("root2d = root2^2 =  1 + p*p*(alpha)*v0*v0 + p*p*p*p*(beta)*v0*v0*v0*v0\n");
			    sf_warning("[ t = %3.3f ; x = %3.3f ; p = %3.3f]\n",tstr[ix][ip][it],xstr[ix][ip][it],pstr[ix][ip][it]);

                            sf_warning("kappa1=%3.3f kappa2=%3.3f kappa3=%3.3f ds=%3.3f\n",kappa1,kappa2,kappa3,ds);  

			    sf_warning("dt/ds=(kappa3 - kappa1*p)*ds = %3.5f\n",(kappa3 - kappa1*p)*ds);  
                            sf_warning("dp/dt=(kappa2*p)*ds = %3.5f\n",(kappa2*p)*ds);  
			    sf_warning("dx/dt=(-kappa1)*ds = %3.5f\n",(-kappa1)*ds);  

                          }// output



		}//t
	    }//x
	}//p

	warp3(slice0,tstr,pstr,xstr,slice);	
	sf_floatwrite (slice[0][0],ntpx,out);

    }
    sf_warning(".");

    exit(0);
}

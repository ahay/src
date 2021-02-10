/* Slope-based tau-p velocity transform for VTI media. */
/*
  Copyright (C) 2010 Politecnico di Milano
 
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
#include <float.h>
#include <rsf.h>


int main (int argc, char* argv[])
{
    bool map;
    sf_map4 nmo;
    char *method;
    int it,ix,ip,iv, nt,nx, np, nv, ntv, nw, nvh, ntvh, ne, ie, nte;
    float dt, t0, p, p0, t, dp, v0, dv, vh0, dvh , e0 , de;
    float Vn2_eff, eta_eff,S_eff,dtau0_dtau,dVn2_dtau,deta_dtau, A, B, C; 
    float N, D, Nt, Dt;
    float *v1, *v2, *v3; /* to compute velocity attribute panels */
    float *R, *Rt, *Q, *Qt, **coord1, **coord2, **coord3, *ord;
    float *vN, *vN2, *vH, *vH2, *e, *e2, *TAU0, *TAU0t;
    sf_file cmp=NULL, velN=NULL, velH=NULL, dip=NULL, dipt=NULL, curv=NULL;
    sf_file tau0=NULL, tau0t=NULL, curvt=NULL, eta=NULL;
    	
    sf_init (argc,argv);
    tau0 = sf_input("in");

    velN = sf_output("out");
    velH = sf_output("velH"); /*horizontal velocity*/
    eta = sf_output("eta"); /*eta*/

    if (!sf_getbool("map",&map)) map=false;
    /* output maps instead of coherency panels */

    if (!map) cmp = sf_input("cmp");

    if (NULL == (method = sf_getstring("method"))) method="stripping";
    /* method to use (stripping,dix,fowler,effective) */
    
    if (SF_FLOAT != sf_gettype(tau0)) sf_error("Need float input");

    if (!sf_histint(tau0,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(tau0,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(tau0,"o1",&t0)) sf_error("No o1= in input");
    
    if (!sf_histint(tau0,"n2",&np)) sf_error("No n2= in input");
    if (!sf_histfloat(tau0,"d2",&dp)) sf_error("No d2= in input");
    if (!sf_histfloat(tau0,"o2",&p0)) sf_error("No o2= in input");
    p0 /= dp; /* slope axis is normalized: no need for slope normalization dp=1*/

    v1  = sf_floatalloc(nt);
    v2  = sf_floatalloc(nt);
    v3  = sf_floatalloc(nt);
    
    if (map) {    
	nmo = sf_stretch4_init (nt, t0, dt, nt, 0.01);

	ord = NULL;
	coord1 = NULL;
	coord2 = NULL;
	coord3 = NULL;

	ntv = 0;
	ntvh= 0;
	nte = 0;	

	vN  = NULL;
	vN2 = NULL;
	vH  = NULL;
	vH2 = NULL;
	e   = NULL;
	e2  = NULL;
    } else {
	nmo = NULL;

	ord = sf_floatalloc(nt);	
	coord1 = sf_floatalloc2(2,nt);
	coord2 = sf_floatalloc2(2,nt);
	coord3 = sf_floatalloc2(2,nt);    

	/* first velocity... VNMO */
	if (!sf_getint("nv",&nv)) sf_error("Need nv=");     /* number of velocities */
	if (!sf_getfloat("v0",&v0)) sf_error("Need v0=");     /* velocity origin */
	if (!sf_getfloat("dv",&dv)) sf_error("Need dv=");     /* velocity sampling */
	/* second velocity... VH */
	if (!sf_getint("nvh",&nvh)) nvh=nv;     /* number of HOR velocities  */
	if (!sf_getfloat("vh0",&vh0)) vh0=v0;     /* HOR velocity origin */
	if (!sf_getfloat("dvh",&dvh)) dvh=dv;     /* HOR velocity sampling */ 
	/* anellipticity */
	if (!sf_getint("ne",&ne))   ne = 101;     /* number of etas */
	if (!sf_getfloat("e0",&e0)) e0 = -0.5;    /* eta origin */
	if (!sf_getfloat("de",&de)) de=0.01;      /* eta sampling */ 
    
    
	/* adding the second dimension to the ouput files*/
	sf_putint(velN,"n2",nv);
	sf_putfloat(velN,"o2",v0);
	sf_putfloat(velN,"d2",dv);
	
	sf_putint(velH,"n2",nvh);
	sf_putfloat(velH,"o2",vh0);
	sf_putfloat(velH,"d2",dvh);
	
	sf_putint(eta,"n2",ne);
	sf_putfloat(eta,"o2",e0);
	sf_putfloat(eta,"d2",de);

	ntv = nt*nv;
	ntvh= nt*nvh;
	nte = nt*ne;	

	vN  = sf_floatalloc(ntv);
	vN2 = sf_floatalloc(ntv);
	vH  = sf_floatalloc(ntvh);
	vH2 = sf_floatalloc(ntvh);  
	e  = sf_floatalloc(nte);
	e2 = sf_floatalloc(nte);
    }
    
    /* reading the number of cmp in the data */
    nx = sf_leftsize(tau0,2);
    
    if (!sf_getint("nw",&nw)) nw=4;
    /* interpolator size (2,3,4,6,8) */

    /* Auxiliary inputs */

    TAU0 = sf_floatalloc(nt);

    if ('e' == method[0] || 'E' == method[0] ||
	'd' == method[0] || 'D' == method[0]) {
	if (NULL != sf_getstring("dip")) dip = sf_input("dip"); 
        /*slope field (required for method=e and method=d) */
	else sf_error("Need dip field");
	if (NULL != sf_getstring("curv")) curv = sf_input("curv"); 
        /*curvature field (required for method=e and method=d) */
	else sf_error("Need curvature field");
	R = sf_floatalloc(nt);
	Q = sf_floatalloc(nt);
    } else {
	R = NULL;
	Q = NULL;
    }

    if ('e' == method[0] || 'E' == method[0]) {
	Rt = NULL;
	Qt = NULL;
	TAU0t = NULL;
    } else {
	if (NULL != sf_getstring("dipt")) dipt = sf_input("dipt"); 
        /*time derivative of slope field*/
	else sf_error("Need dipt field=");
        Rt = sf_floatalloc(nt);    

	if ('f' == method[0] || 'F' == method[0]) {		
	    if (NULL != sf_getstring("tau0t")) tau0t = sf_input("tau0t"); 
	    /*tau0 tau derivative field (required for method=f) */
	    else sf_error("Need time derivative of tau0 field");
	    TAU0t = sf_floatalloc(nt);
	    Qt = NULL;
	} else {
	    if (NULL != sf_getstring("curvt")) curvt = sf_input("curvt"); 
            /*time derivative of curvature field (required for method=d and method=s) */
	    else sf_error("Need curvt field=");
	    TAU0t = NULL;
	    Qt = sf_floatalloc(nt);
	}
    }
    
    for (ix = 0; ix < nx; ix++) { /* CMP loop*/
	
	if (!map) {
	    for (iv=0; iv < ntv; iv++) {
		vN[iv]=0.;
	    }
	    for (iv=0; iv < ntvh; iv++) {
		vH[iv]=0.;
	    }
	    for (ie=0; ie< nte; ie++) {
		e[ie]=0.;
	    }
	}
	
        for (ip = 0; ip < np; ip++) { /* slope p loop*/
            p=p0+ip;

	    sf_floatread (TAU0, nt,tau0);
            if (!map) sf_floatread (ord, nt, cmp);	    

	    if ('e' == method[0] || 'E' == method[0] ||
		'd' == method[0] || 'D' == method[0]) {
		sf_floatread (  R, nt, dip);
		sf_floatread (  Q, nt, curv);
	    } 

	    if ('e' != method[0] && 'E' != method[0]) {
                sf_floatread (Rt, nt, dipt);
		
		if ('f' == method[0] || 'F' == method[0]) {	
		    sf_floatread (TAU0t, nt, tau0t);	
		} else {
		    sf_floatread (Qt, nt, curvt);
		}
	    }

            for (it=0; it < nt; it++) { /* time tau loop*/
                t = t0 + it*dt;
                
		if (TAU0[it] < 0.) { /* TODO: FIX THIS CHECK */
		    v1[it] = 0.;
		    v2[it] = 0.;
		    v3[it] = 0.;
                } else {
		    switch(method[0]) {
			case 's': case 'S': /* STRIPPING EQUATIONS*/            			    
			    
			    Nt=(3*Rt[it]+p*Qt[it]-3*Rt[it]*Rt[it]*p);
			    Dt=(3*Rt[it]+p*Qt[it]+1*Rt[it]*Rt[it]*p);
			    v1[it] =
				sqrtf(fabsf( 
					  (-1/p)*(16*Rt[it]*Rt[it]*Rt[it]*(1/dp/dp) )
					  /(Nt*(Dt+FLT_EPSILON) )
					  ) );
			    
                        
			    v2[it] =
				sqrtf(fabsf( (1/(p*p*dp*dp) )*
					     ( (Nt-4*Rt[it]) /(Nt) )
					  ) );
			    v3[it] = .5* ( (v2[it]*v2[it] / v1[it]/v1[it]) -1 );
			    break;
			case 'd': case 'D': /* DIX EQUATIONS*/

			    N=(3*t*R[it]*(dt)+t*p*Q[it]*(dt)-3*R[it]*R[it]*(dt*dt)*p);
			    D=(3*t*R[it]*(dt)+t*p*Q[it]*(dt)+1*R[it]*R[it]*(dt*dt)*p);
			    
			    Nt = (3*t-6*p*R[it]*dt)*Rt[it]+p*t*Qt[it]+(3*R[it]+p*Q[it])*dt;
			    Dt = (3*t+2*p*R[it]*dt)*Rt[it]+p*t*Qt[it]+(3*R[it]+p*Q[it])*dt;
                            
			    Vn2_eff = (-1/p)*(16*t*R[it]*R[it]*R[it]*(dt*dt*dt/dp/dp) ) /(N*D);
			    eta_eff = ( 1/p) * ( N * (4*t*R[it]*dt-D) ) / (32*t*R[it]*R[it]*R[it]*dt*dt*dt );        		
			    S_eff = 1+8*eta_eff;
			    
			    /*dtau0_dtau*/	
			    dtau0_dtau = TAU0[it]/t+(
				(t*t)/(2*TAU0[it])*
				(Nt*D-N*Dt)/(D*D) ); 
			    
			    
			    /*dVn2_dtau*/	
			    
			    /*dVn2_dtau_n = (-16/((P/dp)*(Nn*Dn)^2))*(ts^3/dp^2)*(...
			      ( Rn^3+3* Rn^2*TAU/ts*Rtn)*(Nn*Dn)-TAU*Rn^3*(Ntn*Dn+Nn*Dtn));*/
			    
			    dVn2_dtau = (-16/(p* (N*N*D*D) ) )*(dt*dt*dt)/(dp*dp)*( 
				(R[it]*R[it]*R[it]+3*R[it]*R[it]*t/dt*Rt[it])*(N*D)
				-t*R[it]*R[it]*R[it]*(Nt*D+N*Dt) );
			    
			    /*1/(32*P/dp*(TAU*Rn^3)^2)*(1/ts^2)*(...
			      (Ntn*(4*TAU*Rn-Dn/ts)+Nn*(4*(Rn+TAU*Rtn/ts)-Dtn/ts))*(TAU*Rn^3)-...
			      ( (Nn*(4*TAU*Rn-Dn/ts)) * (Rn^3+3*Rn^2*TAU*Rtn/ts) )... 
			      );*/
			    
			    deta_dtau = 1/( 32*p*(t*t*R[it]*R[it]*R[it]*R[it]*R[it]*R[it]) )/(dt*dt)*(
				(Nt*(4*t*R[it]-D/dt)+N*(4*(R[it]+t*Rt[it]/dt)-Dt/dt))*(t*R[it]*R[it]*R[it])-
				( (N*(4*t*R[it]-D/dt)) * (R[it]*R[it]*R[it]+3*R[it]*R[it]*t*Rt[it]/dt) )
				);
			    
			    /* compute Vn int = dtau (tau0 Vn2)/ dtau0_dtau*/
			    /* sqrt(abs( (dtau0_dtau*Vn2_eff+tau0*dVn2_dtau)/dtau0_dtau_n) ) */
			    v1[it] =
				sqrtf(fabsf( (dtau0_dtau*Vn2_eff+TAU0[it]*dVn2_dtau) / dtau0_dtau ) ) ;
			    
			    /* tempvec[it] = v1[it]; */
			    
			    
			    /*A = S_eff*Vn2_eff/(Vn(it,ip)^2);
			      B = tau0(it,ip)*(Vn2_eff/(Vn(it,ip)^4))/dtau0_dtau;
			      C = (Vn2_eff*8*deta_dtau+S_eff*dVn2_dtau);*/
			    
			    /* compute eta int */
			    A = (S_eff*Vn2_eff/v1[it]/v1[it]);
			    B = TAU0[it]*Vn2_eff/v1[it]/v1[it]/v1[it]/v1[it]/dtau0_dtau;
			    C =  Vn2_eff*8*deta_dtau+S_eff*dVn2_dtau;
			    /* sf_warning("\nA %f B %f C %f R %f",A,B,C, ((A + B * C )-1)/8 ); */
			    v3[it] = ((A + B * C )-1)/8 ;
			    v2[it] = sqrtf(fabsf( v1[it]*v1[it] * (1+2*v3[it]) 	) );
			    /* sf_warning("\nVn %f eta %f",v1[it],v3[it]); */
			    break;
			case 'e': case 'E':
			    N=(3*t*R[it]*(dt)+t*p*Q[it]*(dt)-3*R[it]*R[it]*(dt*dt)*p);
			    D=(3*t*R[it]*(dt)+t*p*Q[it]*(dt)+1*R[it]*R[it]*(dt*dt)*p);

			    v1[it] = sqrtf(fabsf( 
						      (-1/p)*(16*t*R[it]*R[it]*R[it]*(dt*dt*dt/dp/dp) )
						      /(N*D)
						      ) );
			    v2[it] = sqrtf(fabsf( (1/(p*p*dp*dp)) * (
							     (N-4*t*R[it]*(dt) ) /(N) )
						      ) );
			    v3[it] = .5* ( (v2[it]*v2[it] / v1[it]/v1[it]) -1 );
			    break;
			case 'f': case 'F':
			    Rt[it] *= (-dt/dp);

			    v1[it] =
				sqrtf(fabsf( 
					  ((TAU0t[it]*TAU0t[it])-(dt*dt))*((TAU0t[it]*TAU0t[it])-(dt*dt))/
					  (TAU0t[it]*TAU0t[it]*p*p*p*Rt[it]*dt*dp*dp*dp )
					  ) );
                        
			    /* tempvec[it]= v1[it]; */
			    
			    v2[it] =
				sqrtf(fabsf( 
					  ( (TAU0t[it]*TAU0t[it]*(p*dp*Rt[it]-dt)) + (dt*dt*dt) )/
					  ( p*p*p*TAU0t[it]*TAU0t[it]*Rt[it]*dp*dp*dp )	
					  ) );
             
			    v3[it] = .5* ( (v2[it]*v2[it] / v1[it]/v1[it]) -1 );
			    break;			    
			default:
			    sf_error("Unknown method %s",method);
			    break;
		    }   /* end of method switch */  
		} /* end of it tau0 < 0 */                  
 
		if (!map) {
                    coord1[it][0] = TAU0[it];
                    coord2[it][0] = TAU0[it];
                    coord3[it][0] = TAU0[it];

		    coord1[it][1] = v1[it];
		    coord2[it][1] = v2[it];
		    coord3[it][1] = v3[it];
		}		
            } /* END tau t loop */

	    if (map) {
		sf_stretch4_define (nmo,TAU0,false);
	    
		sf_stretch4_apply (false,nmo,v1,v1); sf_floatwrite (v1,nt,velN);	    
		sf_stretch4_apply (false,nmo,v2,v2); sf_floatwrite (v2,nt,velH);
		sf_stretch4_apply (false,nmo,v3,v3); sf_floatwrite (v3,nt,eta);
	    } else {
		sf_int2_init (coord1, t0,v0, dt,dv, nt,nv, sf_spline_int, nw, nt);
		sf_int2_lop (true,true,ntv,nt,vN,ord);
		
		sf_int2_init (coord2, t0,vh0, dt,dvh, nt,nvh, sf_spline_int, nw, nt);
		sf_int2_lop (true,true,ntvh,nt,vH,ord);
            
		sf_int2_init (coord3, t0,e0, dt,de, nt,ne, sf_spline_int, nw, nt);
		sf_int2_lop (true,true,nte,nt,e,ord);
	    }
        } /* END slope p loop */
        
	if (!map) {
	    /* from spline coefficients to model */
	    if (nw > 2) { /* loop spline */
		for (iv=0; iv < nv; iv++) {
		    sf_spline_post (nw, iv*nt, 1, nt, vN, vN2);
		}
		for (iv=0; iv < nvh; iv++) {
		    sf_spline_post (nw, iv*nt, 1, nt, vH, vH2);
		}
		
		for (ie=0; ie < ne; ie++) {
		    sf_spline_post (nw, ie*nt, 1, nt, e, e2);
		}
		
		for (it=0; it < nt; it++) {
		    sf_spline_post (nw, it, nt, nv, vN2, vN);
		    sf_spline_post (nw, it, nt, nvh, vH2, vH);
		    sf_spline_post (nw, it, nt, ne, e2, e);
		}
	    } /* END loop spline */
        
	    sf_floatwrite (vN,ntv,velN);
	    sf_floatwrite (vH,ntvh,velH);
	    sf_floatwrite (e,nte,eta);
	}
    } /* END CMP loop */
 
    exit (0);
} /* END of main */

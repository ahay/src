/* Kirchhoff integral modeling using Newton's method*/
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
#include <math.h>

#include "vectorops.h"

#include "kirmod.h"
#include "general_traveltime.h"
/*^*/

#include "setvelocity.h"

void kirmodnewton_table(int vstatus /* Type of model (vconstant(0) or vgradient(1))*/,
			bool debug /* Debug Newton */,
			float xs /* Source */,
			float xr /* Receiver */,
			float bmin /* Min value in x-direction */,
			float bmax /* Max value in x-direction */,
			int niter /* Number of iteration for Newton */,
			double tol /* Tolerance level for Newton */,
			int n /* Number of reflection */,
			float *xinitial /* Initial guess of points of intersection */,
			float *xref /* x-coordinate reference points */,
			float *zref /* z-coorditnate reference points */,
			float *v /* Velocities at reference points */,
			float *gx /* x-gradient at the reference points */,
			float *gz /* z-gradeint at the reference points */,
			float **aniso /*anisotropy parameters*/, 
			func1 z /* z(k,x) */,
			func1 zder /* z'(k,x) */,
			func1 zder2 /* z''(k,x) */,
			float **oldans /* To export old answer*/,
			bool skip /* To tell the kirmodnewton2 to generate xinitial in the case where the previous ray can't be traced (can't fwdxini)*/,
			ktable table /* [5] output table */)
/*< Compute traveltime attributes >*/
{
    float  *xx, *xxnew, *F, *dk, *xxtem, *zk,*ck_inv, **ans;
    float tt, tx_r=0.0f, ty, tz_s=0.0f, tz_r=0.0f, tn=0.0f, at, an, v_1r=0.0f, v_r=0.0f, alpha, theta, dip=0.0f;
    int q=0; /* Counter for big loop*/
    int b3=0; /* To check whether the ray is tracable*/
	
    float dktemp=0.0f;
    int t,a,b1,b2,b4; /* Counter*/

    float value;
    func3 f;

    int i,j1,i3,i4; /* Counter*/
    float Ftem=0;

    int i2,j2,i5; /* Counter*/
    int w = niter; /* Number of loops for yk=yk-dk*/
    int l; /* Counter*/
    int m; /* Counter*/

    float ck_in, ck_in_temp=0.0,S1,S3;
    int c1,c2,c3,c4,c5;
    int c; /* Counter*/
	   
    /* Allocate space-------------------------------------------------------------------------------------*/
	
    xx = sf_floatalloc(n+2); /* Positions of intersection points to be iteratively perturbed*/
    xxnew = sf_floatalloc(n+2); /* Temporary array for xx*/
    ans = sf_floatalloc2(2,n+2);
	
    F = sf_floatalloc(n+2); /* Array of first derivative of traveltime (Snell's law at each interface)*/
    dk = sf_floatalloc(n+2); /* Changes to xx*/
    xxtem = sf_floatalloc(n+2); /* Pre-calculated xx to check the boundary condition*/
    ck_inv = sf_floatalloc(n+2);
    zk = sf_floatalloc(n+2);	

    /* Switch bmin and bmax if necessary */
	
    if(bmin>bmax) {
	value = bmin;
	bmin = bmax;
	bmax = value;
    }
	
    /* Set vconstant or vgradient----------------------------------------------------------------------*/
	
    f.T_k = 0; /* Initialize structure f*/
    f.T_k_k = 0;
    f.T_k_k1 = 0;
    f.T_k_k_k = 0;
    f.T_k_k1_k1 = 0;
    f.T_k_k_k1 = 0;
    f.T_k_zk = 0;
    f.T_k_zk1 = 0;
    f.T_k_zk_zk = 0;
    f.T_k_zk1_zk1 = 0;
    f.T_k_zk_zk1 = 0;
    f.T_k_k_zk = 0;
    f.T_k_k1_zk1 = 0;
    f.T_k_k_zk1 = 0;
    f.T_k_k1_zk = 0;
	
    setfunc(vstatus,&f); /* Set value of structure f*/ 

    /* Step 1: Calculate F(y) to see if it is sufficiently close to zero----------------------------------*/
	
    xx[0] = xs;
    xx[n+1] = xr;
	
    if (n!=0) {
	for (i3=0; i3<n; i3++) {
	    xx[i3+1] = xinitial[i3]; /* Create an array of points of intersection from soure to receiver*/
	}
    }
    else {
	b3 = 0;
	goto mark;
    }
    
    for (i=0; i<n; i++) {
	initialize(i+1,n,xx,v,xref,zref,gx,gz,aniso,z,zder,zder2); /*Initialize y_1k, y_k and y_k1*/
	F[i+1] = T_hat_1k_k(f.T_k_k1,f.T_k_zk1) + T_hat_k_k(f.T_k_k,f.T_k_zk);
    }	
	
    for (i4=0; i4<n; i4++) { /* Check the tolerance*/
	Ftem = Ftem+fabsf(F[i4+1]);
	if (Ftem<n*tol && i4 == n-1) {
	    for (j1=0; j1<n; j1++) {
		if (debug) {
		    sf_warning("F(%d) is sufficeintly close to zero. y[%d] = %g \n",j1+1,j1+1,xx[j1+1]);
		}
	    }
	    goto mark; /* Exit the loop to the part for writing the result*/
	}	
    }
	
    /* MAIN LOOP through the output for repeating yk=yk-dk-----------------------------------------------*/
	
    for (q=0; q<w; q++) {
	Ftem=0; /* Reset Ftem to zero*/
		
	for (i2=0; i2<n; i2++) { /* Recalculate F for new y (Repeat Step 1)*/
	    initialize(i2+1,n,xx,v,xref,zref,gx,gz,aniso,z,zder,zder2);
	    F[i2+1] = T_hat_1k_k(f.T_k_k1,f.T_k_zk1) + T_hat_k_k(f.T_k_k,f.T_k_zk);
	    
	    /*Zero thickness layer case*/
	    if (fabsf(z(i2+1,xx[i2+1])-z(i2+2,xx[i2+2]))<0.00001 || fabsf(z(i2,xx[i2])-z(i2+1,xx[i2+1]))<0.00001){ 
		/*Cheap trick to get it to exit. Not respect Snell's law. Traveltime derivative is 0/0*/
		F[i2+1] = 0.00;
	    }
	}
		
	for (i5=0; i5<n; i5++) { /* Check the tolerance*/
	    Ftem = Ftem+fabsf(F[i5+1]);
	    if (Ftem<n*tol && i5 == n-1) {
		for (j2=0; j2<n; j2++) {
		    if (debug) {
			sf_warning("F(%d) is sufficeintly close to zero. y[%d] = %g \n",j2+1,j2+1,xx[j2+1]);
		    }
		}
		goto mark; /* Exit the loop to the part for writing the result*/
	    }
	}
		
	/* Step 2: Forward recursion-------------------------------------------------------------------------*/
		
	for (l=0; l<n; l++) {
	    initialize(l+1,n,xx,v,xref,zref,gx,gz,aniso,z,zder,zder2);
	    if (l==0) {
		ck_inv[1]= 1/(T_hat_1k_k_k(f.T_k_k1_k1,f.T_k_k1_zk1,f.T_k_zk1,f.T_k_zk1_zk1) + T_hat_k_k_k(f.T_k_k_k,f.T_k_k_zk,f.T_k_zk,f.T_k_zk_zk));
		zk[1] = T_hat_1k_k(f.T_k_k1,f.T_k_zk1) +T_hat_k_k(f.T_k_k,f.T_k_zk);
	    }
	    else {
		ck_inv[l+1]= 1/(T_hat_1k_k_k(f.T_k_k1_k1,f.T_k_k1_zk1,f.T_k_zk1,f.T_k_zk1_zk1) + T_hat_k_k_k(f.T_k_k_k,f.T_k_k_zk,f.T_k_zk,f.T_k_zk_zk) - T_hat_1k_1k_k(f.T_k_k_k1,f.T_k_k1_zk,f.T_k_k_zk1,f.T_k_zk_zk1)*ck_inv[l]*T_hat_1k_1k_k(f.T_k_k_k1,f.T_k_k1_zk,f.T_k_k_zk1,f.T_k_zk_zk1));
		zk[l+1] = T_hat_1k_k(f.T_k_k1,f.T_k_zk1) + T_hat_k_k(f.T_k_k,f.T_k_zk) - T_hat_1k_1k_k(f.T_k_k_k1,f.T_k_k1_zk,f.T_k_k_zk1,f.T_k_zk_zk1)*ck_inv[l]*zk[l];
	    }
			
	    if (isnan(1/ck_inv[l+1]) != 0 || !isfinite(1/ck_inv[l+1]) != 0) {
		sf_warning("ck_inv doesn't exist. The solutions do not converge.\n");
		b3 = -1;
		goto mark;
	    }
	}	
		
	/* Step 3: Backward recursion-----------------------------------------------------------------------*/
		
	for (m=n-1; m>=0; m--) { 
	    initialize(m+1,n,xx,v,xref,zref,gx,gz,aniso,z,zder,zder2);
	    if (m==n-1) {
		dk[m+1] = ck_inv[m+1]*zk[m+1];
	    }
	    else {
		dk[m+1] = ck_inv[m+1]*(zk[m+1]-T_hat_k_k_k1(f.T_k_k_k1,f.T_k_k1_zk,f.T_k_k_zk1,f.T_k_zk_zk1)*dk[m+2]);
	    }
	}
		
	/*Apply boundary---------------------------------------------------------------------------*/
	dk[0] = 0.0; /*To prevent trash data*/
	dk[n+1] = 0.0;
	xxtem[0] = xx[0]; /*Fixed source*/
	xxtem[n+1] = xx[n+1]; /*Fixed receiver*/
	
	for (a=0; a<n; a++) {
	    b1=0;
	    b2=0;
	    xxtem[a+1] = xx[a+1]-dk[a+1];
	    while (xxtem[a+1]<bmin-tol && b1<20) {/* Maximum times to multiply is 20*/
				
		dk[a+1]=0.5*dk[a+1]; /* Decrease the change by half*/
		if (debug) {
		    sf_warning("The new x value exceeds the minimum boundary. dk[%d] is reduced to %g\n",a+1,dk[a+1]);
		}
		xxtem[a+1] = xx[a+1]-dk[a+1]; /* Recompute xxtem to see if it is still exceed the boundary*/
		b1++;
	    }
	    while(xxtem[a+1]>bmax+tol && b2<20) {/* Maximum times to multiply is 20*/
				
		dk[a+1]=0.5*dk[a+1];
		if (debug) {
		    sf_warning("The new x value exceeds the maximum boundary. dk[%d] is reduced to %g\n",a+1,dk[a+1]);
		}
		xxtem[a+1] = xx[a+1]-dk[a+1];
		b2++;
	    }
	    if (b1>=20) {
		if (debug) sf_warning("The ray can't be traced at Reflector %d s=%g and r=%g", n+1,xx[0],xx[n+1]);
		b3 = -1;
		goto mark;
	    }
	    if (b2>=20) {
		if (debug) sf_warning("The ray can't be traced at Reflector %d s=%g and r=%g", n+1,xx[0],xx[n+1]);
		b3= -1;
		goto mark;
	    }
	}
	
	/*Zero thickness---we modify the dk from newton to converge to the same point*/
	for(b4=0; b4<n+1; b4++) {
	    if (fabsf(z(b4,xxtem[b4])-z(b4+1,xxtem[b4+1]))<0.00001){
		if(fabsf(dk[b4])>fabsf(dk[b4+1])) dktemp = dk[b4+1];
		if(fabsf(dk[b4])<fabsf(dk[b4+1])) dktemp = dk[b4];
			
		if (xx[b4]<xx[b4+1]){
		    if(b4==0){ /*First layer = zero thickness*/
			dk[b4+1] = fabsf(xx[b4]-xx[b4+1])/2 +dktemp;
			dk[b4] = dktemp;
		    }
		    else if(b4==n) { /*Last layer = zero thickness*/
			dk[b4+1] = dktemp;
			dk[b4] = (-1)*fabsf(xx[b4]-xx[b4+1])/2 +dktemp;
		    }
		    else { /*Any other layer*/
			dk[b4+1] = fabsf(xx[b4]-xx[b4+1])/2 +dktemp;
			dk[b4] = (-1)*fabsf(xx[b4]-xx[b4+1])/2 +dktemp;
		    }
		}
		else {
		    if(b4==0){ /*First layer = zero thickness*/
			dk[b4+1] = (-1)*fabsf(xx[b4]-xx[b4+1])/2 +dktemp;
			dk[b4] = dktemp;
		    }
		    else if(b4==n) { /*Last layer = zero thickness*/
			dk[b4+1] = dktemp;
			dk[b4] = fabsf(xx[b4]-xx[b4+1])/2 +dktemp;
		    }
		    else { /*Any other layer*/
			dk[b4] = fabsf(xx[b4]-xx[b4+1])/2 +dktemp;
			dk[b4+1] = (-1)*fabsf(xx[b4]-xx[b4+1])/2 +dktemp;
		    }
		}
	    }
	}
	
		
	/* Step 4: Update xx------------------------------------------------------------------------------*/
		
	vector_sub(n+2,xx,n+2,dk,xxnew,0); /* Update xx (Newton)*/
	for (t=0; t<n+2;t++) {
	    if (debug) {
		sf_warning("The original value of y[%d] is %g, d[%d] is %g and the new value of y[%d] is %g\n",t,*(xx+t),t,*(dk+t),t,*(xxnew+t));
		if (t==n+1) {
		    sf_warning("Iteration:%d\n",q+1);
		}					
	    }
	    *(xx+t) = *(xxnew+t); /* Update xx values*/
	}
    }
	
    /* END OF MAIN LOOP--------------------------------------------------------------------------------*/	
	
    /* Write Compute traveltime-----------------------------------------------------*/
	
mark: /* Mark point for goto*/
	
    if (b3!=0) { /* output the data to table in the case where the ray can't be traced*/
	skip = true;
	table->t = -100; 
	table->tx = 0;
	table->ty = 0;
	table->tn = 0;
	table->a = 0;
	table->an = 0;
	table->ar = 0;
		
	b3 = 0; /* Return value to 0*/
    }
    else {
	skip = false;
	for (c2=0; c2<n+2; c2++) {
			
	    for (c1=0; c1<2; c1++) {
		if (c1==0) {
		    ans[c2][c1] = xx[c2];
		}
		else {
		    ans[c2][c1] = z(c2,xx[c2]);
		}
	    }
	}
	
	if (q == niter){
	    for (c5=0; c5<n; c5++) {
		if (debug) sf_warning("F(%d) is sufficiently close to zero. y[%d] = %g %d\n",c5+1,c5+1,xx[c5+1]);
	    }
	}
	
	/* To export the old ans for fwdxini-------------------------------------------*/
		
	if (NULL != oldans) {
	    for (c4=0; c4<n; c4++) {
		oldans[c4][n-1] = xx[c4+1];
	    }
	}
		
	tt=0; /* Initialize traveltime tt*/
	at=1; /* Initializa at*/
		
	for (c=0; c<n+1; c++) {
	    half_initialize(c,n,xx,v,xref,zref,gx,gz,aniso,z,zder,zder2);
	    tt = tt + T_hat_k(f.T_k);
			
	    if (c==0) {
		at = fabsf(T_hat_k_k_k1(f.T_k_k_k1,f.T_k_k1_zk,f.T_k_k_zk1, f.T_k_zk_zk1)); 
		/* tx_s = T_hat_k_k(f.T_k_k,f.T_k_zk); */ /* x-direction on the surface*/
		tz_s = T_hat_k(f.T_k_zk);
	    }
			
			
	    if (c==n) {
		
			if(vstatus != 2){
			    v_r = v[c]+gx[c]*(xx[c+1]-xref[c])+gz[c]*(z(c+1,xx[c+1])-zref[n]);
			    v_1r = v[c]+gx[c]*(xx[c]-xref[c])+gz[c]*(z(c,xx[c])-zref[n]);
			}
			else { /*anisotropy case*/
			    S1 = (aniso[c][0]*(aniso[c][0]-aniso[c][1])*pow((aniso[c][2]-1),2)*(aniso[c][3]-1))/(2*(aniso[c][1]*aniso[c][1]*(aniso[c][3]-1)*(aniso[c][2]-aniso[c][3])+aniso[c][0]*aniso[c][0]*(aniso[c][2]-1)*(aniso[c][2]*aniso[c][2]*(aniso[c][2]-1)-aniso[c][3]+1) +aniso[c][0]*aniso[c][1]*(1-aniso[c][3]*(aniso[c][2]*pow((aniso[c][2]-1),2) -aniso[c][3] +2)) ));
			    S3 = (aniso[c][1]*(aniso[c][1]-aniso[c][0])*pow((aniso[c][3]-1),2)*(aniso[c][2]-1))/(2*(aniso[c][0]*aniso[c][0]*(aniso[c][2]-1)*(aniso[c][3]-aniso[c][2])+aniso[c][1]*aniso[c][1]*(aniso[c][3]-1)*(aniso[c][3]*aniso[c][3]*(aniso[c][3]-1)-aniso[c][2]+1) +aniso[c][1]*aniso[c][0]*(1-aniso[c][2]*(aniso[c][3]*pow((aniso[c][3]-1),2) -aniso[c][2] +2)) ));
			    v_r = 	1/sqrt(((pow(-xx[c] + xx[c+1],2)/aniso[c][0] + pow(-z(c,xx[c]) + z(c+1,xx[c+1]),2)/aniso[c][1])*(1 - 
																	 (S1*pow(-xx[c] + xx[c+1],2) + S3*pow(-z(c,xx[c]) + z(c+1,xx[c+1]),2))/(pow(-xx[c] + xx[c+1],2) + pow(-z(c,xx[c]) + z(c+1,xx[c+1]),2))))/(pow(-xx[c] + xx[c+1],2) + pow(-z(c,xx[c]) + z(c+1,xx[c+1]),2)) + 
					       ((S1*pow(-xx[c] + xx[c+1],2) + S3*pow(-z(c,xx[c]) + z(c+1,xx[c+1]),2))*sqrt(pow(pow(-xx[c] + xx[c+1],2)/aniso[c][0] + pow(-z(c,xx[c]) + z(c+1,xx[c+1]),2)/aniso[c][1],2)/
															   pow(pow(-xx[c] + xx[c+1],2) + pow(-z(c,xx[c]) + z(c+1,xx[c+1]),2),2) + 
															   (2*pow(-xx[c] + xx[c+1],2)*pow(-z(c,xx[c]) + z(c+1,xx[c+1]),2)*(-1 + (aniso[c][2]*pow(-xx[c] + xx[c+1],2) + aniso[c][3]*pow(-z(c,xx[c]) + z(c+1,xx[c+1]),2))/(pow(-xx[c] + xx[c+1],2) + pow(-z(c,xx[c]) + z(c+1,xx[c+1]),2))))/
															   (aniso[c][0]*aniso[c][1]*(pow(-xx[c] + xx[c+1],2) + pow(-z(c,xx[c]) + z(c+1,xx[c+1]),2))*(S1*pow(-xx[c] + xx[c+1],2) + S3*pow(-z(c,xx[c]) + z(c+1,xx[c+1]),2)))))/
					       (pow(-xx[c] + xx[c+1],2) + pow(-z(c,xx[c]) + z(c+1,xx[c+1]),2)));
			    v_1r = v_r;
			}
			tx_r = T_hat_k_k1(f.T_k_k1,f.T_k_zk1); /* x-direction at the reflector*/
			ty = 0; /* For 3D or 2.5D & at the reflector*/
			tz_r = T_hat_k(f.T_k_zk1); /* z-direction  need to return T_k_zk >> use T_hat_k which return the function itself & at the reflector*/
			tn = ((-1)*T_hat_k(f.T_k_k1)*zder(c+1,xx[c+1])+tz_r)/hypotf(1,zder(c+1,xx[c+1])); /* Normal-direction at the reflector*/
			dip = atan(zder(c+1,xx[c+1])); /* The dipping angle of the reflector*/
					
			if (debug) {
			    sf_warning("Traveltime is %g and the total number of iterations is %d where s=%g r=%g \n\n",tt,q,xx[0],xx[n+1]);
			}
	    }
	}
		
	for (c3=0; c3<n; c3++) {
	    initialize(c3+1,n,xx,v,xref,zref,gx,gz,aniso,z,zder,zder2);
			
	    if (c3==0) {
		ck_in= 1/(T_hat_1k_k_k(f.T_k_k1_k1,f.T_k_k1_zk1,f.T_k_zk1,f.T_k_zk1_zk1) + T_hat_k_k_k(f.T_k_k_k,f.T_k_k_zk,f.T_k_zk,f.T_k_zk_zk));
	    }
	    else {
		ck_in= 1/(T_hat_1k_k_k(f.T_k_k1_k1,f.T_k_k1_zk1,f.T_k_zk1,f.T_k_zk1_zk1) + T_hat_k_k_k(f.T_k_k_k,f.T_k_k_zk,f.T_k_zk,f.T_k_zk_zk) - T_hat_1k_1k_k(f.T_k_k_k1,f.T_k_k1_zk,f.T_k_k_zk1,f.T_k_zk_zk1)*ck_in_temp*T_hat_1k_1k_k(f.T_k_k_k1,f.T_k_k1_zk,f.T_k_k_zk1,f.T_k_zk_zk1));
	    }
			
	    at = at*fabsf(ck_in*T_hat_k_k_k1(f.T_k_k_k1,f.T_k_k1_zk,f.T_k_k_zk1, f.T_k_zk_zk1));
	    ck_in_temp = ck_in;
	}
		
	table->t = tt; /* output the data to table */
/*	sf_warning("Traveltime=%g s=%g r=%g refl=%d\n",tt,xx[0],xx[n+1],n+1);*/
	table->tx = tx_r;
	table->ty = ty;
	table->tn = sqrt(fabsf(1/(v_r*v_r)-tx_r*tx_r)); /* To make it compatible with the convention of the old kirmod (Actually in z-direction)*/
		
	alpha = tz_s*v_1r; /* Cosine of the angle from vertical of the ray on the surface */
	theta = tn*v_r; /* Cosine of the angle from the normal at the reflection point of the last reflection */
		
	if (alpha < 0) alpha *= -1;
		
	if(theta < 0) theta *= -1;
		
	table->a = pow(-1,n)*sqrt(alpha*theta/(at*cos(dip)))/v_r; /* Amplitude = 1/sqrt(geometrical spreading) but we want the reciprocal of amplitude */
	an = -acosf((tz_r-xx[n+1]*zder(n+1,xx[n+1])/(v_r*hypotf(xx[n]-xx[n+1],z(n,xx[n])-z(n+1,xx[n+1]))))*v_r*cos(dip)); /* Angle from vertical*/
	if (isnan(an) != 0) { /*Just in case the argument of acosf exceeds the proper boundary due to error in computation*/
	    an = 0;
	}
	table->an = an;
	table->ar = 0.5;
	
    }
	
}

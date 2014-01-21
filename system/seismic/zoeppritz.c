/* Zoepritz equations */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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

/* 
   G.B. Young and L.W. Braile: A computer program for the application
   of Zoeppritz's amplitude equations and Knott's energy equations,
   Bull. Seis. Soc. Am., Vol. 66, No. 6, pp. 1881-1885, 1976
*/

void zoeppritz (int icoef /* (1, 2, 3) 
			     1: particle displacement
			     2: displacement potential
			     3: energy
			  */,
		float vp1, float vp2, 
		float vs1, float vs2, 
		float rho1, float rho2 /* P and S velocities and densitities */, 
		bool incp /* true if incident P wave, false if S wave */,
		float theta /* incident ray parameter */,
		float *rc  /* [4] reflection and transmission coefficients */,
		float *ang /* [4] phase rotation angle in radians */)
/*< reflection coefficients from Zoeppritz equation >*/
{
    int j;
    float fac1,fac2,fac3,fac4,thetasq,qa,t1,t2,t3,t4;
    float a1,b1,a2,b2,a3,b3,a4,b4,x,y,z,fe[4],fp[4];
    sf_complex p1,p2,p3,p4, d, r[4];

    thetasq = theta * theta;
    qa = 2.0 * (rho2*vs2*vs2 - rho1*vs1*vs1);

    t1 = vp1*vp1 * thetasq;
    t2 = vs1*vs1 * thetasq;
    t3 = vp2*vp2 * thetasq;
    t4 = vs2*vs2 * thetasq;

    /* test for critical P reflection */
    if (theta > 1.0/vp1) {
	b1 = -sqrtf(t1-1.0);
	a1 = 0.0;
	fac1 = 0.0;
    } else {
	a1 = sqrt(1.0-t1);
	b1 = 0.0;
	fac1 = 1.0;
    }

    /* test for critical S reflection */
    if (theta > 1.0/vs1) {
	b2 = -sqrtf(t2-1.0);
	a2 = 0.0;
	fac2 = 0.0;
    } else {
	a2 = sqrtf(1.0-t2);
	b2 = 0.0;
	fac2 = 1.0;
    }

    /* test for critical P refraction */
    if (theta > 1.0/vp2) {
	b3 = -sqrtf(t3-1.0);
	a3 = 0.0;
	fac3 = 0.0;
    } else {
	a3 = sqrtf(1.0-t3);
	b3 = 0.0;
	fac3 = 1.0;
    }

    /* test for critical S refraction */
    if (theta > 1.0/vs2) {
	b4 = -sqrtf(t4-1.0);
	a4 = 0.0;
	fac4 = 0.0;
    } else {
	a4 = sqrtf(1.0-t4);
	b4 = 0.0;
	fac4 = 1.0;
    }

    x = rho2 - qa*thetasq;
    y = rho1 + qa*thetasq;
    z  = rho2 - rho1 - qa*thetasq;
/*    z1 = rho1 - rho2 + qa*thetasq; */

    p1 = sf_cmplx(a1,b1);
    p2 = sf_cmplx(a2,b2);
    p3 = sf_cmplx(a3,b3);
    p4 = sf_cmplx(a4,b4);

#ifdef SF_HAS_COMPLEX_H	 
    d = vp1*vp2*vs1*vs2*thetasq*z*z + vp2*vs2*p1*p2*x*x + 
	vp1*vs1*p3*p4*y*y + rho1*rho2*(vs1*vp2*p1*p4+vp1*vs2*p2*p3) + 
	qa*qa*thetasq*p1*p2*p3*p4;
#else
    d = sf_cmplx(0.,0.);
    sf_error("No complex support yet");
#endif

    /* compute the coefficients */
    if (incp) {
#ifdef SF_HAS_COMPLEX_H	 
	r[0] = -1.0 + 2.0*p1*(vp2*vs2*p2*x*x + vs1*vp2*rho1*rho2*p4 + qa*qa*thetasq*p2*p3*p4)/d;
	r[1] = -2.0*vp1*theta*p1*(qa*p3*p4*y + vp2*vs2*x*z)*fac2/d;
	r[2] =  2.0*vp1*rho1*p1*(vs2*p2*x + vs1*p4*y)*fac3/d;
	r[3] = -2.0*vp1*rho1*theta*p1*(qa*p2*p3 - vs1*vp2*z)*fac4/d;
#else
	sf_error("No complex support yet");
#endif

	if (icoef > 1) {
	    fp[0] = 1.0;
	    fp[1] = vs1 / vp1;
	    fp[2] = vp2 / vp1;
	    fp[3] = vs2 / vp1;
	    if (icoef > 2) {
		fe[0] = 1.0;
#ifdef SF_HAS_COMPLEX_H	 
		fe[1] = p2*vp1 / (p1*vs1);
		fe[2] = rho2*p3*vp1 / (p1*vp2*rho1);
		fe[3] = rho2*p4*vp1 / (rho1*p1*vs2);
#else
		sf_error("No complex support yet");
#endif
	    }
	}
    } else {
#ifdef SF_HAS_COMPLEX_H	 
	r[0] = -2.0*vs1*theta*p2*(qa*p3*p4*y + vp2*vs2*x*z)*fac1/d;
	r[1] = 1.0 - 2.0*p2*(vp2*vs2*p1*x*x + vp1*vs2*rho1*rho2*p3 + qa*qa*thetasq*p1*p3*p4)/d;
	r[2] = 2.0*vs1*rho1*theta*p2*(qa*p1*p4 - vp1*vs2*z)*fac3/d;
	r[3] = 2.0*vs1*rho1*p2*(vp1*p3*y + vp2*p1*x)*fac4/d;
#else
	sf_error("No complex support yet");
#endif

	if (icoef > 1) {
	    fp[0] = vp1 / vs1;
	    fp[1] = 1.0;
	    fp[2] = vp2 / vs1;
	    fp[3] = vs2 / vs1;
	    if (icoef > 2) {
#ifdef SF_HAS_COMPLEX_H	 
		fe[0] = vs1*p1 / (vp1*p2);
		fe[1] = 1.0;
		fe[2] = rho2*vs1*p3 / (rho1*vp2*p2);
		fe[3] = rho2*vs1*p4 / (rho1*vs2*p2);
#else
		sf_error("No complex support yet");
#endif
	    }
	}
    }
    
    for (j=0; j < 4; j++) {
	switch (icoef) {
	    case 1:
		rc[j] = cabsf(r[j]);
		break;
	    case 2:
#ifdef SF_HAS_COMPLEX_H
		rc[j] = cabsf(r[j]*fp[j]);
#else
		rc[j] = cabsf(sf_crmul(r[j],fp[j]));
#endif
		break;
	    case 3:
#ifdef SF_HAS_COMPLEX_H
		rc[j] = crealf(r[j]*fp[j]*conjf(r[j]*fp[j]))*fe[j];
#else
		rc[j] = crealf(sf_cmul(sf_crmul(r[j],fp[j]),
				       conjf(sf_crmul(r[j],fp[j]))))*fe[j];
#endif
		break;
	    case 4:
		rc[j] = crealf(r[j]);
		break;
	    default:
		sf_error("%s: wrong icoef",__FILE__);
		break;
	}

	ang[j] = cargf(r[j]);
    }
}


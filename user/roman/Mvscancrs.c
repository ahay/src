/* Velocity analysis.

Inverse of sfvelmod
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
#include "fint1.h"

static float v, x, h, v1, s;
static float i2V0, xx0, xx0sq, b, cs2, sn, kcre, kcee, h2; // a1x, a2x2;

void beta_cee_cre_2_a(float t0, float cs2, float sn, float rcre, float rcee, 
		      float * a1, float * a2, float * a3)
{
    const float tmp = i2V0 * t0 * cs2;

    *a1 = sn * i2V0;
    *a2 = tmp * rcee;
    *a3 = tmp * rcre;
}
static float crsrk(float t, int it) 
{
    float a1, a2, a3, a1x, a2x2, a3h2;

    beta_cee_cre_2_a(t, cs2, sn, kcre, kcee, &a1, &a2, &a3);		

    a1x = a1*xx0;
    a2x2 = a2*xx0sq;
    a3h2 = a3*h2; 

    return sqrtf( (t + a1x)*(t + a1x) + a2x2 + a3h2);
}

static float hyperb(float t, int it) 
{
    return hypotf(t,v);
}

static float nonhyperb(float t, int it) 
{
    /* shifted hyperbola */
    return t*(1.0-1.0/s) + sqrtf(t*t+s*v*v)/s;
}

static float hyperb1(float t, int it) 
{
    return sqrtf(t*t+v*v-v1*v1*h*h);
}

static float nonhyperb1(float t, int it) 
{
    return t*(1.0-1.0/s) + sqrtf(t*t+s*(v*v-v1*v1*h*h))/s;
}

static float curved(float t, int it) 
{
    return sqrtf(t*t+v*h);
}

static float noncurved(float t, int it) 
{
    return t*(1.0-1.0/s) + sqrtf(t*t+s*v*h)/s;
}

static float curved1(float t, int it) 
{
    return sqrtf(t*t+v*h-v1*h*h);
}

static float noncurved1(float t, int it) 
{
    return t*(1.0-1.0/s) + sqrtf(t*t+s*(v*h-v1*h*h))/s;
}

int main(int argc, char* argv[])
{
    fint1 crs/*nmo*/;
    bool sembl, half, slow, dsembl, asembl, weight, squared;
    int it,ih,ix,nt,nh,nx,ib,ie,nb,i, nw, CDPtype, mute, *mask=NULL;
    float amp, amp2, dt, dh, t0, h0, smax, num, den, dy, str, sh=0., sh2=0.;
    float *trace=NULL, *** *stack=NULL, *** *stack2=NULL, *** *stackh=NULL, *hh=NULL;
    char *time=NULL, *space=NULL, *unit=NULL;
    size_t len;
    sf_file cmp=NULL, scan=NULL, offset=NULL, msk=NULL;
    mapfunc crsfunc/*nmofunc*/;

    float iNth, x0, dx, x_central; // a10, a20, a30, da1, da2, da3 ;
    float V0, b0, db, rcre0, drcre, rcee0, drcee, rcee, rcre;
    //int   na1, na2, na3, ia1, ia2, ia3, Nta1a2a3;
    int   nrcre, nrcee, ia1, ia2, ia3;

    sf_init (argc,argv);
    cmp = sf_input("in");
    scan = sf_output("out");

    if (!sf_histint(cmp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(cmp,"n2",&nh)) sf_error("No n2= in input");
    if (!sf_histint(cmp,"n3",&nx)) sf_error("No n3= in input");
    // nx = sf_leftsize(cmp,2);

    if (!sf_getbool("semblance",&sembl)) sembl=false;
    /* if y, compute semblance; if n, stack */
    if (sembl || !sf_getbool("diffsemblance",&dsembl)) dsembl=false;
    /* if y, compute differential semblance */
    if (sembl || dsembl || !sf_getbool("avosemblance",&asembl)) asembl=false;
    /* if y, compute AVO-friendly semblance */
    if (!sf_getint("nb",&nb)) nb=2;
    /* semblance averaging */
    if (!sf_getbool("weight",&weight)) weight=true;
    /* if y, apply pseudo-unitary weighting */

    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");


    if (!sf_histfloat(cmp,"o3",&x0)) sf_error("No o3= in input");
    if (!sf_histfloat(cmp,"d3",&dx)) sf_error("No d3= in input");
    x_central = x0 + (nx-1)/2.f*dx;


    if (!sf_getbool("half",&half)) half=true;
    /* if y, the second axis is half-offset instead of full offset */

    CDPtype=1;
    if (NULL != sf_getstring("offset")) {
	offset = sf_input("offset");
	hh = sf_floatalloc(nh);

	h0 = dh = 0.;
    } else {
	if (!sf_histfloat(cmp,"o2",&h0)) sf_error("No o2= in input");
	if (!sf_histfloat(cmp,"d2",&dh)) sf_error("No d2= in input");
	
	sf_putfloat(scan,"h0",h0);
	sf_putfloat(scan,"dh",dh);
	sf_putint(scan,"nh",nh);

	if (sf_histfloat(cmp,"d3",&dy)) {
	    CDPtype=half? 0.5+dh/dy: 0.5+0.5*dh/dy;
	    if (0 == CDPtype) CDPtype=1;
	    if (1 != CDPtype) {
		sf_histint(cmp,"CDPtype",&CDPtype);
		sf_warning("CDPtype=%d",CDPtype);
	    }
	}

	offset = NULL;
	hh = NULL;
    }

    if (NULL != sf_getstring("mask")) {
	msk = sf_input("mask");
	mask = sf_intalloc(nh);
    } else {
	msk = NULL;
	mask = NULL;
    }
    
    /* b(betha), Rn(Rcee),Rnip(Rcre) => a2, a3 CRS params */
    if (!sf_getfloat("V0",&V0)) 
	sf_error("Need V0=");

    i2V0 = 2.f / V0;

    if (!sf_getfloat("b0",&b0)) 
	sf_error("Need b0=");
    if (!sf_getfloat("db",&db))
	sf_error("Need da1=");
    if (!sf_getint("nb",&nb)) 
	sf_error("Need nb=");
    b0 *= SF_PI / 180.f;
    db *= SF_PI / 180.f;

    if (!sf_getfloat("rcee0",&rcee0)) 
	sf_error("Need rcee0=");
    if (!sf_getfloat("drcee",&drcee))
	sf_error("Need drcee=");
    if (!sf_getint("nrcee",&nrcee)) 
	sf_error("Need nrcee=");

    if (!sf_getfloat("rcre0",&rcre0)) 
	sf_error("Need rcre0=");
    if (!sf_getfloat("drcre",&drcre))
	sf_error("Need drcre=");
    if (!sf_getint("nrcre",&nrcre)) 
	sf_error("Need nrcre=");



	/* SCAN */
    sf_putfloat(scan,"o2",rcre0);
    sf_putfloat(scan,"d2",drcre);
    sf_putint(scan,"n2",nrcre);

    sf_putfloat(scan,"o4",rcee0);
    sf_putfloat(scan,"d4",drcee);
    sf_putint(scan,"n4",nrcee);

    sf_putfloat(scan,"o3",b0);
    sf_putfloat(scan,"d3",db);
    sf_putint(scan,"n3",nb); 



    if (!sf_getfloat("smax",&smax)) smax=2.0;
    /* maximum heterogeneity */
    /*if (!sf_getint("ns",&ns)) ns=1;
     number of heterogeneity scans  
    ds = ns>1? (smax-1.0)/(ns-1): 0.;

    
    if (ns > 1) {
	sf_putfloat(scan,"o3",1.0);
	sf_putfloat(scan,"d3",ds);
	sf_putint(scan,"n3",ns);

	sf_shiftdim(cmp, scan, 3);
	} */


    if (!sf_getbool("slowness",&slow)) slow=false;
    /* if y, use slowness instead of velocity */
    sf_putstring(scan,"label2","Rcre");
    sf_putstring(scan,"label3","Betha");
    sf_putstring(scan,"label4","Rcee");

    if (!sf_getbool("squared",&squared)) squared=false;
    /* if y, the slowness or velocity is squared */

    crsfunc = crsrk; /*  hyperb; */

    /*( v1 reference velocity )
    if (!sf_getfloat("v1",&v1)) {
	if (ns > 1) {
	    nmofunc = squared? noncurved: nonhyperb;
	} else {
	    nmofunc = squared? curved: hyperb;
	}
    } else {
	if (ns > 1) {
	    nmofunc = squared? noncurved1: nonhyperb1;
	} else {
	    nmofunc = squared? curved1: hyperb1;
	}
	if (!slow) v1 = 1./v1;
    }
    */

    if (NULL != (time = sf_histstring(cmp,"unit1")) &&
	NULL != (space = sf_histstring(cmp,"unit2"))) {
	len = strlen(time)+strlen(space)+2;
	unit = sf_charalloc(len);
	if (slow) {
	    snprintf(unit,len,"%s/%s",time,space);
	} else {
	    snprintf(unit,len,"%s/%s",space,time);
	}
	sf_putstring(scan,"unit2",unit);
    }

    /* stack =  sf_floatalloc3(nt,nv,ns); 
    stack2 = (sembl || dsembl || asembl)? sf_floatalloc3(nt,nv,ns) : NULL;
    stackh = asembl? sf_floatalloc3(nt,nv,ns) : NULL; */
    stack = sf_floatalloc4(nt, nrcre, nb, nrcee);
    stack2 = (sembl || dsembl || asembl)? sf_floatalloc4(nt, nrcre, nb, nrcee) : NULL;
    stackh = asembl? sf_floatalloc4(nt, nrcre, nb, nrcee) : NULL; 

    if (!sf_getint("extend",&nw)) nw=4;
    /* trace extension */

    if (!sf_getint("mute",&mute)) mute=12;
    /* mute zone */

    if (!sf_getfloat("str",&str)) str=0.5;
    /* maximum stretch allowed */

    trace = sf_floatalloc(nt);
    crs/*nmo*/ = fint1_init(nw,nt,mute);

    iNth = 1.f/(nt*nh);

    for (ix=0; ix < nx; ix++) {

	x = x0 + ix * dx - x_central;
	xx0 = x - x_central;
	xx0sq = xx0*xx0;

	sf_warning("cmp %d of %d x=%f",ix+1,nx,x);
	
	for (it=0; it < nt*nrcre*nrcee*nb; it++) {
	    stack[0][0][0][it] = 0.;
	    if (sembl || asembl) stack2[0][0][0][it] = 0.;
	    if (asembl) stackh[0][0][0][it] = 0.;
	}

	if (NULL != offset) 
	    sf_floatread(hh,nh,offset);
	if (NULL != msk) 
	    sf_intread(mask,nh,msk);

	if (asembl) sh = sh2 = 0.;

	for (ih=0; ih < nh; ih++) {

	    sf_floatread(trace,nt,cmp); 

	    if (NULL != msk && 0==mask[ih]) continue;

	    h = (NULL != offset)? hh[ih]: 
		h0 + ih * dh + (dh/CDPtype)*(ix%CDPtype);
	    h2 = h*h;
	    if (half) h *= 2.;

	    if (asembl) {
		sh  += h;
		sh2 += h*h;
	    }

	    for (it=0; it < nt; it++) {
		trace[it] *= iNth; // /= nt*nh;
	    }
	    
	    fint1_set(crs/*nmo*/,trace);

	    /*for (is=0; is < ns; is++) {
		s = 1.0 + is*ds;

		for (iv=0; iv < nv; iv++) {
		    v = v0 + iv * dv;
		    v = slow? h*v: h/v;*/

	    for (ia3=0; ia3 < nrcee/*na2*/; ia3++) {
		rcee = rcee0 + ia3 * drcee; // da2;
		kcee = 1.f / rcee;

	    for (ia2=0; ia2 < nb/*na1*/; ia2++) {
		b = b0 + ia2 * db;//da1;
		sn = sinf(b);
		cs2 = cosf(b);
		cs2 *= cs2;

	    for (ia1=0; ia1 < nrcre/*na3*/; ia1++) {
		rcre = rcre0 + ia1 * drcre;		    
		kcre = 1.f / rcre;
		/* GloBAL ZEVEL !!! 
		   a3 *= h;
		   v = a3; */


		    stretch(crs/*nmo*/,crsfunc/*nmofunc*/,nt,dt,t0,nt,dt,t0,trace,str);

		    for (it=0; it < nt; it++) {
			
			const float a3 = i2V0 * cs2 * kcre * (t0 + it*dt);

			amp = weight? fabsf(a3)*trace[it]: trace[it];

			if (dsembl) {

			    if (ih > 0) {
				// amp2 = amp - stack2[is][iv][it];
				amp2 = amp - stack2[ia3][ia2][ia1][it];
				//stack[is][iv][it] += amp2*amp2;
				stack[ia3][ia2][ia1][it] += amp2*amp2;
			    }
			    stack[ia3][ia2][ia1][it] = amp;

			} else {
			    if (sembl || asembl) stack2[ia3][ia2][ia1][it]  += amp*amp;
			    if (asembl) stackh[ia3][ia2][ia1][it] += amp*h;

			    stack[ia3][ia2][ia1][it] += amp;
			}
		    }

		} /* a3 */
		} /* a2 */
		} /* a1 */

		/* } s */
	} /* h */
    } /* x */
    
	if (sembl || asembl) {
	    int ia1, ia2, ia3, na3 = nrcee, na2=nb, na1=nrcre;
	
	    /*for (is=0; is < ns; is++) {
	      for (iv=0; iv < nv; iv++) {*/
	    for (ia3=0; ia3 < na3; ia3++) {
	    for (ia2=0; ia2 < na2; ia2++) {
	    for (ia1=0; ia1 < na1; ia1++) {

		    for (it=0; it < nt; it++) {

			ib = it-nb;

			ie = it+nb+1;

			if (ib < 0) ib=0;

			if (ie > nt) ie=nt;

			num = 0.;
			den = 0.;

			for (i=ib; i < ie; i++) {
			    num += stack[ia3][ia2][ia1][i]*stack[ia3][ia2][ia1][i];
			    den += stack2[ia3][ia2][ia1][i];

			    /* (h2 s^2 - 2 h s sh + N sh^2)/((-h^2 + h2 N) s2) */
			    if (asembl) 
				num += (nh*stackh[ia3][ia2][ia1][i]*stackh[ia3][ia2][ia1][i] - 
					2.*sh*stack[ia3][ia2][ia1][i]*stackh[ia3][ia2][ia1][i])/sh2;
			}
			if (asembl) {
			    den *= (nh-sh*sh/sh2);
			} else {
			    den *= nh;
			}
			trace[it] = (den > 0.)? num/den: 0.;
		    }
		    sf_floatwrite(trace,nt,scan);

		} /* a1 */
		} /* a2 */
	    } /* a3 */
	} else {
	    sf_floatwrite (stack[0][0][0],nt*nrcre*nrcee*nb,scan);
	}

    exit(0);
}

/* Velocity analysis using similarity-weighted semblance.*/

/*
  Copyright (C) 2013 University of Texas at Austin

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

static float v, h, v1, s;

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
    fint1 nmo;
    bool sembl, half, slow, dsembl, asembl, weight, squared, trend, ratio;
    int it,ih,ix,iv, nt,nh,nx,nv, ib,ie,nb,i, nw, is, ns, CDPtype, mute, *mask;
    float amp, amp2, dt, dh, t0, h0, v0, dv, ds, smax, num, den, dy, str, sh=0., sh2=0.;
    float *trace, ***stack, ***stack2, ***stack2h, ***stackh, *hh, **bb;
    char *time, *space, *unit;
    const char *type;
    size_t len;
    sf_file cmp, scan, offset, msk, grd,ref;
    mapfunc nmofunc;
    sf_init (argc,argv);

	/*added codes down */
		float eps=0.01f,*traceref,thr;
		float *rat1,*rat2,*rat,*one,*two;
		ref=sf_input("ref");
	/*added codes up */


    cmp = sf_input("in");
    scan = sf_output("out");

    if (!sf_histint(cmp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(cmp,"n2",&nh)) sf_error("No n2= in input");
    nx = sf_leftsize(cmp,2);

    if (!sf_getbool("semblance",&sembl)) sembl=false;
    /* if y, compute semblance; if n, stack */
    if (sembl || !sf_getbool("diffsemblance",&dsembl)) dsembl=false;
    /* if y, compute differential semblance */
    if (sembl || dsembl || !sf_getbool("avosemblance",&asembl)) asembl=false;
    /* if y, compute AVO-friendly semblance */

    if (NULL == (type = sf_getstring("type"))) {
	/* type of semblance (avo,diff,sembl,power,weighted) */
	if (asembl) {
	    type="avo";
	} else if (dsembl) {
	    type="diff";
	} else if (sembl) {
	    type="sembl";
	} else {
	    type="power";
	}
    }

    trend = (bool) ('a' == type[0] || 'w' == type[0]);
    ratio = (bool) ('p' != type[0] && 'd' != type[0]);

    if (!sf_getint("nb",&nb)) nb=2;
    /* semblance averaging */
    if (!sf_getbool("weight",&weight)) weight=true;
    /* if y, apply pseudo-unitary weighting */

    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");

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
	/* optional mask file */ 
	msk = sf_input("mask");
	mask = sf_intalloc(nh);
    } else {
	msk = NULL;
	mask = NULL;
    }

    if (!sf_getfloat("v0",&v0) && !sf_histfloat(cmp,"v0",&v0)) 
	sf_error("Need v0=");
    /*(v0 first scanned velocity )*/
    if (!sf_getfloat("dv",&dv) && !sf_histfloat(cmp,"dv",&dv)) 
	sf_error("Need dv=");
    /*(dv step in velocity )*/
    if (!sf_getint("nv",&nv) && !sf_histint(cmp,"nv",&nv)) 
	sf_error("Need nv=");
    /*(nv number of scanned velocities )*/

    sf_putfloat(scan,"o2",v0);
    sf_putfloat(scan,"d2",dv);
    sf_putint(scan,"n2",nv);

    if (!sf_getfloat("smax",&smax)) smax=2.0;
    /* maximum heterogeneity */
    if (!sf_getint("ns",&ns)) ns=1;
    /* number of heterogeneity scans */ 
    ds = ns>1? (smax-1.0)/(ns-1): 0.;

    if (ns > 1) {
	sf_putfloat(scan,"o3",1.0);
	sf_putfloat(scan,"d3",ds);
	sf_putint(scan,"n3",ns);

	sf_shiftdim(cmp, scan, 3);
    }

    if (!sf_getbool("slowness",&slow)) slow=false;
    /* if y, use slowness instead of velocity */
    sf_putstring(scan,"label2",slow? "Slowness": "Velocity");

    if (!sf_getbool("squared",&squared)) squared=false;
    /* if y, the slowness or velocity is squared */

    if (!sf_getfloat("v1",&v1)) {
	/*( v1 reference velocity )*/
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

    if (NULL != sf_getstring("grad")) {
	grd = sf_input("grad");

	bb = sf_floatalloc2(nt,nv);
	sf_floatread(bb[0],nt*nv,grd);

	sf_fileclose(grd);
    } else {
	bb = NULL;
    }

    stack =  sf_floatalloc3(nt,nv,ns);
    stack2 = ('p' != type[0])? sf_floatalloc3(nt,nv,ns): NULL;
    stackh = trend? sf_floatalloc3(nt,nv,ns): NULL;
    stack2h = ('w' == type[0])? sf_floatalloc3(nt,nv,ns): NULL;

    if (!sf_getint("extend",&nw)) nw=4;
    /* trace extension */

    if (!sf_getint("mute",&mute)) mute=12;
    /* mute zone */

    if (!sf_getfloat("str",&str)) str=0.5;
    /* maximum stretch allowed */

    trace = sf_floatalloc(nt);

	/* added codes below */
    traceref = sf_floatalloc(nt);
    one = sf_floatalloc(nt);
    two = sf_floatalloc(nt);
    rat1 = sf_floatalloc(nt);
    rat2 = sf_floatalloc(nt);
    rat = sf_floatalloc(nt);

	int n[2]; n[0]=nt;n[1]=1;	
	bool verb=false;	
	int rect[2]; rect[0]=10;rect[1]=1;
	int niter=10; 
	if(!sf_getfloat("thr",&thr)) thr=0.3;
	if(!sf_getfloat("eps",&eps)) eps=0.01;
	if(!sf_getint("rect",rect)) rect[0]=5;
	if(!sf_getint("niter",&niter)) niter=10;
    sf_divn_init(1, nt, n, rect, niter, verb);
	/* added codes up */

    nmo = fint1_init(nw,nt,mute);

    for (ix=0; ix < nx; ix++) {
    	sf_floatread(traceref,nt,ref);
	sf_warning("cmp %d of %d;",ix+1,nx);

	for (it=0; it < nt*nv*ns; it++) {
	    stack[0][0][it] = 0.;
	    if (ratio) stack2[0][0][it] = 0.;
	    if (trend) stackh[0][0][it] = 0.;
	}

	if (NULL != offset) sf_floatread(hh,nh,offset);
	if (NULL != msk) sf_intread(mask,nh,msk);

	if (trend) sh = sh2 = 0.;

	for (ih=0; ih < nh; ih++) {
	    sf_floatread(trace,nt,cmp); 
	    if (NULL != msk && 0==mask[ih]) continue;

	    h = (NULL != offset)? hh[ih]: 
		h0 + ih * dh + (dh/CDPtype)*(ix%CDPtype);
	    if (half) h *= 2.;

	    if (trend) {
		sh  += h;    /* sf  */
		sh2 += h*h;  /* sf2 */
	    }

	    for (it=0; it < nt; it++) {
		trace[it] /= nt*nh;
	    }
	    fint1_set(nmo,trace);

	    for (is=0; is < ns; is++) {
		s = 1.0 + is*ds;

		for (iv=0; iv < nv; iv++) {
		    v = v0 + iv * dv;
		    v = slow? h*v: h/v;

		    stretch(nmo,nmofunc,nt,dt,t0,nt,dt,t0,trace,str);

		/* Add codes here below */
		for(it=0;it<nt;it++){one[it]=trace[it];two[it]=traceref[it];}
		sf_divne(one,two,rat1,eps);
        sf_divne(two,one,rat2,eps);
		sf_divn_combine (rat1,rat2,rat);	
		
		for(it=0;it<nt;it++)
			{if(rat[it]< thr) trace[it]=0; 	}	//else trace[it]=trace[it]*rat[it]; else sf_warning("trace[it]=%g",it,trace[it]);
		
		/* Add codes here up */
		    for (it=0; it < nt; it++) {
			amp = weight? fabsf(v)*trace[it]: trace[it];
			if (NULL != bb) amp *= (1.0-bb[iv][it]*h);
			
			switch(type[0]) {
			    case 'd':
				if (ih > 0) {
				    amp2 = amp - stack2[is][iv][it];
				    stack[is][iv][it] += amp2*amp2;
				}
				stack2[is][iv][it] = amp;
				break;
			    case 's':
				stack2[is][iv][it] += amp*amp;
				stack[is][iv][it] += amp;
				break;
			    case 'a': 
				stackh[is][iv][it] += amp*h;    /* saf */
				stack2[is][iv][it] += amp*amp;  /* sa2 */
				stack[is][iv][it] += amp;       /* sa1 */
				break;
			    case 'w':
				stackh[is][iv][it] += amp*h;       /* saf */
				stack2h[is][iv][it] += amp*amp*h;  /* sfa2 */
				stack2[is][iv][it] += amp*amp;     /* sa2 */
				stack[is][iv][it] += amp;          /* sa1 */
				break;
			    case 'p':
			    default:
				stack[is][iv][it] += amp;
				break;				
			} 
		    } /* t */
		} /* v */
	    } /* s */
	} /* h */
	
	if (ratio) {
	    for (is=0; is < ns; is++) {
		for (iv=0; iv < nv; iv++) {
		    for (it=0; it < nt; it++) {
			ib = it-nb;
			ie = it+nb+1;
			if (ib < 0) ib=0;
			if (ie > nt) ie=nt;
			num = 0.;
			den = 0.;
			for (i=ib; i < ie; i++) {
			    switch(type[0]) {
				case 'a':
				    /* (N*saf^2 + sa1^2*sf2 - 2*sa1*saf*sf)/((N*sf2 - sf^2)*sa2) */

				    num += nh*stackh[is][iv][i]*stackh[is][iv][i] + 
					sh2*stack[is][iv][i]*stack[is][iv][i] - 
					2.*sh*stack[is][iv][i]*stackh[is][iv][i];
				    den += stack2[is][iv][i];
				    break;
				case 'w':
				    /* 4*(sa1*sfa2 - sa2*saf)*(N*saf - sa1*sf)/(N*sfa2 - sa2*sf)^2 */

				    num += 
					(stack[is][iv][i]*stack2h[is][iv][i]-
					 stack2[is][iv][i]*stackh[is][iv][i])*
					(nh*stackh[is][iv][i]-stack[is][iv][i]*sh);
				    den += 
					(nh*stack2h[is][iv][i]-stack2[is][iv][i]*sh)*
					(nh*stack2h[is][iv][i]-stack2[is][iv][i]*sh);
				    break;
				case 's':
				default:
				    num += stack[is][iv][i]*stack[is][iv][i];
				den += stack2[is][iv][i];
				break;
			    }
			}
			    
			switch(type[0]) {
			    case 'a':
				den *= (nh*sh2-sh*sh);
				break;
			    case 'w':
				num *= 4.0f;
				break;
			    case 's':
				den *= nh;
				break;
			}

			trace[it] = (den > 0.)? num/den: 0.;
		    }
		    sf_floatwrite(trace,nt,scan);
		} /* v */
	    } /* s */
	} else {
	    sf_floatwrite (stack[0][0],nt*nv*ns,scan);
	}
    } /* x */
    sf_warning(".");
    
    exit(0);
}

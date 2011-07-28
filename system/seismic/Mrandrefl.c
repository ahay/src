/* Simple synthetics with random reflectivity. */
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
#include <stdio.h>

#include <rsf.h>

static int compare_float (const void *a, const void *b)
{
    float fa = * ((float*) a);
    float fb = * ((float*) b);
 
    return (fa < fb)? -1: (fa > fb)? 1: 0;
}

static void cumsum(int n, const float *a, float *b)
{
    float t;
    int i;

    t = 0.;
    for (i = 0; i < n; i++) {
	t += a[i];
	b[i] = t;
    }
}

/* computes synthetic for Rs at Ts(sec) for times t */
static void synf(float fo, int n, float* s, int ns, 
		 const float* t, const float *ts, const float *rs)
{
    int i, is;
    float rn, tn, a;
    
    for (i=0; i < n; i++) {
	s[i]=0.;
    }
  
    for (is=0; is < ns; is++) {
	rn=rs[is];
	tn=ts[is];
	/* Ricker wavelet at tn with amplitude Rn */
	/* accumulate the sum of all reflections */ 
	for (i=0; i < n; i++) {
	    a = SF_PI*fo*(t[i]-tn);
	    a *= a;
	    s[i] += rn*(1-2*a)*expf(-a);
	}
    }
}

int main (int argc, char* argv[])
{
    int nr, nt, it;
    float t0, dt, fo[3]={20.,8.,5.}, tscale;
    float *tim, *pp, *ps, *ss, *tpp, *tss, *tps;
    float *ts, *dtpp, *tmean, *p2ss, *p2ps, *dtss, *dtps, *rs;
    const char* func;
    sf_file mod, vpvs;

    sf_init (argc,argv);
    mod = sf_output("out");

    if (!sf_getint("nr",&nr)) sf_error("Need nr=");
    /* number of reflectors */
    if (!sf_getint("n1",&nt))   nt=3501;  
    /* time length */
    sf_putint(mod,"n1",nt);
    if (!sf_getfloat("d1",&dt)) dt=0.001; 
    /* time sampling */
    sf_putfloat(mod,"d1",dt);
    if (!sf_getfloat("o1",&t0)) t0=0.0;   
    /* time origin */
    sf_putfloat(mod,"o1",t0);
    sf_putint(mod,"n2",3);
    sf_putstring(mod,"label1","Time");
    sf_putstring(mod,"unit1","s");
    sf_setformat(mod,"native_float");

    if (!sf_getfloat("tscale",&tscale)) tscale=1.;
    /* maximum time */

    sf_getfloats ("fo",fo,3); 

    vpvs = sf_output("vpvs");
    sf_putint(vpvs,"n1",nr);
    sf_putint(vpvs,"n2",2);
    sf_setformat(vpvs,"native_float");

    tim = sf_floatalloc (nt);
    pp = sf_floatalloc (nt);
    ps = sf_floatalloc (nt);
    ss = sf_floatalloc (nt);
    ts = sf_floatalloc (nr);

    for (it=0; it < nt; it++) {
	tim[it] = t0 + it*dt;
    }

    dtpp = sf_floatalloc (nr);
    tmean = sf_floatalloc (nr);
    rs = sf_floatalloc (nr);
    tpp = sf_floatalloc (nr);
    tss = sf_floatalloc (nr);
    tps = sf_floatalloc (nr);
    p2ss = sf_floatalloc (nr);
    p2ps = sf_floatalloc (nr);
    dtss = sf_floatalloc (nr);
    dtps = sf_floatalloc (nr);

    if (NULL == (func = sf_getstring("func"))) func="const";
    /* type of vpvs function */

    init_genrand(2003);

    /* ts - reflector positions */    
    for (it=0; it < nr; it++) {
	ts[it] = tscale*(0.05+0.95*genrand_real1());
    }
    qsort(ts,nr,sizeof(float),compare_float);

    /* dtpp - layer thickness in PP time */    
    dtpp[0] = ts[0];
    for (it=1; it < nr; it++) {
	dtpp[it] = ts[it] - ts[it-1];
    }

    /* tmean - time in the middle of the layer */
    /* p2ss - Vp/Vs ratio as a function of time */
    for (it=0; it < nr; it++) {
	tmean[it]=ts[it]-0.5*dtpp[it];
	switch (func[0]) {
	    case 'h': /* Hamilton */
		p2ss[it]=1./(0.12+(0.5*tmean[it]/tscale));
		break;
	    case 's': /* sinusoid */
		p2ss[it]=2.+0.2*cosf(4.*SF_PI*tmean[it]/tscale);
		break;
	    case 'c': /* constant */
		p2ss[it]=2.;
		break;
	    case 'l': /* linear */
		p2ss[it]=3.-tmean[it];
		break;
	    case 'm': /* medium */
		p2ss[it]=1.+0.5/(0.12+(0.5*tmean[it]/tscale));
		break;
	    default:
		sf_error("Unknown case %s",func);
		break;
	}
    }

    sf_floatwrite(tmean,nr,vpvs);
    sf_floatwrite(p2ss,nr,vpvs);
    sf_fileclose(vpvs);

    for (it=0; it < nr; it++) {
	p2ps[it]=0.5*(1+p2ss[it]);
	dtss[it]=p2ss[it]*dtpp[it]; /* SS thickness */
	dtps[it]=p2ps[it]*dtpp[it]; /* PS thickness */
    }

    /* rs - reflection coefficient */
    sf_randn (nr, rs);
    for (it=0; it < nr; it++) {
	rs[it] = 0.1/nr + 0.05*rs[it];
/*	fprintf(stderr,"%g ",rs[it]); */
    }
/*    fprintf(stderr,"\n"); */

    cumsum(nr,dtpp,tpp);
    cumsum(nr,dtps,tps);
    cumsum(nr,dtss,tss);

    synf(fo[0],nt,pp,nr,tim,tpp,rs);
    synf(fo[1],nt,ps,nr,tim,tps,rs);
    synf(fo[2],nt,ss,nr,tim,tss,rs);

    sf_floatwrite(pp,nt,mod);
    sf_floatwrite(ps,nt,mod);
    sf_floatwrite(ss,nt,mod);


    exit (0);
}

/* 	$Id$	 */

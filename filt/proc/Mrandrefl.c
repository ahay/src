#include <math.h>

#include <rsf.h>

#include "randn.h"

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
    float t0, dt, fo[3]={20.,8.,5.};
    float *tim, *pp, *ps, *ss, *tpp, *tss, *tps;
    float *ts, *dtpp, *tmean, *p2ss, *p2ps, *dtss, *dtps, *rs;
    sf_file mod, vpvs;

    sf_init (argc,argv);
    mod = sf_output("out");

    if (!sf_getint("nr",&nr)) sf_error("Need nr=");
    if (!sf_getint("n1",&nt))   nt=3501;  sf_putint(mod,"n1",nt);
    if (!sf_getfloat("d1",&dt)) dt=0.001; sf_putfloat(mod,"d1",dt);
    if (!sf_getfloat("o1",&t0)) t0=0.0;   sf_putfloat(mod,"o1",t0);
    sf_putint(mod,"n2",3);
    sf_setformat(mod,"native_float");

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

    /* ts - reflector positions */    
    srand(2003);

    for (it=0; it < nr; it++) {
	ts[it] = 0.1+0.9*random_one();
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
	p2ss[it]=1./(0.12+(0.5*tmean[it]));
    }

    sf_write(tmean,sizeof(float),nr,vpvs);
    sf_write(p2ss,sizeof(float),nr,vpvs);
    sf_fileclose(vpvs);

    for (it=0; it < nr; it++) {
	p2ps[it]=0.5*(1+p2ss[it]);
	dtss[it]=p2ss[it]*dtpp[it]; /* SS thickness */
	dtps[it]=p2ps[it]*dtpp[it]; /* PS thickness */
    }

    /* rs - reflection coefficient */
    randn (nr, rs);
    for (it=0; it < nr; it++) {
	rs[it] = 0.1/nr + 0.05*rs[it];
    }

    cumsum(nr,dtpp,tpp);
    cumsum(nr,dtps,tps);
    cumsum(nr,dtss,tss);

    synf(fo[0],nt,pp,nr,tim,tpp,rs);
    synf(fo[1],nt,ps,nr,tim,tps,rs);
    synf(fo[2],nt,ss,nr,tim,tss,rs);

    sf_write(pp,sizeof(float),nt,mod);
    sf_write(ps,sizeof(float),nt,mod);
    sf_write(ss,sizeof(float),nt,mod);

    exit (0);
}





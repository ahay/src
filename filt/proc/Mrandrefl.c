#include <stdlib.h>
#include <math.h>

#include <rsf.h>

#include "randn.h"

static const float pi = 3.1415926535898;

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
	    a = pi*fo*(t-tn);
	    a *= a;
	    s[i] += rn*(1-2*a)*expf(-a);
	}
    }
}

int main (int argc, char* argv[])
{
    int nr, nt, it;
    float t0, dt, fo[3]={20.,8.,5.};
    float *t, *pp, *ps, *ss, *tpp, *tss, *tps;
    float *ts, *dtpp, *tm, *p2ss, *p2ps, *dtss, *dtps, *rs;
    sf_file mod, vpvs;

    sf_init (argc,argv);
    mod = sf_output("out");

    if (!sf_getint("nr",&nr)) sf_error("Need nr=");
    if (!sf_getint("n1",&nt))   nt=3501;  sf_putint(mod,"n1",nt);
    if (!sf_getfloat("d1",&dt)) dt=0.001; sf_putfloat(mod,"d1",dt);
    if (!sf_getfloat("o1",&t0)) dt=0.001; sf_putfloat(mod,"o1",t0);
    sf_putint(mod,"n2",3);
    sf_setformat(mod,"native_float");

    sf_getfloats ("fo",3,fo); 

    vpvs = sf_output("vpvs");
    sf_putint(vpvs,"n1",nr);
    sf_putint(vpvs,"n2",2);
    sf_setformat(vpvs,"native_float");

    t = sf_floatalloc (nt);
    pp = sf_floatalloc (nt);
    ps = sf_floatalloc (nt);
    ss = sf_floatalloc (nt);
    ts = sf_floatalloc (nt);

    for (it=0; it < nt; it++) {
	t[it] = t0 + it*dt;
    }

    dtpp = sf_floatalloc (nr);
    tm = sf_floatalloc (nr);
    rs = sf_floatalloc (nr);
    tpp = sf_floatalloc (nr);
    tss = sf_floatalloc (nr);
    tps = sf_floatalloc (nr);
    p2ss = sf_floatalloc (nr);
    p2ps = sf_floatalloc (nr);
    dtss = sf_floatalloc (nr);
    dtps = sf_floatalloc (nr);

    /* ts - reflector positions */
    

  call random_number (ts)
  ts = 0.1+0.9*ts
  call qsort_init (ts)
  call qsort ()

  ! dtpp - layer thickness in PP time
  dtpp(1) = ts(1)
  dtpp(2:nr) = ts(2:nr) - ts(1:nr-1)

  ! tm - time in the middle of the layer
  tm=ts-0.5*dtpp
  ! p2ss - Vp/Vs ratio as a function of time
  p2ss=1./(0.12+(0.5*tm))

  call sep_write (tm,"vpvs")
  call sep_write (p2ss,"vpvs")

  p2ps=0.5*(1+p2ss)

  dtss=p2ss*dtpp ! SS thickness
  dtps=p2ps*dtpp ! PS thickness

  ! rs - reflection coefficient
  call randn_number (nr, rs)
  rs = 0.1/nr + 0.05*rs

  call cum_sum(dtpp,tpp)
  call cum_sum(dtps,tps)
  call cum_sum(dtss,tss)

  call synf1(pp,t,fo(1),tpp,rs)
  call synf1(ps,t,fo(2),tps,rs)
  call synf1(ss,t,fo(3),tss,rs)

  call sep_write (pp)
  call sep_write (ps)
  call sep_write (ss)

  call exit (0)
end program Randrefl




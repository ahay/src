/* 1-D finite-difference wave extrapolation */
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
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
int main(int argc, char* argv[]) 
{
    int nx, nt, nk, nft, ix, it, ik ;
    float dt, dx, x, dk, k, tmp, tmpdt, pi=SF_PI;
    float *sig,  *nxt,  *old, *cur;
    sf_complex  *uk; 
    sf_complex  tmpex;
    kiss_fftr_cfg cfg;
    float  *v, *vx; 
    sf_file inp, out, vel, grad;
    bool opt;    /* optimal padding */
#ifdef _OPENMP
    int nth;
#endif
     

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");
    vel = sf_input("vel");   /* velocity */
    grad = sf_input("grad"); /* velocity gradient */

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(vel)) sf_error("Need float input");
    if (!sf_histint(vel,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(vel,"d1",&dx)) sf_error("No d1= in input");
    /*  if (!sf_histint(inp,"n2",&nt)) sf_error("No n2= in input"); */
    /*  if (!sf_histfloat(inp,"d2",&dt)) sf_error("No d2= in input"); */
    if (!sf_getbool("opt",&opt)) opt=true;
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt input");
    if (!sf_getint("nt",&nt)) sf_error("Need nt input");
    /* if y, determine optimal size for efficiency */

    sf_putint(out,"n1",nx);
    sf_putfloat(out,"d1",dx);
/*    sf_putfloat(out,"o1",x0); */
    sf_putint(out,"n2",nt);
    sf_putfloat(out,"d2",dt);
    sf_putfloat(out,"o2",0.0); 

    nft = opt? 2*kiss_fft_next_fast_size((nx+1)/2): nx;
    if (nft%2) nft++;
    nk = nft/2+1;
    dk = 1./(nft*dx);

    uk = sf_complexalloc(nk);

    sig    =  sf_floatalloc(nx);
    old    =  sf_floatalloc(nx);
    //cur    =  sf_complexalloc(nx);
    cur    =  sf_floatalloc(nx);
    nxt    =  sf_floatalloc(nx);

    v = sf_floatalloc(nx);
    vx = sf_floatalloc(nx);

    sf_floatread(v,nx,vel);
    sf_floatread(vx,nx,grad);

    cfg = kiss_fftr_alloc(nft,0,NULL,NULL); 

    sf_floatread(sig,nx,inp);		
    sf_floatwrite(sig,nx,out);

    for (ix=0; ix < nx; ix++) {
	cur[ix] =  sig[ix];
	old[ix] =  0.0; 
    }

#ifdef _OPENMP
#pragma omp parallel
    {nth = omp_get_num_threads();}
    sf_warning("using %d threads",nth);
#endif

    /* dt=0.0;  DEBUG */
    /* propagation in time */
    for (it=1; it < nt; it++) {


//	kiss_fft_stride(cfg,(kiss_fft_cpx *)cur,(kiss_fft_cpx *)uk,1);/*compute  u(k) */
	kiss_fftr(cfg,cur,(kiss_fft_cpx*)uk);/*compute  u(k) */
#ifdef _OPENMP
#pragma omp parallel for private(ik,ix,x,k,tmp,tmpex,tmpdt) 
#endif

	for (ix=0; ix < nx; ix++) {
	    nxt[ix] = 0.0;
	    // x = x0 + ix * dx;
	    x =  ix * dx;
#ifdef SF_HAS_COMPLEX_H
	    for (ik=0; ik < nk; ik++) {
		// k = (k0 + ik * dk)*2.0*pi;
		k =  ik * dk*2.0*pi;
		tmpdt = v[ix]*fabs(k)*dt;
		tmp = x*k +0.5*v[ix]*(vx[ix]*k)*dt*dt;
		tmpex = sf_cmplx(cosf(tmp),sinf(tmp));
		if (ik == 0 || ik == nk-1) nxt[ix] += creal(uk[ik]*tmpex)*cosf(tmpdt);
		else  nxt[ix] += creal(uk[ik]*tmpex)*cosf(tmpdt)*2.0;
	    }

#else
	    for (ik=0; ik < nk; ik++) {
		//k = (k0 + ik * dk)*2.0*pi;
		k =  ik * dk*2.0*pi;
		tmpdt = v[ix]*fabs(k)*dt;
		tmp = x*k +0.5*v[ix]*(vx[ix]*k)*dt*dt;
		tmpex = sf_cmplx(cosf(tmp),sinf(tmp));
                if (ik == 0 || ik == nk-1) nxt[ix] += sf_crealf(sf_crmul(sf_cmul(uk[ik],tmpex),cosf(tmpdt)));
                else nxt[ix] += sf_crealf(sf_crmul(sf_cmul(uk[ik],tmpex),cosf(tmpdt)*2.0));
	    }
#endif
	    nxt[ix] /= (nk-1);
	    nxt[ix] -= old[ix];
	}
	sf_floatwrite(nxt,nx,out);

	for(ix=0; ix<nx; ix++){
	    old[ix] = cur[ix];
	    cur[ix] = nxt[ix];
	}
    }

    free(v);     
    free(vx);     
    free(nxt);     
    free(cur);     
    free(old);     
    free(uk);     
    free(sig);


    exit(0); 
}           
           

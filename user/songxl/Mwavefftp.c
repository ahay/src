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
    int nx, nt, nk, ix, it, ik, sl;
    float dt, dx, x,  dk, k0, k, tmp, tmpdt, pi=SF_PI;
    float *sig, *wav;
    sf_complex   *new, *cur,*old,*uk; 
    sf_complex  tmpex;
    kiss_fft_cfg cfg;
    float  *v, *vx; 
    sf_file inp, out, vel, grad;
    bool opt;    /* optimal padding */
    int npad;  /* padding */
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
    if (!sf_histint(inp,"n1",&nt)) sf_error("No n2= in input");
    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d2= in input");
    if (!sf_getbool("opt",&opt)) opt=true;
    /* if y, determine optimal size for efficiency */
    if (!sf_getint("pad",&npad)) npad=1;
    if (!sf_getint("sl",&sl)) sl=nx/2;

    sf_putint(out,"n1",nx);
    sf_putfloat(out,"d1",dx);
    sf_putint(out,"n2",nt);
    sf_putfloat(out,"d2",dt);
    sf_putfloat(out,"o2",0.0); 
    nk = opt? kiss_fft_next_fast_size(nx*npad): nx*npad;
    if (nk != nx) sf_warning("padded to %d",nk);
    
    uk = sf_complexalloc(nk);

    new    =  sf_complexalloc(nx);
    cur    =  sf_complexalloc(nx);
    old    =  sf_complexalloc(nx);
    sig    =  sf_floatalloc(nt);
    wav    =  sf_floatalloc(nx);

    v = sf_floatalloc(nx);
    vx = sf_floatalloc(nx);

    sf_floatread(v,nx,vel);
    sf_floatread(vx,nx,grad);
    sf_floatread(sig,nt,inp);		

    dk = 1./(nk*dx);
    k0 = -0.5/dx;
    cfg = kiss_fft_alloc(nk,0,NULL,NULL); 

    for (ix=0; ix < nx; ix++) {
	new[ix] =  sf_cmplx(0.0,0.0);
	cur[ix] =  sf_cmplx(0.0,0.0);
    }

#ifdef _OPENMP
#pragma omp parallel
    {nth = omp_get_num_threads();}
    sf_warning("using %d threads",nth);
#endif
    /* dt=0.0;  DEBUG */
    /* propagation in time */
    for (it=0; it < nt; it++) {

	//     new[sl] += sig[it];
        
        for (ix=0; ix < nx; ix++) {
	    old[ix] =  cur[ix];
	    cur[ix] =  new[ix];
	}

	kiss_fft_stride(cfg,(kiss_fft_cpx *)cur,(kiss_fft_cpx *)uk,1);/*compute  u(k) */

#ifdef _OPENMP
#pragma omp parallel for private(ik,ix,x,k,tmp,tmpdt,tmpex) 
#endif
	for (ix=0; ix < nx; ix++) {
	    new[ix] = sf_cmplx(0.0,0.0);
	    for (ik=0; ik < nk; ik++) {
		x = ix * dx;
		k = (k0 + ik * dk)*2.0*pi;
		tmpdt = v[ix]*fabs(k)*dt;
		tmp = x*k +0.5*v[ix]*(vx[ix]*k)*dt*dt;
		tmpex = sf_cmplx(cos(tmp),sin(tmp))*cos(tmpdt)*2.0;
#ifdef SF_HAS_COMPLEX_H
		tmpex = sf_cmplx(cos(tmp),sin(tmp))*cos(tmpdt)*2.0;
		new[ix] += tmpex*uk[ik];
#else
		tmpex = sf_crmul(sf_cmplx(cos(tmp),sin(tmp)),cos(tmpdt)*2.0);
		new[ix] = sf_cadd(new[ix],sf_cmul(tmpex,uk[ik]));
#endif
	    }
#ifdef SF_HAS_COMPLEX_H
	    new[ix] /= (float)nk;
	    new[ix] -= old[ix];
#else
	    new[ix] = sf_crmul(new[ix],1.0/nk);
	    new[ix] = sf_cadd(new[ix],sf_crmul(old[ix],-1.0));
#endif
	}

#ifdef SF_HAS_COMPLEX_H
	new[sl] += sf_cmplx(sig[it],0.0);
#else
	new[sl] = sf_cadd(new[0],sf_cmplx(sig[it],0.0)); 
#endif

#ifdef SF_HAS_COMPLEX_H
	for(ix=0;ix<nx;ix++) wav[ix]=creal(new[ix]);
#else
	for(ix=0;ix<nx;ix++) wav[ix]=sf_crealf(new[ix]);
#endif
	for(ix=0;ix<nx;ix++) new[ix]=sf_cmplx(wav[ix],0.0);
        sf_floatwrite(wav,nx,out);
    }

    free(v);     
    free(vx);     
    free(new);     
    free(cur);     
    free(old);     
    free(uk);     
    free(sig);
    free(wav);

    sf_close();
    exit(0); 
}           
           

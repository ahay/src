/* 1-D Lowrank Fourier finite-difference wave extrapolation */
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

int main(int argc, char* argv[]) 
{
    int nx, nt, nk, ik, ix, it, nft;
    float dt, dx, dk, k, k0;
    float *old, *nxt, *cur, *sig, *v, *nxttmp, v0, pi=SF_PI, tmpk;
    float *a, *b1, *b2, *b3, *b4, *b5;
    sf_file in, out, Gmatrix, vel;
    sf_complex  *uk, *uktmp; 
    kiss_fftr_cfg cfg, cfgi;
    int im,im2,im3,im4,im5,ip,ip2,ip3,ip4,ip5;
    float factor, kmax;

    sf_init(argc,argv);
    in  = sf_input("in");
    vel = sf_input("vel");   /* velocity */
    Gmatrix = sf_input("G");   /* velocity */
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(Gmatrix)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(vel)) sf_error("Need float input");
    if (!sf_histint(vel,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(vel,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_getint("nt",&nt)) sf_error("No nt in input");
    if (!sf_getfloat("dt",&dt)) sf_error("No dt in input");
    if (!sf_getfloat("factor",&factor)) factor=0.5;

    sf_putint(out,"n1",nx);
    sf_putfloat(out,"d1",dx);
    sf_putint(out,"n2",nt);
    sf_putfloat(out,"d2",dt);
    sf_putfloat(out,"o2",0.0); 

    nk = nx/2+1;
    //nk = nx;
    dk = 1./(nx*dx);
    nft = nx;
    kmax = 1./(dx*2.)*2.0*pi;
    k0 = 0.;
    sig = sf_floatalloc(nx);
    old = sf_floatalloc(nx);
    nxt = sf_floatalloc(nx);
    nxttmp = sf_floatalloc(nx);
    cur = sf_floatalloc(nx);
    v = sf_floatalloc(nx);
    a = sf_floatalloc(nx);
    b1 = sf_floatalloc(nx);
    b2 = sf_floatalloc(nx);
    b3 = sf_floatalloc(nx);
    b4 = sf_floatalloc(nx);
    b5 = sf_floatalloc(nx);
    uk = sf_complexalloc(nk);
    uktmp = sf_complexalloc(nk);
    
 
    sf_floatread(v,nx,vel);
    sf_floatread(sig,nx,in);		
    sf_floatread(a,nx,Gmatrix);
    sf_floatread(b1,nx,Gmatrix);
    sf_floatread(b2,nx,Gmatrix);
    sf_floatread(b3,nx,Gmatrix);
    sf_floatread(b4,nx,Gmatrix);
    sf_floatread(b5,nx,Gmatrix);

    v0=0.0;
    for (ix=0; ix < nx; ix++) v0 += v[ix]*v[ix];
    v0 /= (float)nx;
    /* v0 RMS velocity*/
    v0 = sqrtf(v0);

	/* initial conditions */
    for (ix=0; ix < nx; ix++){
        cur[ix] =  sig[ix];
        old[ix] =  0.0; 
	nxt[ix] = 0.;
    }

    free(v);
    free(sig);

    cfg = kiss_fftr_alloc(nft,0,NULL,NULL); 
    cfgi = kiss_fftr_alloc(nft,1,NULL,NULL); 
    /* propagation in time */
    for (it=0; it < nt; it++) {
        sf_floatwrite(cur,nx,out);
	kiss_fftr(cfg,cur,(kiss_fft_cpx*)uk);/*compute  u(k) */
        for (ik=0; ik < nk; ik++) {
            k = k0 + ik * dk*2.0*pi;
            tmpk = v0*k*dt;
#ifdef SF_HAS_COMPLEX_H
            uktmp[ik] = uk[ik]*2.0*(cosf(tmpk));
#else
	    uktmp[ik] = sf_crmul(uk[ik],2.0*(cosf(tmpk)));
#endif
	    if (fabs(k) > factor*kmax) {uktmp[ik]=sf_cmplx(0.,0.); } 
	}
	kiss_fftri(cfgi,(kiss_fft_cpx*)uktmp,nxttmp);

	for (ix=0; ix < nx; ix++) nxttmp[ix] /= (float)nft; 

	/* Stencil */
	for (ix=0; ix < nx; ix++) {
            im  =(ix-1+nx)%nx; 
            im2 =(ix-2+nx)%nx; 
            im3 =(ix-3+nx)%nx; 
            im4 =(ix-4+nx)%nx; 
            im5 =(ix-5+nx)%nx; 
            ip = (ix+1+nx)%nx;
            ip2 =(ix+2+nx)%nx;
            ip3 =(ix+3+nx)%nx;
            ip4 =(ix+4+nx)%nx;
            ip5 =(ix+5+nx)%nx;
	    nxt[ix] = nxttmp[ix]*a[ix]  
                    + (nxttmp[im]+nxttmp[ip])*0.5*b1[ix]
                    + (nxttmp[im2]+nxttmp[ip2])*0.5*b2[ix]
                    + (nxttmp[im3]+nxttmp[ip3])*0.5*b3[ix]
                    + (nxttmp[im4]+nxttmp[ip4])*0.5*b4[ix]
                    + (nxttmp[im5]+nxttmp[ip5])*0.5*b5[ix];
	}
	
	for (ix=0; ix < nx; ix++) {
	    nxt[ix] +=   - old[ix];
	}

	for (ix=0; ix < nx; ix++) {
	    old[ix] = cur[ix];
	    cur[ix] = nxt[ix];
	}
    }


    exit(0);
}

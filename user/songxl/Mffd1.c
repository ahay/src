/* 1-D Fourier finite-difference wave extrapolation */
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
    float dt, dx, dk, k;
    float *old, *nxt, *cur, *sig, *v, *nxttmp, v0, **aa, tv, tv0, dv, pi=SF_PI, tmpk;
    sf_file in, out, vel;
    bool opt;    /* optimal padding */
    sf_complex  *uk, *uktmp; 
    kiss_fftr_cfg cfg, cfgi;

    sf_init(argc,argv);
    in  = sf_input("in");
    vel = sf_input("vel");   /* velocity */
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(vel)) sf_error("Need float input");
    if (!sf_histint(vel,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(vel,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_getbool("opt",&opt)) opt=true;
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt input");
    if (!sf_getint("nt",&nt)) sf_error("Need nt input");

    sf_putint(out,"n1",nx);
    sf_putfloat(out,"d1",dx);
//    sf_putfloat(out,"o1",x0);
    sf_putint(out,"n2",nt);
    sf_putfloat(out,"d2",dt);
    sf_putfloat(out,"o2",0.0); 

    nft = opt? 2*kiss_fft_next_fast_size((nx+1)/2): nx;
    if (nft%2) nft++;
    nk = nft/2+1;
    dk = 1./(nft*dx);

    sig = sf_floatalloc(nx);
    old = sf_floatalloc(nx);
    nxt = sf_floatalloc(nx);
    nxttmp = sf_floatalloc(nx);
    cur = sf_floatalloc(nx);
    v = sf_floatalloc(nx);
    aa = sf_floatalloc2(2,nx);
    uk = sf_complexalloc(nk);
    uktmp = sf_complexalloc(nk);
    
 
    sf_floatread(v,nx,vel);
    sf_floatread(sig,nx,in);		
    sf_floatwrite(sig,nx,out);

    v0=0.0;
    for (ix=0; ix < nx; ix++) v0 += v[ix]*v[ix];
    v0 /= (float)nx;
    /* v0 RMS velocity*/
    v0 = sqrt(v0);
    tv0 = v0*v0*dt*dt;

    for (ix=0; ix < nx; ix++){
        tv = dt*dt*v[ix]*v[ix];
        dv = (v[ix]*v[ix]-v0*v0)*dt*dt/(dx*dx);
        aa[ix][0] = tv*(1.0 - dv/6.0);
        aa[ix][1] = tv*dv/12.0;
    } 
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
    for (it=1; it < nt; it++) {
	kiss_fftr(cfg,cur,(kiss_fft_cpx*)uk);/*compute  u(k) */
#ifdef SF_HAS_COMPLEX_H
        for (ik=0; ik < nk; ik++) {
            k =  ik * dk*2.0*pi;
            tmpk = v0*fabs(k)*dt;
            uktmp[ik] = uk[ik]*2.0*(cosf(tmpk)-1.0)/tv0;
         }

#else
         for (ik=0; ik < nk; ik++) {
             k =  ik * dk*2.0*pi;
             tmpk = v0*fabs(k)*dt;
             uktmp[ik] = sf_crmul(uk[ik],2.0*(cosf(tmpk)-1.0)/tv0);
         }
#endif
	 kiss_fftri(cfgi,(kiss_fft_cpx*)uktmp,nxttmp);

	for (ix=0; ix < nx; ix++) nxttmp[ix] /= (float)nft; 

	/* Stencil */
	nxt[0] = nxttmp[0]*aa[0][0] + nxttmp[0]*aa[0][0] + nxttmp[1]*aa[0][1];
	for (ix=1; ix < nx-1; ix++) {
	    nxt[ix] = nxttmp[ix]*aa[ix][0] + nxttmp[ix+1]*aa[ix][1] + nxttmp[ix-1]*aa[ix][1];
	}
	nxt[nx-1] = nxttmp[nx-1]*aa[nx-1][1] + nxttmp[nx-1]*aa[nx-1][0] + nxttmp[nx-2]*aa[nx-1][1];
	
	for (ix=0; ix < nx; ix++) {
	    nxt[ix] +=  2*cur[ix] - old[ix];
	}

	sf_floatwrite(nxt,nx,out);
	for (ix=0; ix < nx; ix++) {
	    old[ix] = cur[ix];
	    cur[ix] = nxt[ix];
	}
    }


    exit(0);
}

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
#include <limits.h>
#ifdef _OPENMP
#include <omp.h>
#endif
int main(int argc, char* argv[]) 
{
    int nx, nt, nk, nft, ix, it, ik, iv, nv, dv ;
    float dt, dx, dk, k,  tmpdt, tmpdt2, pi=SF_PI, lmdv, lmdvx, wsum;
    float *sig,  *nxt,  *old, *cur, *nxtc;
    sf_complex  *uk,*uktmp; 
    kiss_fftr_cfg cfg,cfgi;
   // float  *v, *vx, *weight; 
    float  *v, *vc, **weight, *vx, *vxc; 
    sf_file inp, out, vel, grad;
    bool opt;    /* optimal padding */
   // #ifdef _OPENMP
   // int nth;
   // #endif
     

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");
    vel = sf_input("vel");   /* velocity */
    grad = sf_input("grad"); /* velocity gradient */

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(vel)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(grad)) sf_error("Need float input");
    if (!sf_histint(vel,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(vel,"d1",&dx)) sf_error("No d1= in input");
  //  if (!sf_histint(inp,"n2",&nt)) sf_error("No n2= in input");
  //  if (!sf_histfloat(inp,"d2",&dt)) sf_error("No d2= in input");
    if (!sf_getbool("opt",&opt)) opt=true;
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt input");
    if (!sf_getint("nt",&nt)) sf_error("Need nt input");
    if (!sf_getint("nv",&nv)) sf_error("Need nv input");
    if (!sf_getfloat("lmdv",&lmdv)) sf_error("Need lmdv input");
    if (!sf_getfloat("lmdvx",&lmdvx)) sf_error("Need lmdvx input");
    /* if y, determine optimal size for efficiency */

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

    uk = sf_complexalloc(nk);
    uktmp = sf_complexalloc(nk);

    sig    =  sf_floatalloc(nx);
    old    =  sf_floatalloc(nx);
    cur    =  sf_floatalloc(nx);
    nxt    =  sf_floatalloc(nx);
    nxtc   =  sf_floatalloc(nx);
    vc     =  sf_floatalloc(nv);
    vxc     =  sf_floatalloc(nv);
    weight =  sf_floatalloc2(nx,nv);
    
    v = sf_floatalloc(nx);
    vx = sf_floatalloc(nx);

    sf_floatread(v,nx,vel);
    sf_floatread(vx,nx,grad);

/*
    vmax = -FLT_MAX;
    vmin = +FLT_MAX;
    for (ix=0; ix < nx; ix++) {
        if (v[ix] > vmax) vmax = v[ix];
        if (v[ix] < vmin) vmin = v[ix];
        }
    dv = (vmax-vmin)/(nv-1);
*/
    vc[0] = v[0];
    vc[nv-1] = v[nx-1];
    vxc[0] = vx[0];
    vxc[nv-1] = vx[nx-1];
    
    if(nx % (nv-1) == 0)  dv = nx/(nv -1);
       else   dv = nx/(nv -1)+1;
       for (iv=1; iv < nv-1; iv++) { 
           vc[iv] = v[dv*iv]; 
           vxc[iv] = vx[dv*iv];
           }
    for (ix=0; ix < nx; ix++) {
         wsum = 0.0;
         for (iv=0; iv < nv; iv++) {
             weight[iv][ix] = 1.0/((v[ix]-vc[iv])*(v[ix]-vc[iv])*lmdv+(vx[ix]-vxc[iv])*(vx[ix]-vxc[iv])*lmdvx+1.0); // Weight Function
             wsum += weight[iv][ix];
            }
         for (iv=0; iv < nv; iv++) weight[iv][ix] /= wsum;
        }
    cfg = kiss_fftr_alloc(nft,0,NULL,NULL); 
    cfgi = kiss_fftr_alloc(nft,1,NULL,NULL); 

    sf_floatread(sig,nx,inp);		
    sf_floatwrite(sig,nx,out);

    for (ix=0; ix < nx; ix++) {
     cur[ix] =  sig[ix];
     old[ix] =  0.0; 
              }
/*
    #ifdef _OPENMP
    #pragma omp parallel
   {nth = omp_get_num_threads();}
    sf_warning("using %d threads",nth);
    #endif
*/
    /* propagation in time */
    for (it=1; it < nt; it++) {

        for (ix=0; ix < nx; ix++) nxt[ix] = 0.0; 

	kiss_fftr(cfg,cur,(kiss_fft_cpx*)uk);/*compute  u(k) */

/*    #ifdef _OPENMP
    #pragma omp parallel for private(ik,ix,x,k,tmp,tmpex,tmpdt) 
    #endif
*/
        for (iv=0; iv < nv; iv++){
#ifdef SF_HAS_COMPLEX_H
              for (ik=0; ik < nk; ik++) {
                  k =  ik * dk*2.0*pi;
                  tmpdt = vc[iv]*fabs(k)*dt;
                  tmpdt2 = 0.5*vc[iv]*vxc[iv]*k*dt*dt;
                  uktmp[ik] = uk[ik] *2.0*cosf(tmpdt)*sf_cmplx(cosf(tmpdt2),sinf(tmpdt2));
                  }

#else
              for (ik=0; ik < nk; ik++) {
                  k =  ik * dk*2.0*pi;
                  tmpdt = vc[iv]*fabs(k)*dt;
                  uktmp[ik] = sf_cmul(sf_crmul(uk[ik],2.0*cosf(tmpdt)),sf_cmplx(cosf(tmpdt2),sinf(tmpdt2)));
                  }
#endif
	      kiss_fftri(cfgi,(kiss_fft_cpx*)uktmp,nxtc);/*compute  u(k) */
	      for (ix=0; ix < nx; ix++) {  
                   nxtc[ix] /= nft; 
		   nxtc[ix] -= old[ix];
		   nxtc[ix] *= weight[iv][ix];
                   nxt[ix] += nxtc[ix];
                   }
               }  
            for(ix=0; ix<nx; ix++){
	     old[ix] = cur[ix];
	     cur[ix] = nxt[ix];
              }
         sf_floatwrite(nxt,nx,out);
         }

   free(v);     
   free(nxt);     
   free(nxtc);     
   free(cur);     
   free(old);     
   free(uk);     
   free(uktmp);     
   free(sig);

 
   exit(0); 
}           
           

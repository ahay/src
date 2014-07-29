/* 2-D least-squares Kirchhoff pre-stack time migration with different regul..
regularization (or preconditioning) operator: 
reg=0: no regularization; 
reg=1: regularization => first derivative along offset axis;
reg=2: precondition => causual integration along offset axis;
reg=3: precondition => triangle smoother along offset axis, 
reg=4: precondition => local slope constraints along t-x plane and smoothing along offset axis 
*/
/*
  Copyright (C) 2012 China University of Petrolum (East China)
  
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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "tkirmig.h"
#include "causinth.h"
#include "tcaih.h"
#include "triangleh.h"
#include "predicth.h"
#include "predicth2.h"

int main(int argc, char* argv[])
{
    int nt, ncmp, ncdp, nh, nh2, nm, nd, memsize, niter, reg, ix, ih, i3, i2, i1, iter, filt, nw, np;
    float t0, cmp0, cdp0, h0, dt, dcmp, dcdp, dh, apt, rho, aal, norm;
    bool verb, half, amp;
    float ***data, ***modl, **vrms, **mask, *off, *error=NULL;
    float **pp, **qq, *aa;
    char *errfile;
    sf_file in, out, vel, offset, err=NULL;
    sf_file fdip;

    int ompchunk = 1;
    int ompnth = 1;

#ifdef _OPENMP
    int ompath=1;
#endif

    /*------------------------------------------------------------*/
    sf_init(argc,argv);

    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;
    /* OpenMP data chunk size */
#ifdef _OPENMP
    if(! sf_getint("ompnth",  &ompnth))     ompnth=0;
    /* OpenMP available threads */

#pragma omp parallel
    ompath=omp_get_num_threads();
    if(ompnth<1) ompnth=ompath;
    omp_set_num_threads(ompnth);
    sf_warning("using %d threads of a total of %d",ompnth,ompath);
#endif

    in = sf_input("in");
    vel = sf_input("vel");
    out = sf_output("out");

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    if (!sf_getbool("half",&half)) half = true;
    /* if y, the third axis is half-offset instead of full offset */

    if (!sf_getbool("amp",&amp)) amp = true;
    /* if y, use amplitue factor */

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&t0)) sf_error("No o1= in input");

    if (!sf_histint(in,"n2",&ncmp)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&dcmp)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o2",&cmp0)) sf_error("No o2= in input");

    if (!sf_getint("ncdp",&ncdp)) ncdp = ncmp;
    if (!sf_getfloat("dcdp",&dcdp)) dcdp = dcmp;
    if (!sf_getfloat("cdp0",&cdp0)) cdp0 = cmp0;

    sf_putint(out,"n2",ncdp);
    sf_putfloat(out,"d2",dcdp);
    sf_putfloat(out,"o2",cdp0);

    if (!sf_histint(in,"n3",&nh)) sf_error("No n3= in input");

    if (NULL != sf_getstring("offset")) {
        offset = sf_input("offset");
        nh2 = sf_filesize(offset);

        if (nh2 != nh*ncmp) sf_error("Wrong dimensions in offset, it should be %d",nh*ncmp);

        off = sf_floatalloc(nh2);
        sf_floatread (off,nh2,offset);
        sf_fileclose(offset);
        if (!half) {
           for (ih = 0; ih < nh2; ih++) {
               off[ih] *= 0.5;
            }
        }
    } else {
        if (!sf_histfloat(in,"o3",&h0)) sf_error("No o3=");
        if (!sf_histfloat(in,"d3",&dh)) sf_error("No d3=");

        if (!half) dh *= 0.5,h0 *= 0.5;

        off = sf_floatalloc(nh*ncmp);
        for (ix = 0; ix < ncmp; ix++) {
            for (ih = 0; ih < nh; ih++) {
                off[ih*ncmp+ix] = h0 + ih*dh;
            }
        }
        offset = NULL;
    }

    if (!sf_getint("reg",&reg)) reg=0;
    /* regularization type */ 

    if (!sf_getfloat("antialias",&aal)) aal = 1.0;
    /* antialiasing */

    if (!sf_getfloat("apt",&apt)) apt=ncmp;
    /* migration aperture */

    if (!sf_getfloat("rho",&rho)) rho = 1.-1./nt;
    /* Leaky integration constant */

    if (!sf_getint("niter",&niter)) niter=5;
    /* number of iterations */


    nm = nt*ncdp*nh;
    nd = nt*ncmp*nh;

    vrms = sf_floatalloc2(nt,ncdp);
    mask = sf_floatalloc2(ncmp,nh);
    data = sf_floatalloc3(nt,ncmp,nh);
    modl = sf_floatalloc3(nt,ncdp,nh);

    /* read velocity file */
    sf_floatread(vrms[0],nt*ncdp,vel);
    sf_fileclose(vel);

    memsize = nm+nd+nt*ncdp+ncmp*nh;
    if (verb) sf_warning("memory needs: %f G (%f M)",4.*memsize/1024/1024/1024,4.*memsize/1024/1024);

    if (niter > 0) {
        errfile = sf_getstring("err");
        /* output file for error */
        if (NULL != errfile) {
            err = sf_output(errfile);
            sf_putint(err,"n1",niter);
            sf_putfloat(err,"d1",1);
            sf_putfloat(err,"o1",1);
            sf_putstring(err,"label1","Iteration Number");
            sf_putstring(err,"label2","Relative Squared Error");
            sf_putint(err,"n2",1);
            sf_putint(err,"n3",1);
        }
        error = sf_floatalloc(niter);
    }

    sf_floatread(data[0][0],nd,in);

    for (i3=0; i3 < nh; i3++) {
        for (i2=0; i2 < ncmp; i2++) {
            mask[i3][i2]=cblas_sdot(nt,data[i3][i2],1,data[i3][i2],1);
        }
    }

    tkirmig_init(ompnth,ompchunk,nt,dt,t0,ncmp,dcmp,cmp0,ncdp,dcdp,cdp0,nh,dh,h0,apt,aal,rho,vrms,off,mask,amp,verb);

    sf_cdstep_init();

    if (verb) sf_warning("Iteration begin...");

    if (reg == 0)
       sf_solver(tkirmig_lop,sf_cdstep,nm,nd,modl[0][0],data[0][0],
                 niter,"nmem",0,"nfreq",niter,"err",error,"end");

    else if (reg == 1) {
       filt=2;
       aa=sf_floatalloc(filt);
       aa[0]=1.;
       aa[1]=-1.;
       tcaih_init(filt,aa,nt,ncdp,nh);
       sf_solver_reg(tkirmig_lop,sf_cdstep,tcaih_lop,nm+filt*nt*ncdp,nm,nd,
                    modl[0][0],data[0][0],niter,0.01,"nmem",0,"nfreq",niter,
                    "err",error,"end");
    }
    else if (reg == 2) {
       sf_causinth_init(nt,ncdp,nh);
       sf_solver_prec(tkirmig_lop,sf_cdstep,sf_causinth_lop,nm,nm,nd,
                      modl[0][0],data[0][0],niter,0.01,"nmem",0,"nfreq",niter,
                      "err",error,"end");
    }
    else if (reg == 3) {
       sf_triangleh_init(3,nt,ncdp,nh);
       sf_solver_prec(tkirmig_lop,sf_cdstep,sf_triangleh_lop,nm,nm,nd,
                      modl[0][0],data[0][0],niter,0.01,"nmem",0,"nfreq",niter,
                      "err",error,"end");
    }

    else if (reg == 4) {
       sf_warning("pwd constraints along t-x plane and smoothing along offset axis");
       if (!sf_getstring("fdip")) sf_error("Need input dip file!");
       if (!sf_getint("nw",&nw)) nw=3;
       fdip = sf_input("fdip");

       if (!sf_histint(fdip,"n3",&np)) np=1;
       sf_warning("np=%d",np);
       pp = sf_floatalloc2(nt,ncdp);

       if (np > 1) {
          qq = sf_floatalloc2(nt,ncdp);
       } else {
          qq = NULL;
       }

       if (NULL != qq) {
          predicth2_init(nt,ncdp,nh,0.1,nw,pp,qq);
       } else {
          predicth_init(nt,ncdp,nh,0.1,nw,1,false);
          predict_set(pp);
       }

       sf_floatread(pp[0],nt*ncdp,fdip);

       if (NULL != qq) {
          sf_floatread(qq[0],nt*ncdp,fdip);
          sf_solver_prec(tkirmig_lop,sf_cdstep,predicth2_lop,nm,nm,nd,
                      modl[0][0],data[0][0],niter,0.01,"nmem",0,"nfreq",niter,
                      "err",error,"end");
         predict2_close();
       } else {
         sf_solver_prec(tkirmig_lop,sf_cdstep,predicth_lop,nm,nm,nd,
                      modl[0][0],data[0][0],niter,0.01,"nmem",0,"nfreq",niter,
                      "err",error,"end");
         predict_close();
      }
    }

    sf_cdstep_close();
    sf_floatwrite(modl[0][0],nm,out);

    if (NULL != err) {
       for (i3=0; i3 < nh; i3++) {
           for (i2=0; i2 < ncmp; i2++) {
               for (i1=0; i1 < nt; i1++) {
                   norm += data[i3][i2][i1]*data[i3][i2][i1];
               }
           }
        }
        
        for (iter=0; iter < niter; iter++) error[iter] /=norm;
        sf_floatwrite(error,niter,err);
    }

    sf_warning("iter/niter=%d/%d, err=%f",iter,niter,error);

    exit(0);
}

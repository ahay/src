/* 2-D passive seismic RTM and its adjoint */
/*
  Copyright (C) 2016 University of Texas at Austin
  
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

#include "triutil.h"

int main(int argc, char* argv[])
{
    bool verb,adj,abc,inv,prec,sw,ctr;  /* execution flags */
    int nt, nx, nz, depth, nb, n2;      /* dimensions */
    int niter, ngrp, size;              /* # of iters, groups, sw size */
    int rectz, rectx, rectt,repeat;     /* smoothing pars */
    int stack,is,it,ix,iz,tsize;        /* local stacking length */
    float ox, oz, dx, dz, dt, cb;       /* intervals */
    float perc, hard;                   /* hard thresholding and division padding */
    float **dd, **vv, ***ww, ***mwt;    /* arrays */
    float ***ww2=NULL;
    sf_file in, out, vel, weight;       /* I/O files */

    /* initialize Madagascar */
    sf_init(argc,argv);
    
    /* initialize OpenMP support */
#ifdef _OPENMP
    omp_init();
#endif

    if(!sf_getbool("verb", &verb)) verb=false;  /* verbosity flag */
    if(!sf_getbool("adj", &adj)) adj=false;     /* adjoint flag, 0: modeling, 1: migration */
    if(!sf_getbool("abc",&abc)) abc = false;    /* absorbing boundary condition */
    if(!sf_getbool("inv", &inv)) inv=false;     /* inversion flag */
    if(!sf_getbool("prec", &prec)) prec=false;  /* use ctr as precondioner */
    if(!sf_getbool("sw", &sw)) sw=false;        /* inversion flag */
    if(!sf_getbool("ctr", &ctr)) ctr=false;     /* CTR IC flag */
    if(!sf_getint("depth", &depth)) depth=0;    /* acquisition surface */
    if(!sf_getint("niter", &niter)) niter=0;    /* number of iterations */
    if(!sf_getint("ngrp", &ngrp)) ngrp=1;       /* number of groups of receivers */
    if(!sf_getint("size", &size)) size=0;       /* sliding window size */
    if(!sf_getint("rectz", &rectz)) rectz=1;    /* smoothing radius in z */
    if(!sf_getint("rectx", &rectx)) rectx=1;    /* smoothing radius in x */
    if(!sf_getint("rectt", &rectt)) rectt=1;    /* smoothing radius in t */
    if(!sf_getint("repeat", &repeat)) repeat=1; /* smoothing repeatation */
    if(!sf_getint("stack", &stack)) stack=1;    /* local stacking length */
    if(!sf_getfloat("perc", &perc)) perc=SF_EPS;/* stable division padding percentage (of max) */
    if(!sf_getfloat("hard", &hard)) hard=0.0f;  /* hard thresholding */

    if (inv) adj = true; /* inv requires the same output dimension as adj */

    /* setup I/O files */
    in  = sf_input("in");
    out = sf_output("out");
    vel = sf_input("velocity");
    
    /* Dimensions */
    if(!sf_histint  (vel, "n1", &nz)) sf_error("No n1= in velocity");
    if(!sf_histint  (vel, "n2", &nx)) sf_error("No n2= in velocity");
    if(!sf_histfloat(vel, "o1", &oz)) sf_error("No o1= in velocity");
    if(!sf_histfloat(vel, "o2", &ox)) sf_error("No o2= in velocity");
    if(!sf_histfloat(vel, "d1", &dz)) sf_error("No d1= in velocity");
    if(!sf_histfloat(vel, "d2", &dx)) sf_error("No d2= in velocity");

    if (adj){ /* migration */
        if(!sf_histint(in, "n1", &nt)) sf_error("No n1= in data");
        if(!sf_histfloat(in, "d1", &dt)) sf_error("No d1= in data");
        if(!sf_histint(in, "n2", &n2) || n2!=nx) sf_error("Need n2=%d in data", nx);

        sf_putint   (out, "n1", nz);
        sf_putfloat (out, "o1", oz);
        sf_putfloat (out, "d1", dz);
        sf_putstring(out, "label1", "Depth");
        sf_putstring(out, "unit1" , "km");
        sf_putint   (out, "n2", nx);
        sf_putfloat (out, "o2", ox);
        sf_putfloat (out, "d2", dx);
        sf_putstring(out, "label2", "Distance");
        sf_putstring(out, "unit2" , "km");
        sf_putint   (out, "n3", nt/stack);
        sf_putfloat (out, "d3", dt);
        sf_putfloat (out, "o3", 0.0f);
        sf_putstring(out, "label3", "Time");
        sf_putstring(out, "unit3" , "s");
    }else{ /* modeling */
        if(!sf_histint(in, "n3", &nt)) sf_error("No n3= in src");
        if(!sf_histfloat(in, "d3", &dt)) sf_error("No d3= in src");
        
        sf_putint   (out, "n1", nt);
        sf_putfloat (out, "d1", dt);
        sf_putfloat (out, "o1", 0.0);
        sf_putstring(out, "label1", "Time");
        sf_putstring(out, "unit1" , "s");
        sf_putint   (out, "n2", nx);
        sf_putfloat (out, "o2", ox);
        sf_putfloat (out, "d2", dx);
        sf_putstring(out, "label2", "Distance");
        sf_putstring(out, "unit2" , "km");
        sf_putint   (out, "n3", 1);
    }

    if (inv && prec) {
        if (NULL!=sf_getstring("weight")) {
//            weight = sf_input("weight");
            weight = sf_output("weight");

            sf_putint   (weight, "n1", nz);
            sf_putfloat (weight, "o1", oz);
            sf_putfloat (weight, "d1", dz);
            sf_putstring(weight, "label1", "Depth");
            sf_putstring(weight, "unit1" , "km");
            sf_putint   (weight, "n2", nx);
            sf_putfloat (weight, "o2", ox);
            sf_putfloat (weight, "d2", dx);
            sf_putstring(weight, "label2", "Distance");
            sf_putstring(weight, "unit2" , "km");
            sf_putint   (weight, "n3", nt);
            sf_putfloat (weight, "d3", dt);
            sf_putfloat (weight, "o3", 0.0f);
            sf_putstring(weight, "label3", "Time");
            sf_putstring(weight, "unit3" , "s");
        }
    } else weight=NULL;
    
    /* padding and abc */
    if(!sf_getint("nb", &nb) || nb<NOP) nb = NOP;
    if(!sf_getfloat("cb", &cb)) cb = 0.0f;
   
    /* allocate arrays */
    vv = sf_floatalloc2(nz, nx);
    dd = sf_floatalloc2(nt, nx);
    ww = sf_floatalloc3(nz, nx, nt);
    if (inv && prec) mwt = sf_floatalloc3(nz, nx, nt);
    if (stack > 1) {
        ww2= sf_floatalloc3(nz, nx, (int)(nt/stack));
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(is,ix,iz)
#endif
        for (is=0; is<nt/stack; is++)
            for (ix=0; ix<nx; ix++)
                for (iz=0; iz<nz; iz++)
                    ww2[is][ix][iz] = 0.;
    }
    /* read velocity */
    sf_floatread(vv[0], nz*nx, vel);
    if (adj) sf_floatread(dd[0], nt*nx, in);
    else sf_floatread(ww[0][0], nz*nx*nt, in);

    /* initialize time-reversal imaging */
    timerev_init(verb, abc, nt, nx, nz, nb, depth, dt, dx, dz, cb, vv);

    /* calculate model weighting using correlative imaging condition */
    if (inv && prec) { 
        //sf_floatread(mwt[0][0], nz*nx*nt, weight);
        if (ctr) {
            ctimerev(ngrp,mwt,dd);
            absval(nz*nx*nt,mwt[0][0]);
        } else {
            timerev_lop(adj, false, nz*nx*nt, nt*nx, mwt[0][0], dd[0]);
            autopow(nz*nx*nt,(float)ngrp,mwt[0][0]);
        }
        /* smoothing */
        smooth(nz, nx, nt, rectz, rectx, rectt, repeat, mwt[0][0]);
        /* local normalizaiton */
        swnorm(verb, sw, nz, nx, nt, size, perc, mwt[0][0]);
        /* hard thresholding */
        if (hard>0) threshold(false, nz*nx*nt, hard, mwt[0][0]);
    }

    /* apply time-reversal imaging linear operator */
    if (inv) {
        if (prec) sf_solver(timerev_lop,sf_cgstep,nz*nx*nt,nt*nx,ww[0][0],dd[0],niter,"mwt",mwt[0][0],"verb",verb,"end");
        else sf_solver(timerev_lop,sf_cgstep,nz*nx*nt,nt*nx,ww[0][0],dd[0],niter,"verb",verb,"end");
    } else {
        if (adj && ctr) ctimerev(ngrp,ww,dd);
        else timerev_lop(adj, false, nz*nx*nt, nt*nx, ww[0][0], dd[0]);
    }

    if (stack > 1) {
        tsize = stack;
        sf_warning("nt/stack=%d",nt/stack);
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(is,it,ix,iz)
#endif
        for (is=0; is<nt/stack; is++) {
            if (is == nt/stack-1) tsize=stack+(nt-(nt/stack)*stack);
            //sf_warning("is=%d,tsize=%d",is,tsize);
            for (it=0; it<tsize; it++)
                for (ix=0; ix<nx; ix++)
                    for (iz=0; iz<nz; iz++)
                        ww2[is][ix][iz] += ww[is*stack+it][ix][iz];
        }
    }

    if (adj) {
        if (stack > 1) sf_floatwrite(ww2[0][0], nz*nx*nt/stack, out);
        else sf_floatwrite(ww[0][0], nz*nx*nt, out);
    } else sf_floatwrite(dd[0], nt*nx, out);

    if (NULL!=weight) sf_floatwrite(mwt[0][0], nz*nx*nt, weight);

    /* close */
    timerev_close();
    free(*dd); free(dd); 
    free(**ww); free(*ww); free(ww);
    if (NULL!=ww2) { free(**ww2); free(*ww2); free(ww2); }
    if (inv && prec) { free(**mwt); free(*mwt); free(mwt); }

    exit (0);
}


/* 2-D inversion to common scattering-point gathers
The axes in the data space are {time,offset,cmp}
The axes in the image space are {time,equiv_offset,csp}
*/
/*
  Copyright (C) 2013 China University of Petrolum (East China)
  
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

#include "csp2d.h"
#include "radon2d.h"

int main(int argc, char* argv[])
{
    int nhe, nxs, nh, nxm, nt, nv;
    float dhe, he0, dxs, xs0, dh, h0, dxm, xm0, dt, t0, v, apt, dv, ov, eps;
    int nm, nd, niter;
    bool half, verb, linear, weight, reg, sparse;

    float *cmp, *csp, *error;

    char *errfile;

    sf_file in, out, err;

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
    out = sf_output("out");

    if (!sf_getint("niter",&niter)) niter=true;
    /* iteration number */

    if (!sf_getbool("weight",&weight)) weight=false;
    /* weighting flag */

    if (!sf_getbool("linear",&linear)) linear=true;
    /* yes: linear interpolation, no: nearest-neighbor interpolation */

    if (!sf_getfloat("v",&v)) v=2000.;
    /* velocity */

    if (!sf_getbool("reg",&reg)) reg=false;
    /* regularization flag */

    if (!sf_getbool("sparse",&sparse)) sparse=false;
    /* sparse flag */

    if (!sf_getbool("half",&half)) half=true;
    /* half offset flag */

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

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

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&t0)) sf_error("No o1= in input");

    if (!sf_histint(in,"n3",&nxm)) sf_error("No n3= in input");
    if (!sf_histfloat(in,"d3",&dxm)) sf_error("No d3= in input");
    if (!sf_histfloat(in,"o3",&xm0)) sf_error("No o3= in input");

    if (!sf_histint(in,"n2",&nh)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&dh)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o2",&h0)) sf_error("No o2= in input");

    if (!sf_getint("nhe",&nhe)) nhe = nh;
    if (!sf_getfloat("dhe",&dhe)) dhe = dh;
    if (!sf_getfloat("he0",&he0)) he0 = h0;

    if (!sf_getint("nxs",&nxs)) nxs = nxm;
    if (!sf_getfloat("dxs",&dxs)) dxs = dxm;
    if (!sf_getfloat("xs0",&xs0)) xs0 = xm0;

    sf_putint(out,"n1",nt);
    sf_putfloat(out,"d1",dt);
    sf_putfloat(out,"o1",t0);

    sf_putint(out,"n2",nhe);
    sf_putfloat(out,"d2",dhe);
    sf_putfloat(out,"o2",he0);

    sf_putint(out,"n3",nxs);
    sf_putfloat(out,"d3",dxs);
    sf_putfloat(out,"o3",xs0);

    if (!sf_getfloat("apt",&apt)) apt=SF_MAX(fabsf(he0),fabsf(he0+(nhe-1)*dhe));
    /* aperture */

    nm = nhe*nxs*nt;
    nd = nh*nxm*nt;

    /* Allocate 3-D array for the convenience of inversion with regularization */
    cmp = sf_floatalloc(nd);
    csp = sf_floatalloc(nm);

    if (verb) sf_warning("Memory needs: %f G (%f M)", 4.*(nm+nd)/1024./1024./1024., 
              4.*(nm+nd)/1024./1024);
    if (verb) sf_warning("The aperture is %f", apt);

    sf_floatread(cmp,nd,in);

    if (verb) sf_warning("Inversion to CSP...");

    sf_cdstep_init();

    csp2d_init(ompnth,ompchunk,nhe,dhe,he0,nxs,dxs,xs0,nh,dh,h0,nxm,dxm,xm0,nt,dt,t0,apt,v,
               weight, linear,half,verb);


    if (reg) {
       if (!sf_getint("nv",&nv)) sf_error("No nv in parameters");
       /* scanning velocity number */
       if (!sf_getfloat("dv",&dv)) sf_error("No dv in parameters");
       /* scanning velocity increment */
       if (!sf_getfloat("ov",&ov)) sf_error("No ov in parameters");
       /* scanning velocity origin */

       if (!sf_getfloat("eps",&eps)) eps=0.01;
       /* penalty parameter */

       radon2d_init(ompnth,ompchunk,nxs,nt,t0,dt,nhe,he0,dhe,nv,ov,dv,half,verb);

       if (sparse) {
          int miter;
          float *w;
          if (!sf_getint("miter",&miter)) sf_error("No miter in parameters");
          /* number of nonlinear iteration */
          w = sf_floatalloc(nm);
          for (int im=0; im < nm; im++) w[im] = 1.f;
          for (int iter=0; iter < miter; iter++) {
              //sf_solver_reg(csp2d_lop,sf_cdstep,radon2d_lop,nxs*nv*nt,nm,nd,
              //      csp,cmp,niter,0.01,"nmem",0,"nfreq",niter,"err",error,"mwt", w,"end");
              for (int im=0; im < nm; im++) w[im] = fabsf(csp[im]);
          }

       } else {
         sf_solver_reg(csp2d_lop,sf_cdstep,radon2d_lop,nxs*nv*nt,nm,nd,
                    csp,cmp,niter,eps,"nmem",0,"nfreq",niter,"err",error,"end");
       }
    } else {
      sf_solver(csp2d_lop,sf_cdstep,nm,nd,csp,cmp,
                niter,"nmem",0,"nfreq",niter,"err",error,"end");
    }

    sf_cdstep_close();

    sf_floatwrite(csp,nm,out);

    if (NULL != err) {
       for (int iter=1; iter < niter; iter++) error[iter] /= error[0];
       error[0] = 1.;
       sf_floatwrite(error,niter,err);
    }

    exit(0);

}

/*      $Id: Micsp2d.c 777 2013-11-26 18:46:07Z Yujin Liu $       */

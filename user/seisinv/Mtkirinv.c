/* 2-D least-squares Kirchhoff pre-stack time migration. */ 

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

#include "tkirmig.h"
#include "causinth.h"
#include "tcaih.h"
#include "triangleh.h"
#include "predicth.h"
#include "predicth2.h"

int main(int argc, char* argv[])
{
    int n1, n2, n3, n123, niter, i1, i2, i3, iter;
    float o1, o2, o3, d1, d2, d3, dip, norm;
    bool verb, half;
    int  reg;
    float ***data, ***modl, **vrms, **mask, *error=NULL;
    char *errfile;
    sf_file in, out, vel, err=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    vel = sf_input("vel");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&o1)) sf_error("No o1= in input");

    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o2",&o2)) sf_error("No o2= in input");

    if (!sf_histint(in,"n3",&n3)) sf_error("No n3= in input");
    if (!sf_histfloat(in,"d3",&d3)) sf_error("No d3= in input");
    if (!sf_histfloat(in,"o3",&o3)) sf_error("No o3= in input");

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */
    if (!sf_getint("reg",&reg)) reg=0;
    /* regularization type: 
       reg=0,no regularization; 
       reg=1,regularization,smooth along offset axis;
       reg=2,precondition,smooth along offset axis; */
    if (!sf_getfloat("dip",&dip)) dip=45.;
    /* the dip range of migration aperture */
    if (!sf_getbool("half",&half)) half=false;
    /* half offset flag */
    if (!sf_getint("niter",&niter)) niter=10;
    /* number of iterations */


    vrms = sf_floatalloc2(n1,n2);
    mask = sf_floatalloc2(n2,n3);
    data = sf_floatalloc3(n1,n2,n3);
    modl = sf_floatalloc3(n1,n2,n3);
    n123 = n1*n2*n3;
    /* read velocity file */
    sf_floatread(vrms[0],n1*n2,vel);
    sf_fileclose(vel);

//    sf_triangle1_init(10,n1);

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
        }
        error = sf_floatalloc(niter);
    }

    sf_floatread(data[0][0],n123,in);

    for (i3=0; i3 < n3; i3++) {
        for (i2=0; i2 < n2; i2++) {
            mask[i3][i2]=cblas_sdot(n1,data[i3][i2],1,data[i3][i2],1);
        }
    }

    tkirmig_init(n1,d1,o1,n2,d2,o2,n3,d3,o3,dip,vrms,mask,half,verb);

    sf_cdstep_init();

    if (reg == 0)
       sf_solver(tkirmig_lop,sf_cdstep,n123,n123,modl[0][0],data[0][0],
                 niter,"nmem",0,"nfreq",niter,"err",error,"end");

    else if (reg == 1) {
       int filt; /* the filter length can be changed */
       float *aa;
       filt=2;
       aa=sf_floatalloc(filt);
       aa[0]=1.;
       aa[1]=-1.;
       tcaih_init(filt,aa,n1,n2,n3);
       sf_solver_reg(tkirmig_lop,sf_cdstep,tcaih_lop,n123+filt*n1*n2,n123,n123,
                    modl[0][0],data[0][0],niter,0.01,"nmem",0,"nfreq",niter,
                    "err",error,"end");
    }
    else if (reg == 2) {
       sf_causinth_init(n1,n2,n3);
       sf_solver_prec(tkirmig_lop,sf_cdstep,sf_causinth_lop,n123,n123,n123,
                      modl[0][0],data[0][0],niter,0.01,"nmem",0,"nfreq",niter,
                      "err",error,"end");
    }
    else if (reg == 3) {
       sf_triangleh_init(3,n1,n2,n3);
       sf_solver_prec(tkirmig_lop,sf_cdstep,sf_triangleh_lop,n123,n123,n123,
                      modl[0][0],data[0][0],niter,0.01,"nmem",0,"nfreq",niter,
                      "err",error,"end");
    }

    else if (reg == 4) {
       sf_file fdip;
       int np,nw;
       float **pp, **qq;
       sf_warning("pwd constraints along t-x plane and smoothing along offset axis");
       if (!sf_getstring("fdip")) sf_error("Need input dip file!");
       if (!sf_getint("nw",&nw)) nw=3;
       fdip = sf_input("fdip");

       if (!sf_histint(fdip,"n3",&np)) np=1;
       sf_warning("np=%d",np);
       pp = sf_floatalloc2(n1,n2);

       if (np > 1) {
          qq = sf_floatalloc2(n1,n2);
       } else {
          qq = NULL;
       }

       if (NULL != qq) {
          predicth2_init(n1,n2,n3,0.1,nw,pp,qq);
       } else {
          predicth_init(n1,n2,n3,0.1,nw,1,false);
          predict_set(pp);
       }

       sf_floatread(pp[0],n1*n2,fdip);

       if (NULL != qq) {
          sf_floatread(qq[0],n1*n2,fdip);
          sf_solver_prec(tkirmig_lop,sf_cdstep,predicth2_lop,n123,n123,n123,
                      modl[0][0],data[0][0],niter,0.01,"nmem",0,"nfreq",niter,
                      "err",error,"end");
         predict2_close();
       } else {
         sf_solver_prec(tkirmig_lop,sf_cdstep,predicth_lop,n123,n123,n123,
                      modl[0][0],data[0][0],niter,0.01,"nmem",0,"nfreq",niter,
                      "err",error,"end");
         predict_close();
      }
    }

    sf_cdstep_close();
    sf_floatwrite(modl[0][0],n123,out);

    if (NULL != err) {
       for (i3=0; i3 < n3; i3++) {
           for (i2=0; i2 < n2; i2++) {
               for (i1=0; i1 < n1; i1++) {
                   norm += data[i3][i2][i1];
               }
           }
        }
        
        for (iter=0; iter < niter; iter++) error[iter] /=norm;
        sf_floatwrite(error,niter,err);
    }

    exit(0);
}

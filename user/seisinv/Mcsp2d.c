/* 2-D common scattering-point gathers mapping and its adjoint
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

int main(int argc, char* argv[])
{
    int nhe, nxs, nh, nxm, nt;
    float dhe, he0, dxs, xs0, dh, h0, dxm, xm0, dt, t0, v, apt;
    int nm, nd;
    bool adj, half, verb, linear, weight;

    float *data, *modl;

    sf_file in, out;

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

    if (!sf_getbool("adj",&adj)) adj=true;
    /* yes: CSP mapping, no: CMP building */

    if (!sf_getbool("weight",&weight)) weight=false;
    /* weighting flag */

    if (!sf_getbool("linear",&linear)) linear=false;
    /* yes: linear interpolation, no: nearest-neighbor interpolation */

    if (!sf_getfloat("v",&v)) v=2000.;
    /* velocity */

    if (!sf_getbool("half",&half)) half=false;
    /* half offset flag */

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&t0)) sf_error("No o1= in input");

    sf_putint(out,"n1",nt);
    sf_putfloat(out,"d1",dt);
    sf_putfloat(out,"o1",t0);

    if (adj) {
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

       sf_putint(out,"n2",nhe);
       sf_putfloat(out,"d2",dhe);
       sf_putfloat(out,"o2",he0);

       sf_putint(out,"n3",nxs);
       sf_putfloat(out,"d3",dxs);
       sf_putfloat(out,"o3",xs0);


    } else {
       if (!sf_histint(in,"n3",&nxs)) sf_error("No n3= in input");
       if (!sf_histfloat(in,"d3",&dxs)) sf_error("No d3= in input");
       if (!sf_histfloat(in,"o3",&xs0)) sf_error("No o3= in input");

       if (!sf_histint(in,"n2",&nhe)) sf_error("No n2= in input");
       if (!sf_histfloat(in,"d2",&dhe)) sf_error("No d2= in input");
       if (!sf_histfloat(in,"o2",&he0)) sf_error("No o2= in input");

       if (!sf_getint("nh",&nh)) nh = nhe;
       if (!sf_getfloat("dh",&dh)) dh = dhe;
       if (!sf_getfloat("h0",&h0)) h0 = he0;

       if (!sf_getint("nxm",&nxm)) nxm = nxs;
       if (!sf_getfloat("dxm",&dxm)) dxm = dxs;
       if (!sf_getfloat("xm0",&xm0)) xm0 = xs0;

       sf_putint(out,"n2",nh);
       sf_putfloat(out,"d2",dh);
       sf_putfloat(out,"o2",h0);

       sf_putint(out,"n3",nxm);
       sf_putfloat(out,"d3",dxm);
       sf_putfloat(out,"o3",xm0);

    }

    if (!sf_getfloat("apt",&apt)) apt=SF_MAX(fabsf(he0),fabsf(he0+(nhe-1)*dhe));
    /* aperture */

    nm = nhe*nxs*nt;
    nd = nh*nxm*nt;

    /* Allocate 3-D array for the convenience of inversion with regularization */
    data = sf_floatalloc(nd);
    modl = sf_floatalloc(nm);

    if (verb) sf_warning("Memory needs: %f G (%f M)", 4.*(nm+nd)/1024./1024./1024., 
              4.*(nm+nd)/1024./1024);
    if (verb) sf_warning("The aperture is %f", apt);

    csp2d_init(ompnth,ompchunk,nhe,dhe,he0,nxs,dxs,xs0,nh,dh,h0,nxm,dxm,xm0,nt,dt,t0,apt,v,
               weight, linear,half,verb);

    if (adj) {
       if (verb) sf_warning("CSP mapping...");
       sf_floatread(data,nd,in);
       csp2d_lop(true,false,nm,nd,modl,data);
       sf_floatwrite(modl,nm,out);
    } else {
       if (verb) sf_warning("CMP modeling...");
       sf_floatread(modl,nm,in);
       csp2d_lop(false,false,nm,nd,modl,data);
       sf_floatwrite(data,nd,out);
    }

    exit(0);
}

/*      $Id: Mcsp2d.c 744 2013-07-11 18:46:07Z Yujin Liu $       */

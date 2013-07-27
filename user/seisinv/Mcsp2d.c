/* 2-D common scattering-point gathers mapping and its adjoint */

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

#include "csp2d.h"

int main(int argc, char* argv[])
{
    int nhe, nxs, nh, nxm, nt;
    float dhe, he0, dxs, xs0, dh, h0, dxm, xm0, dt, t0, v;
    int nm, nd;
    bool adj, verb, linear, weight;

    float *data, *modl; /* data set [offset,cmp,t] */

    sf_file in, out;

    sf_init(argc,argv);

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

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    if (!sf_histint(in,"n3",&nt)) sf_error("No n3= in input");
    if (!sf_histfloat(in,"d3",&dt)) sf_error("No d3= in input");
    if (!sf_histfloat(in,"o3",&t0)) sf_error("No o3= in input");

    sf_putint(out,"n3",nt);
    sf_putfloat(out,"d3",dt);
    sf_putfloat(out,"o3",t0);

    if (adj) {
       if (!sf_histint(in,"n1",&nh)) sf_error("No n1= in input");
       if (!sf_histfloat(in,"d1",&dh)) sf_error("No d1= in input");
       if (!sf_histfloat(in,"o1",&h0)) sf_error("No o1= in input");

       if (!sf_histint(in,"n2",&nxm)) sf_error("No n2= in input");
       if (!sf_histfloat(in,"d2",&dxm)) sf_error("No d2= in input");
       if (!sf_histfloat(in,"o2",&xm0)) sf_error("No o2= in input");

       if (!sf_getint("nhe",&nhe)) sf_error("No nhe in parameters");
       if (!sf_getfloat("dhe",&dhe)) sf_error("No dhe in parameters");
       if (!sf_getfloat("he0",&he0)) sf_error("No he0 in parameters");

       if (!sf_getint("nxs",&nxs)) sf_error("No nxs in parameters");
       if (!sf_getfloat("dxs",&dxs)) sf_error("No dxs in parameters");
       if (!sf_getfloat("xs0",&xs0)) sf_error("No xs0 in parameters");

       sf_putint(out,"n1",nhe);
       sf_putfloat(out,"d1",dhe);
       sf_putfloat(out,"o1",he0);

       sf_putint(out,"n2",nxs);
       sf_putfloat(out,"d2",dxs);
       sf_putfloat(out,"o2",xs0);


    } else {
       if (!sf_histint(in,"n1",&nhe)) sf_error("No n1= in input");
       if (!sf_histfloat(in,"d1",&dhe)) sf_error("No d1= in input");
       if (!sf_histfloat(in,"o1",&he0)) sf_error("No o1= in input");

       if (!sf_histint(in,"n2",&nxs)) sf_error("No n2= in input");
       if (!sf_histfloat(in,"d2",&dxs)) sf_error("No d2= in input");
       if (!sf_histfloat(in,"o2",&xs0)) sf_error("No o2= in input");

       if (!sf_getint("nh",&nh)) sf_error("No nh in parameters");
       if (!sf_getfloat("dh",&dh)) sf_error("No dh in parameters");
       if (!sf_getfloat("h0",&h0)) sf_error("No h0 in parameters");

       if (!sf_getint("nxm",&nxm)) sf_error("No nxm in parameters");
       if (!sf_getfloat("dxm",&dxm)) sf_error("No dxm in parameters");
       if (!sf_getfloat("xm0",&xm0)) sf_error("No xm0 in parameters");

       sf_putint(out,"n1",nh);
       sf_putfloat(out,"d1",dh);
       sf_putfloat(out,"o1",h0);

       sf_putint(out,"n2",nxm);
       sf_putfloat(out,"d2",dxm);
       sf_putfloat(out,"o2",xm0);

    }

    nm = nhe*nxs*nt;
    nd = nh*nxm*nt;

    /* Allocate 3-D array for the convenience of inversion with regularization */
    data = sf_floatalloc(nd);
    modl = sf_floatalloc(nm);

    if (verb) sf_warning("memory needs: %f G (%f M)", 4.*(nm+nd)/1024./1024./1024., 
              4.*(nm+nd)/1024./1024);

    csp2d_init(nhe,dhe,he0,nxs,dxs,xs0,nh,dh,h0,nxm,dxm,xm0,nt,dt,t0,v,
               weight, linear,verb);

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

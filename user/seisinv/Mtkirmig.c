/* 2-D Kirchhoff pre-stack time migration/demigration. 
The axes in the data are {time,cmp,offset}
The axes in the offset are {1,cmp,offset}
The axes in the image are {time,cdp,offset}
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

int main(int argc, char* argv[])
{
    int nt, ncmp, ncdp, nh, nh2, nm, nd, memsize, ix, ih, i3, i2;
    float t0, cmp0, cdp0, h0, dt, dcmp, dcdp, dh, apt, rho, aal;
    bool adj, verb, half, amp;
    float *data, *modl, *off, **vrms, **mask;
    sf_file in, out, vel, offset;

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

    if (!sf_getbool("adj",&adj)) adj=true;
    /* yes: migration, no: modeling */

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    if (!sf_getbool("half",&half)) half = true;
    /* if y, the third axis is half-offset instead of full offset */

    if (!sf_getbool("amp",&amp)) amp = true;
    /* if y, use amplitue factor */

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&t0)) sf_error("No o1= in input");

    if (adj) {
    if (!sf_histint(in,"n2",&ncmp)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&dcmp)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o2",&cmp0)) sf_error("No o2= in input");

    if (!sf_getint("ncdp",&ncdp)) ncdp = ncmp;
    if (!sf_getfloat("dcdp",&dcdp)) dcdp = dcmp;
    if (!sf_getfloat("cdp0",&cdp0)) cdp0 = cmp0;

    sf_putint(out,"n2",ncdp);
    sf_putfloat(out,"d2",dcdp);
    sf_putfloat(out,"o2",cdp0);

    } else {

    if (!sf_histint(in,"n2",&ncdp)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&dcdp)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o2",&cdp0)) sf_error("No o2= in input");
    
    if (!sf_getint("ncmp",&ncmp)) ncmp = ncdp;
    if (!sf_getfloat("dcmp",&dcmp)) dcmp = dcdp;
    if (!sf_getfloat("cmp0",&cmp0)) cmp0 = cdp0;

    sf_putint(out,"n2",ncmp);
    sf_putfloat(out,"d2",dcmp);
    sf_putfloat(out,"o2",cmp0);

    }

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

        if (!half) dh *= 0.5, h0 *= 0.5;

        off = sf_floatalloc(nh*ncmp);
        for (ix = 0; ix < ncmp; ix++) {
            for (ih = 0; ih < nh; ih++) {
                off[ih*ncmp+ix] = h0 + ih*dh;
            }
        }
        offset = NULL;
    }

    if (!sf_getfloat("antialias",&aal)) aal = 1.0;
    /* antialiasing */

    if (!sf_getfloat("apt",&apt)) apt=ncmp;
    /* migration aperture */

    if (!sf_getfloat("rho",&rho)) rho = 1.-1./nt;
    /* Leaky integration constant */

    nm = nt*ncdp*nh;
    nd = nt*ncmp*nh;

    vrms = sf_floatalloc2(nt,ncdp);
    mask = sf_floatalloc2(ncmp,nh);
    data = sf_floatalloc(nd);
    modl = sf_floatalloc(nm);

    /* read velocity file */
    sf_floatread(vrms[0],nt*ncdp,vel);

    for (i3=0; i3 < nh; i3++) {
        for (i2=0; i2 < ncmp; i2++) {
            mask[i3][i2]=1.;
        }
    }

    tkirmig_init(ompnth,ompchunk,nt,dt,t0,ncmp,dcmp,cmp0,ncdp,dcdp,cdp0,nh,dh,h0,apt,aal,rho,vrms,off,mask,amp,verb);

    memsize = nd+nm+nt*ncdp+ncmp*nh;
    if (verb) sf_warning("memory needs %f G (%f M)",4.*memsize/1024/1024/1024,4.*memsize/1024/1024);

    if (adj) {
       sf_floatread(data,nd,in);
       tkirmig_lop(true,false,nm,nd,modl,data);
       sf_floatwrite(modl,nm,out);
       
    } else {
      sf_floatread(modl,nm,in);
      tkirmig_lop(false,false,nm,nd,modl,data);
      sf_floatwrite(data,nd,out);
    }

    sf_close();

    exit(0);
}

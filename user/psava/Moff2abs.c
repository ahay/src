/* Transform vector-offset to absolute-offset 
             h = sqrt(hx^2+hy^2+hz^2)
   pcs 2005 
*/

/*
  Copyright (C) 2004 University of Texas at Austin
  
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
    bool  verb;
    sf_axis  az,ah,ahx,ahy,ahz;
    int  iz,   ihx,ihy,ihz;
    int  nz,   nhx,nhy,nhz;
    float       hx, hy, hz;
    float oh,dh,ohx,dhx,ohy,dhy,ohz,dhz;
    sf_bands spl;

    sf_file Fd; /*  data =   vector offset (hx,hy,hz)-z */
    sf_file Fm; /* model = absolute offset     h     -z */

    int nw;    /* spline order */
    int nd,id; /*  data size (nd=nhx*nhy*nhz) */
    int nh;    /* model size (nm=nh) */

    float *dat=NULL;
    float *mod=NULL;
    float *map=NULL;

    float *mwt=NULL;
    float *dwt=NULL;

    int i;
/*    int im;*/
/*------------------------------------------------------------*/

    sf_init(argc,argv);

    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getint(   "nw",&nw))     nw=4;     /* spline order */

    Fd = sf_input ("in");
    ahx = sf_iaxa(Fd,1); sf_setlabel(ahx,"hx"); if(verb) sf_raxa(ahx);
    ahy = sf_iaxa(Fd,2); sf_setlabel(ahy,"hy"); if(verb) sf_raxa(ahy);
    ahz = sf_iaxa(Fd,3); sf_setlabel(ahz,"hz"); if(verb) sf_raxa(ahz);
    az  = sf_iaxa(Fd,4); sf_setlabel(az,"z");   if(verb) sf_raxa(az);

    nhx = sf_n(ahx); ohx = sf_o(ahx); dhx = sf_d(ahx);
    nhy = sf_n(ahy); ohy = sf_o(ahy); dhy = sf_d(ahy);
    nhz = sf_n(ahz); ohz = sf_o(ahz); dhz = sf_d(ahz);
    nz = sf_n(az);

    if(!sf_getint  ("nh",&nh)) nh=nhx + ohx/dhx;
    if(!sf_getfloat("oh",&oh)) oh=0;
    if(!sf_getfloat("dh",&dh)) dh=dhx;
    ah = sf_maxa(nh,oh,dh); sf_setlabel(ah,"h"); if(verb) sf_raxa(ah);

    Fm = sf_output("out");
    sf_oaxa(Fm,ah,1);
    sf_oaxa(Fm,az,2);
    sf_putint(Fm,"n3",1);
    sf_putint(Fm,"n4",1);

/*------------------------------------------------------------*/
    nd = nhx*nhy*nhz;  /*  data size */

    map = sf_floatalloc(nd); /* mapping */

    mod = sf_floatalloc(nh); /* model vector */
    dat = sf_floatalloc(nd); /*  data vector */

    mwt = sf_floatalloc(nh); /* model weight */
    dwt = sf_floatalloc(nd); /*  data weight */

    spl = sf_spline_init(nw,nd);

    for(ihz=0;ihz<nhz;ihz++) {
	hz = ohz + ihz * dhz;         hz*=hz;
	for(ihy=0;ihy<nhy;ihy++) {
	    hy = ohy + ihy * dhy;     hy*=hy;
	    for(ihx=0;ihx<nhx;ihx++) {
		hx = ohx + ihx * dhx; hx*=hx;
		
		i = ihz * nhx*nhy + 
		    ihy * nhx     +
		    ihx;
		
		map[i] = sqrtf(hx+hy+hz);
	    }
	}
    }

    sf_int1_init( map, 
		  oh, dh, nh, 
		  sf_spline_int, 
		  nw, 
		  nd, 
		  0.0);

    for(id=0;id<nd;id++) {
	dwt[id]=1;
    }
    sf_banded_solve(spl,dwt);  

    for(iz=0;iz<nz;iz++) {
	sf_warning("iz=%d of %d",iz+1,nz);

	sf_floatread(dat,nd,Fd);

	sf_banded_solve(spl,dat);  
	sf_int1_lop( true,   /* adj */
		     false,  /* add */
		     nh,     /* n model */
		     nd,     /* n data */
		     mod,   
		     dat);

	sf_floatwrite(mod,nh,Fm);
    }

/*------------------------------------------------------------*/

    sf_int1_close();

    free(map);
    free(mod);
    free(dat);
    free(mwt);
    free(dwt);

    exit(0);
}


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
    axa  az,ah,ahx,ahy,ahz,aj;
    int  iz,   ihx,ihy,ihz;
    float       hx, hy, hz;

    sf_file Fd; /*  data =   vector offset (hx,hy,hz)-z */
    sf_file Fm; /* model = absolute offset     h     -z */

    int nw;    /* spline order */
    int nd,id; /*  data size (nd=nhx*nhy*nhz) */
    int nm,im; /* model size (nm=nh) */

    float *dat=NULL;
    float *mod=NULL;
    float *map=NULL;

    float *mwt=NULL;
    float *dwt=NULL;

    int i;

/*------------------------------------------------------------*/

    sf_init(argc,argv);

    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getint(   "nw",&nw))     nw=4;     /* spline order */

    Fd = sf_input ("in");
    iaxa(Fd,&ahx,1); ahx.l="hx"; if(verb) raxa(ahx);
    iaxa(Fd,&ahy,2); ahy.l="hy"; if(verb) raxa(ahy);
    iaxa(Fd,&ahz,3); ahz.l="hz"; if(verb) raxa(ahz);
    iaxa(Fd,&az, 4);  az.l="z";  if(verb) raxa(az);

    if(!sf_getint  ("nh",&ah.n)) ah.n=ahx.n + ahx.o/ahx.d;
    if(!sf_getfloat("oh",&ah.o)) ah.o=0;
    if(!sf_getfloat("dh",&ah.d)) ah.d=ahx.d;
    ah.l="h";
    if(verb) raxa(ah);

    aj.n=1; aj.o=0; aj.d=1.; aj.l="";

    Fm = sf_output("out");
    oaxa(Fm,&ah,1);
    oaxa(Fm,&az,2);
    oaxa(Fm,&aj,3);
    oaxa(Fm,&aj,4);

/*------------------------------------------------------------*/
    nm = ah.n;               /* model size */
    nd = ahx.n*ahy.n*ahz.n;  /*  data size */

    map = sf_floatalloc(nd); /* mapping */

    mod = sf_floatalloc(nm); /* model vector */
    dat = sf_floatalloc(nd); /*  data vector */

    mwt = sf_floatalloc(nm); /* model weight */
    dwt = sf_floatalloc(nd); /*  data weight */


    sf_prefilter_init(nw,    /* spline order */
		      nd,    /* temporary storage */
		      2*nd); /* padding */

    for(ihz=0;ihz<ahz.n;ihz++) {
	hz = ahz.o + ihz * ahz.d;         hz*=hz;
	for(ihy=0;ihy<ahy.n;ihy++) {
	    hy = ahy.o + ihy * ahy.d;     hy*=hy;
	    for(ihx=0;ihx<ahx.n;ihx++) {
		hx = ahx.o + ihx * ahx.d; hx*=hx;
		
		i = ihz * ahx.n*ahy.n + 
		    ihy * ahx.n       +
		    ihx;
		
		map[i] = sqrtf(hx+hy+hz);
	    }
	}
    }

    sf_int1_init( map, 
		  ah.o, ah.d, ah.n, 
		  sf_spline_int, 
		  nw, 
		  nd);

    for(id=0;id<nd;id++) {
	dwt[id]=1;
    }
    sf_prefilter_apply(nd,dwt);  
    sf_int1_lop(true,false,nm,nd,mwt,dwt);

    for(iz=0;iz<az.n;iz++) {
	sf_warning("iz=%d of %d",iz+1,az.n);

	sf_floatread(dat,nd,Fd);

	sf_prefilter_apply(nd,dat);  
	sf_int1_lop( true,   // adj
		     false,  // add
		     nm,     // n model
		     nd,     // n data
		     mod,   
		     dat);

/*	for(im=0;im<nm;im++) {*/
/*	    if(mwt[im] >0) mod[im] /= mwt[im];*/
/*	}*/

	sf_floatwrite(mod,nm,Fm);
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


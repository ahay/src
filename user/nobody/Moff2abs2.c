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
    axa  az,ah,aa,ahx,ahz,aj;
    float    h, a, hx, hz;
    int  iz,ih,ia;

    sf_file Fd; /*  data =   vector offset (hx,hz)-z */
    sf_file Fm; /* model = absolute offset ( h, a)-z */

    int nw;    /* spline order */
    int nd,id; /*  data size (nd=nh *na ) */
    int nm;    /* model size (nm=nhx*nhz) */

    float  *dat=NULL;
    float  *mod=NULL;
    float **map=NULL;

/*------------------------------------------------------------*/

    sf_init(argc,argv);

    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getint(   "nw",&nw))     nw=4;     /* spline order */

    Fm = sf_input ("in");
    iaxa(Fm,&ahx,1); ahx.l="hx"; if(verb) raxa(ahx);
    iaxa(Fm,&ahz,2); ahz.l="hz"; if(verb) raxa(ahz);
    iaxa(Fm,&az, 3);  az.l="z";  if(verb) raxa(az);

    if(!sf_getint  ("nh",&ah.n)) ah.n=ahx.n + ahx.o/ahx.d;
    if(!sf_getfloat("oh",&ah.o)) ah.o=0;
    if(!sf_getfloat("dh",&ah.d)) ah.d=ahx.d;
    ah.l="h";
    if(verb) raxa(ah);

    aa.n = 180;
    aa.d = 2;
    aa.o = 0;
    aa.l ="a";
    if(verb) raxa(aa);

    aj.n=1; aj.o=0; aj.d=1.; aj.l="";

    Fd = sf_output("out");
    oaxa(Fd,&ah,1);
    oaxa(Fd,&aa,2);
    oaxa(Fd,&az,3);
    oaxa(Fd,&aj,4);
    oaxa(Fd,&aj,5);

/*------------------------------------------------------------*/
    nm = ahx.n*ahz.n;  /* model size */
    nd =  ah.n* aa.n;  /*  data size */
    sf_warning("nm=%d nd=%d",nm,nd);

    map = sf_floatalloc2(2,nd); /* mapping */

    mod = sf_floatalloc(nm); /* model vector */
    dat = sf_floatalloc(nd); /*  data vector */

    id=0;
    for(ia=0;ia<aa.n;ia++) {
	a = aa.o + ia * aa.d;
	a*= SF_PI/180;
	
	for(ih=0;ih<ah.n;ih++) {
	    h = ah.o + ih * ah.d;
	    
	    hx = h * cos(a);
	    hz = h * sin(a);

	    map[id][0] = hx;
	    map[id][1] = hz;

	    id++;
	}
    }

    sf_prefilter_init(nw,    /* spline order */
		      nd,    /* temporary storage */
		      2*nd); /* padding */

    sf_int2_init( map, 
		  ahx.o, ahz.o,
		  ahx.d, ahz.d,
		  ahx.n, ahz.n,
		  sf_spline_int, nw, 
		  nd);

    for(iz=0;iz<az.n;iz++) {
	sf_warning("iz=%d of %d",iz+1,az.n);

	sf_floatread (mod,nm,Fm);

	sf_prefilter_apply(nd,dat);
	sf_int2_lop(false,false,nm,nd,mod,dat);

	sf_floatwrite(dat,nd,Fd);
    }

/*------------------------------------------------------------*/

    free(*map); free(map);
    free(mod);
    free(dat);

    exit(0);
}


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
    sf_axis  az,ah,aa,ab,ahx,ahy,ahz;
    sf_bands spl;
    float    h, a, b, hx, hy, hz;
    int  iz,ih,ia,ib;
    int  nz,nh,na,nb,nhx,nhy,nhz;
    float   dh,da,db,dhx,dhy,dhz;
    float   oh,oa,ob,ohx,ohy,ohz;

    sf_file Fd; /*  data =   vector offset (hx,hy,hz)-z */
    sf_file Fm; /* model = absolute offset ( h, a, b)-z */

    int nw;    /* spline order */
    int nd,id; /*  data size (nd=nh *na *nb ) */
    int nm;    /* model size (nm=nhx*nhy*nhz) */

    float  *dat=NULL;
    float  *mod=NULL;
    float **map=NULL;

/*------------------------------------------------------------*/

    sf_init(argc,argv);

    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getint(   "nw",&nw))     nw=4;     /* spline order */

    Fm = sf_input ("in");
    ahx = sf_iaxa(Fm,1); sf_setlabel(ahx,"hx"); if(verb) sf_raxa(ahx);
    ahy = sf_iaxa(Fm,2); sf_setlabel(ahy,"hy"); if(verb) sf_raxa(ahy);
    ahz = sf_iaxa(Fm,3); sf_setlabel(ahz,"hz"); if(verb) sf_raxa(ahz);
    az  = sf_iaxa(Fm,4); sf_setlabel(az  ,"z"); if(verb) sf_raxa(az);

    nhx = sf_n(ahx); ohx = sf_o(ahx); dhx = sf_d(ahx);
    nhy = sf_n(ahy); ohy = sf_o(ahy); dhy = sf_d(ahy);
    nhz = sf_n(ahz); ohz = sf_o(ahz); dhz = sf_d(ahz);
    nz  = sf_n(az);

    /* h = absolute offset */
    if(!sf_getint  ("nh",&nh)) nh=nhx + ohx/dhx;
    if(!sf_getfloat("oh",&oh)) oh=0;
    if(!sf_getfloat("dh",&dh)) dh=dhx;
    ah = sf_maxa(nh,oh,dh); sf_setlabel(ah,"h"); if(verb) sf_raxa(ah);
 
    /* a = angle in x-z plane */
    if(!sf_getint  ("na",&na)) na=180; 
    if(nhz==1) na=1;
    if(!sf_getfloat("oa",&oa)) oa=0.;
    if(!sf_getfloat("da",&da)) da=2.;
    aa = sf_maxa(na,oa,da); sf_setlabel(aa,"a"); if(verb) sf_raxa(aa);

    /* b = angle in x-y plane */
    if(!sf_getint  ("nb",&nb)) nb=180; 
    if(nhy==1) nb=1;
    if(!sf_getfloat("ob",&ob)) ob=0.;
    if(!sf_getfloat("db",&db)) db=2.;
    ab = sf_maxa(nb,ob,db); sf_setlabel(ab,"b"); if(verb) sf_raxa(ab);

    Fd = sf_output("out");
    sf_oaxa(Fd,ah,1);
    sf_oaxa(Fd,ab,2);
    sf_oaxa(Fd,aa,3);
    sf_oaxa(Fd,az,4);
    sf_putint(Fd,"n5",1);

/*------------------------------------------------------------*/
    nm = nhx*nhy*nhz;  /* model size */
    nd =  nh* na* nb;  /*  data size */
    sf_warning("nm=%d nd=%d",nm,nd);

    map = sf_floatalloc2(3,nd); /* mapping */

    mod = sf_floatalloc(nm); /* model vector */
    dat = sf_floatalloc(nd); /*  data vector */

    id=0;
    for(ia=0;ia<na;ia++) {
	a = oa + ia * da;
	a*= SF_PI/180;

	for(ib=0;ib<nb;ib++) {
	    b = ob + ib * db;
	    b*= SF_PI/180;
	    	    
	    for(ih=0;ih<nh;ih++) {
		h = oh + ih * dh;
		
		hx = h * cos(a) * cos(b);
		hy = h * cos(a) * sin(b);
		hz = h * sin(a);
		
		map[id][0] = hx;
		map[id][1] = hy;
		map[id][2] = hz;
		
		id++;
	    }
	}
    }
    
    spl = sf_spline_init(nw,nd);

    sf_int3_init( map, 
		  ohx, ohy, ohz,
		  dhx, dhy, dhz,
		  nhx, nhy, nhz,
		  sf_spline_int, nw, 
		  nd);

    for(iz=0;iz<nz;iz++) {
	sf_warning("iz=%d of %d",iz+1,nz);

	sf_floatread (mod,nm,Fm);

	sf_banded_solve(spl,dat);
	sf_int3_lop(false,false,nm,nd,mod,dat);

	sf_floatwrite(dat,nd,Fd);
    }

/*------------------------------------------------------------*/

    free(*map); free(map);
    free(mod);
    free(dat);

    exit(0);
}


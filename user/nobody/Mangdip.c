/* Dip correction for kh/km from kh/kz
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
    axa  az,aa,ax;
    int  iz,ia,ix;

    sf_file Fd; /*  data = angle z-a (input)  */
    sf_file Fm; /* model = angle z-a (output) */
    sf_file Fa;
    
    int nw;    /* spline order */

    float **mmm=NULL;
    float **ddd=NULL;

    float  *dat=NULL;
    float  *mod=NULL;
    float  *map=NULL;
    float  *dip=NULL;

    float d,t;

/*------------------------------------------------------------*/

    sf_init(argc,argv);

    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getint(   "nw",&nw))     nw=4;     /* spline order */

    Fd = sf_input ("in");
    iaxa(Fd,&az,1); az.l="h"; if(verb) raxa(az);
    iaxa(Fd,&aa,2); aa.l="a"; if(verb) raxa(aa);
    iaxa(Fd,&ax,3); ax.l="x"; if(verb) raxa(ax);

    Fm = sf_output("out");
    oaxa(Fm,&az,1);
    oaxa(Fm,&aa,2);
    oaxa(Fm,&ax,3);

    Fa = sf_input("dip");
    
/*------------------------------------------------------------*/
    
    mmm = sf_floatalloc2(az.n,aa.n);
    ddd = sf_floatalloc2(az.n,aa.n);

    map = sf_floatalloc(aa.n); /* mapping */
    mod = sf_floatalloc(aa.n); /* model vector */
    dat = sf_floatalloc(aa.n); /*  data vector */
    dip = sf_floatalloc(az.n);      /*   dip vector */

    sf_prefilter_init(nw,      /* spline order */
		      aa.n,    /* temporary storage */
		      2*aa.n); /* padding */

    sf_int1_init( map, 
		  aa.o, aa.d, aa.n, 
		  sf_spline_int, 
		  nw, 
		  aa.n);

    for(ix=0;ix<ax.n;ix++) {
	sf_floatread(dip   ,az.n     ,Fa);
	sf_floatread(ddd[0],az.n*aa.n,Fd);

	for(iz=0;iz<az.n;iz++) {
	    sf_warning("iz=%d of %d",iz+1,az.n);

	    d = dip[iz];
	    d = 0;
	    d*= d;

	    for(ia=0;ia<aa.n;ia++){
		dat[ia] = ddd[ia][iz];

		t = aa.o + ia*aa.d;
		map[ia] = t / sqrt(1+d);
	    }
	    
	    sf_prefilter_apply( aa.n,
				dat);  

	    sf_int1_init(map, 
			 aa.o, aa.d, aa.n, 
			 sf_spline_int, 
			 nw, aa.n);

	    sf_int1_lop( true,   /* adj */
			 false,  /* add */
			 aa.n,   /* n model */
			 aa.n,   /* n data */
			 mod,   
			 dat);

	    for(ia=0;ia<aa.n;ia++){
		mmm[ia][iz] = mod[ia];
	    }
	}
	
	sf_floatwrite(mmm[0],az.n*aa.n,Fm);
    }

/*------------------------------------------------------------*/

    sf_int1_close();

    free(map);
    free(mod);
    free(dat);
    free(dip);
    free(*ddd); free(ddd);
    free(*mmm); free(mmm);

    exit(0);
}


/* 3-D imaging conditions for shot-profile WE migration */
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
#include "weimc.h"

int main (int argc, char *argv[])
{
    bool  verb;           /* verbosity */
    bool  hsym;           /* symmetric offset */
    char *itype;          /* imaging type (zero-offset, prestack) */
    axa az, ax, ay, aw,aj;
    axa ahz,ahx,ahy;      /* offset */

    sf_file Fus; /* source   wavefield D(x,y,z,w) */
    sf_file Fur; /* receiver wavefield U(x,y,z,w) */
    sf_file Fi;  /*    image           R(x,y,z)   */

    fslice imag,sdat,rdat;

    /*------------------------------------------------------------*/
    sf_init(argc,argv);

    if (!sf_getbool("verb",&verb)) verb =  true; /* verbosity flag */
    
    if (NULL == (itype = sf_getstring("itype"))) itype = "z";

    Fus = sf_input ( "in");
    Fur = sf_input ("rwf");
    Fi  = sf_output("out"); sf_settype(Fi,SF_FLOAT);

    iaxa(Fus,&ax,1); ax.l="x"; oaxa(Fi,&ax,1);
    iaxa(Fus,&ay,2); ay.l="y"; oaxa(Fi,&ay,2);    
    iaxa(Fus,&az,3); az.l="z"; oaxa(Fi,&az,3);
    iaxa(Fus,&aw,4); aw.l="w";
	    
    sdat = fslice_init( ax.n*ay.n*az.n,aw.n,sizeof(float complex));
    rdat = fslice_init( ax.n*ay.n*az.n,aw.n,sizeof(float complex));

    switch(itype[0]) {
	case 'p':
	    sf_warning("Full offset IC");

	    if(!sf_getint("nhx",&ahx.n)) ahx.n=1;
	    if(!sf_getint("nhy",&ahy.n)) ahy.n=1;
	    if(!sf_getint("nhz",&ahz.n)) ahz.n=1;

	    ahx.o=0;      ahy.o=0;      ahz.o=0;
	    ahx.d=ax.d;   ahy.d=ay.d;   ahz.d=az.d;
	    ahx.l="hx";   ahy.l="hy";   ahz.l="hz";

	    if(!sf_getbool("hsym",&hsym)) hsym = false;
	    if(hsym) {
		if(ahx.n>1) { ahx.o = - ahx.n * ahx.d; ahx.n *=2; }
		if(ahy.n>1) { ahy.o = - ahy.n * ahy.d; ahy.n *=2; }
		if(ahz.n>1) { ahz.o = - ahz.n * ahz.d; ahz.n *=2; }
	    }
	    
	    oaxa(Fi,&ahx,4);
	    oaxa(Fi,&ahy,5);
	    oaxa(Fi,&ahz,6);

	    imag = fslice_init( ax.n* ay.n* az.n,
			       ahx.n*ahy.n*ahz.n,sizeof(float));
	    break;
	case 'z':
	default:
	    sf_warning("Zero offset IC");
	    aj.n=1; aj.o=0; aj.d=1; aj.l=" ";
	    oaxa(Fi,&aj,4);

	    imag = fslice_init( ax.n*ay.n*az.n,1,sizeof(float));
	    break;
    }

    fslice_load(Fus,sdat,SF_COMPLEX);
    fslice_load(Fur,rdat,SF_COMPLEX);
    /*------------------------------------------------------------*/

    weimc_init(verb,ax,ay,az,aw);

    switch(itype[0]) {
	case 'p':
	    hhimc_init(ahx,ahy,ahz);
	    hhimc(sdat,rdat,imag);
	    break;
	case 'z':
	default:
	    zoimc(sdat,rdat,imag);
	    break;
    }

    weimc_close();

    /*------------------------------------------------------------*/    
    /* slice management (temp files) */
    fslice_dump(Fi,imag,SF_FLOAT);

    fslice_close(sdat);
    fslice_close(rdat);
    fslice_close(imag);

    exit (0);
}

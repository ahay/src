/* */
/*
  Copyright (C) 2006 Colorado School of Mines
  
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
    float threshold;

    sf_axis  az,ax, ak,aa;
    int      iz,ix, ik;
    float     z, x;

    int      nn,kk;

    sf_file Fi;
    sf_file Fo;

    float **map=NULL;
    float  *xxx=NULL;
    float  *zzz=NULL;
    float  *rrr=NULL;

/*------------------------------------------------------------*/

    sf_init(argc,argv);

    if(! sf_getbool("verb",&verb)) verb=false;             /* verbosity flag */
    if(! sf_getfloat("threshold",&threshold)) threshold=0; /* pick threshold */

    Fi = sf_input ("in");
    az = sf_iaxa(Fi,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az);
    ax = sf_iaxa(Fi,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax);

/*------------------------------------------------------------*/
    map = sf_floatalloc2( sf_n(az), sf_n(ax) );
    nn = sf_n(az)* sf_n(ax);
    sf_floatread(map[0],nn,Fi);
/*------------------------------------------------------------*/
/* count */
    ik = 0;
    for (ix=0; ix<sf_n(ax); ix++) {
	for (iz=0; iz<sf_n(az); iz++) {
	    if(SF_ABS( map[ix][iz]) >= threshold) ik++;
	}
    }
    kk = ik;
    
/*------------------------------------------------------------*/
    aa = sf_maxa( 3,0.,1.); if(verb) sf_raxa(aa);
    ak = sf_maxa(kk,0.,1.); if(verb) sf_raxa(ak);

    Fo = sf_output("out");
    sf_oaxa(Fo,ak,1);
    sf_oaxa(Fo,aa,2);
    sf_putint(Fo,"n3",1);

/*------------------------------------------------------------*/
    xxx = sf_floatalloc(kk);
    zzz = sf_floatalloc(kk);
    rrr = sf_floatalloc(kk);

    ik = 0;
    for (ix=0; ix<sf_n(ax); ix++) {
	x = sf_o(ax) + ix * sf_d(ax);

	for (iz=0; iz<sf_n(az); iz++) {
	    z = sf_o(az) + iz * sf_d(az);
	    
	    /* -> test threshold */
	    if(SF_ABS( map[ix][iz]) >= threshold) {
		ik++;
		
		xxx[ik] = x;
		zzz[ik] = z;
		rrr[ik] = map[ix][iz];

	    }
	    /* -> end test threshold */
	}
    }

/*------------------------------------------------------------*/
    
    sf_floatwrite(xxx,kk,Fo);
    sf_floatwrite(zzz,kk,Fo);
    sf_floatwrite(rrr,kk,Fo);

/*------------------------------------------------------------*/

    free(*map); free(map);

    exit(0);
}


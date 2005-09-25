/* Fill holes (nearest-neighbor)
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
    int   fill; /* hole-filling parameter */
    float feps;

    axa az,ax,ay; /* Cartesian coordinates */
    int iz,ix,iy;
    int kz,kx,ky;
    int lz,lx,ly;
    int jz,jx,jy;

    int k,n;
    float v;

    sf_file Fh;
    sf_file Ff;

    float ***hh=NULL;
    
    sf_init(argc,argv);

    if(! sf_getbool (    "verb",&verb    ))     verb=false;
    if(! sf_getint  (    "fill",&fill    ))     fill=0;
    if(! sf_getfloat(    "feps",&feps    ))     feps=0.;

    Fh = sf_input ("in");
    iaxa(Fh,&az,1); az.l="z"; if(verb) raxa(az);
    iaxa(Fh,&ax,2); ax.l="x"; if(verb) raxa(ax);
    iaxa(Fh,&ay,3); ay.l="y"; if(verb) raxa(ay);

    Ff = sf_output("out");
    oaxa(Ff,&az,1);
    oaxa(Ff,&ax,2);
    oaxa(Ff,&ay,3);
    
    hh=sf_floatalloc3(az.n,ax.n,ay.n); 
    sf_floatread(hh[0][0],az.n*ax.n*ay.n,Fh);

    /* fill holes */
    n = fill;
    
    for(iy=0;iy<ay.n;iy++) {
	for(ix=0;ix<ax.n;ix++) {
	    for(iz=0;iz<az.n;iz++) {

		if( SF_ABS(hh[iy][ix][iz]) <= feps) {
		    
		    kx=SF_MIN(SF_MAX(ix-n,0),ax.n-1);
		    ky=SF_MIN(SF_MAX(iy-n,0),ay.n-1);
		    kz=SF_MIN(SF_MAX(iz-n,0),az.n-1);
		    
		    lx=SF_MIN(SF_MAX(ix+n,0),ax.n-1);
		    ly=SF_MIN(SF_MAX(iy+n,0),ay.n-1);
		    lz=SF_MIN(SF_MAX(iz+n,0),az.n-1);
		    
		    k=0;
		    v=0.;
		    for(jy=ky;jy<=ly;jy++) {
			for(jx=kx;jx<=lx;jx++) {
			    for(jz=kz;jz<=lz;jz++) {
				
				if(SF_ABS(hh[jy][jx][jz]) > feps) {
				    v +=  hh[jy][jx][jz];
				    k++;
				}
			    }
			}
		    } /* local loop */
		    
		    if(k>0) hh[iy][ix][iz] = v/k;
		} /* end if missing */
		
	    }
	}
    }  /* global loop */

    sf_floatwrite(hh[0][0],az.n*ax.n*ay.n,Ff);

    free(**hh); free(*hh); free(hh);

    exit (0);
}

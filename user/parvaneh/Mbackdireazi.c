/*Background directivity(Azimuth). */
/*
Copyright (C) 2011 University of Texas at Austin

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

int main (int argc, char *argv[])
{
    /*differentiation coefficient*/
    float c0=-1./12., c1=+2./3.;
    sf_file Zo=NULL,Zz=NULL;
    sf_axis at,ay,ax;
    int it,iy,ix;
    int nt,ny,nx;
    float idt,idy,dt,dy;
   
    
    float ***tz,***azimuth,***tzt,***tzy;
    
    /*sf_init(!sf_getbool("verb",&verb)) verb=0;*/
    
    sf_init(argc,argv);
    Zz=sf_input("in");
    Zo=sf_output("out");
    
    at=sf_iaxa(Zz,1); nt=sf_n(at); dt=sf_d(at);
    ay=sf_iaxa(Zz,2); ny=sf_n(ay); dy=sf_d(ay);
    ax=sf_iaxa(Zz,3); nx=sf_n(ax); 

    sf_oaxa(Zo,at,1);
    sf_oaxa(Zo,ay,2);
    sf_oaxa(Zo,ax,3);

    idt=1/dt;
    idy=1/dy;

    tz=sf_floatalloc3(nt,ny,nx); sf_floatread(tz[0][0],nt*ny*nx,Zz);

    azimuth=sf_floatalloc3(nt,ny,nx);
    tzt=sf_floatalloc3(nt,ny,nx);
    tzy=sf_floatalloc3(nt,ny,nx);

    for (it=0; it<nt; it++) {
	for (iy=0; iy<ny; iy++) {
	    for (ix=0; ix<nx; ix++) {
		 azimuth[ix][iy][it]=0;
	    }
	}
    }

    for (it=2; it<nt-2; it++) {
	for (iy=2; iy<ny-2; iy++) {
	    for (ix=2; ix<nx-2; ix++) {
		 tzt[ix][iy][it]=
		     idt*(c0*(tz[ix][iy][it+2]-tz[ix][iy][it-2])+
			  c1*(tz[ix][iy][it+1]-tz[ix][iy][it-1]));

		 tzy[ix][iy][it]=
		     idy*(c0*(tz[ix][iy+2][it]-tz[ix][iy-2][it])+
			  c1*(tz[ix][iy+1][it]-tz[ix][iy-1][it]));
		 
		 azimuth[ix][iy][it]=(atan(-(tzy[ix][iy][it]/tzt[ix][iy][it]))*180)/3.1415;
	    }
	}
    }

    sf_floatwrite(azimuth[0][0],nt*ny*nx,Zo);
    
    sf_close();

    exit (0);
}    

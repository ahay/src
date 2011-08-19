/*Background directivity(Dip). */
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
    bool verb;  /*verbose flag*/
    sf_file Zi=NULL,Zo=NULL,Zz=NULL;
    sf_axis at,ay,ax;
    int it,iy,ix;
    int nt,ny,nx;
    float idt,idx,dx,dt,dy;
    
    float ***tz,***tt,***dip,***tzt,***tzx;
    
    /*sf_init(!sf_getbool("verb",&verb)) verb=0;*/

    sf_init(argc,argv);
    Zi=sf_input("ti");
    Zz=sf_input("in");
    Zo=sf_output("out");
    
    at=sf_iaxa(Zi,1); nt=sf_n(at); dt=sf_d(at);
    ay=sf_iaxa(Zi,2); ny=sf_n(ay); dy=sf_d(ay);
    ax=sf_iaxa(Zi,3); nx=sf_n(ax); dx=sf_d(ax);

    sf_oaxa(Zo,at,1);
    sf_oaxa(Zo,ay,2);
    sf_oaxa(Zo,ax,3);

    idt=1/dt;
    idx=1/dx;

    tt=sf_floatalloc3(nt,ny,nx); sf_floatread(tt[0][0],nt*ny*nx,Zi);
    tz=sf_floatalloc3(nt,ny,nx); sf_floatread(tt[0][0],nt*ny*nx,Zz);

    dip=sf_floatalloc3(nt,ny,nx);
    tzt=sf_floatalloc3(nt,ny,nx);
    tzx=sf_floatalloc3(nt,ny,nx);

   for (it=2; it<nt-2; it++) {
	for (iy=2; iy<ny-2; iy++) {
	    for (ix=2; ix<nx-2; ix++) {
		 dip[ix][iy][it]=0;
	    }
	}
    }

    for (it=0; it<nt; it++) {
	for (iy=0; iy<ny; iy++) {
	    for (ix=0; ix<nx; ix++) {
		 tzt[ix][iy][it]=
		     idt*(c0*(tz[ix][iy][it+2]-tz[ix][iy][it-2])+
			  c1*(tz[ix][iy][it+1]-tz[ix][iy][it-1]));

		 tzx[ix][iy][it]=
		     idx*(c0*(tz[ix+2][iy][it]-tz[ix-2][iy][it])+
			  c1*(tz[ix+1][iy][it]-tz[ix-1][iy][it]));
		 
		 dip[ix][iy][it]=-(tzx[ix][iy][it]/tzt[ix][iy][it]);
	    }
	}
    }

    sf_floatwrite(dip[0][0],nt*ny*nx,Zo);

    sf_close();

    exit (0);
}    

    
		 
	       

/* Directional angle transform for 3-D time image cube I(x,z,t) into G(x,z,d) */
/*
  Copyright (C) 2008 University of Texas at Austin

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

#include <math.h>
#include <rsf.h>


int main (int argc, char *argv[])
{
    bool verb;

    float d,x,z,t,dd,d0,dx,x0,dz,z0,dt;
    int   nd,nx,nz,nt;

    sf_axis ax,az,at,ad;
    int ix,iz,it,id,ixm,izm;

    float xp,zp,imt1,imt2,imt3,imt4,tx,tz;

    float   ***imgt, ***imgd, **slow;
    sf_file   Fimgt,   Fimgd,  slowness;

    sf_init (argc,argv);

    /*------------------------------------------------------------*/

    Fimgt = sf_input("in");
    Fimgd = sf_output("out");

    if (SF_FLOAT != sf_gettype(Fimgt)) sf_error("Need float input");

    ax=sf_iaxa(Fimgt,1); nx=sf_n(ax); x0=sf_o(ax); dx=sf_d(ax);
    az=sf_iaxa(Fimgt,2); nz=sf_n(az); z0=sf_o(az); dz=sf_d(az);
    at=sf_iaxa(Fimgt,3); nt=sf_n(at); dt=sf_d(at);

    if (!sf_getint  ("nd",&nd)) nd=nt;       
    if (!sf_getfloat("dd",&dd)) dd=160.0/(nt-1);
    if (!sf_getfloat("d0",&d0)) d0=-80.0;
    ad = sf_maxa(nd,d0,dd);
    sf_oaxa(Fimgd,ad,3);

    if (!sf_getbool("verb",&verb)) verb=false; /* verbosity flag */

    imgt = sf_floatalloc3(nx,nz,nt);
    imgd = sf_floatalloc3(nx,nz,nd);

    sf_floatread(imgt[0][0],nx*nz*nt,Fimgt);

    slowness = sf_input("slowness");
    slow = sf_floatalloc2(nx,nz);

    sf_floatread(slow[0],nx*nz,slowness);


    /*------------------------------------------------------------*/

   for (id = 0; id < nd; id++) {

        d = d0 + id*dd;
        if (verb) sf_warning("idip=%d",id);

        for (ix = 0; ix < nx; ix++) {

            x = x0 + ix*dx;

            for (iz = 0; iz < nz; iz++) {

                z = z0 + iz*dz;
                imgd[ix][iz][id] = 0.0;           

                for (it = 0; it < nt; it++) {

                    t = it*dt;
                    xp = x - t*sinf(d/180*SF_PI)/(slow[ix][iz]);
                    zp = z + t*cosf(d/180*SF_PI)/(slow[ix][iz]);

                    /* Bilinear interpolation */ 

                    ixm = floor((xp-x0)/dx);
                    izm = floor((zp-z0)/dz);

                    if (ixm >= 0 && izm >= 0 && ixm <= (nx-1) && izm <= (nz-1)){

                       tx = (xp-ixm*dx-x0)/dx;
                       tz = (zp-izm*dz-z0)/dz;

                       imt1 = imgt[ixm][izm][it]*(1.0-tx)*(1.0-tz);
                       imt2 = imgt[ixm+1][izm][it]*tx*(1.0-tz);
                       imt3 = imgt[ixm+1][izm+1][it]*tx*tz;
                       imt4 = imgt[ixm][izm+1][it]*(1.0-tx)*tz;

                       imgd[ix][iz][id] += (imt1 + imt2 + imt3 + imt4);
		    }

		}
            }
	}

    }

    sf_floatwrite(imgd[0][0],nx*nz*nd,Fimgd);

    exit (0);
}

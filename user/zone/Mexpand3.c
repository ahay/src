/* Expand 3D data  */
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

#include <rsf.h>

int expandb(float ***vexpmodel, float ***vmodel, int vnx, int vny, int vnz/*orignal size*/, int vlx, int vrx, int vly, int vry, int vt, int vb)
/*<expand model>*/
{
    int ix, iy, iz;
    
    if (vnx<0 || vny<0 || vnz<0 || vlx<0 || vrx<0 || vly<0 || vry<0 || vt<0 || vb<0) {
	sf_error("Cann't expand model!");	
    }
    
    if ( vexpmodel==NULL || vmodel==NULL) {
	sf_error("Need allocate memory first, expand model fault!");
    }

    // z direction
	for (iy=vly; iy<vly+vny; iy++) {
        for (ix=vlx; ix<vlx+vnx; ix++) {
        	for (iz=vt; iz<vt+vnz; iz++) {
        	    vexpmodel[iy][ix][iz] = vmodel[iy-vly][ix-vlx][iz-vt];
        	}
        	for (iz=0; iz<vt; iz++) {
        	    vexpmodel[iy][ix][iz] = vmodel[iy-vly][ix-vlx][0];
        	}
        	for (iz=vt+vnz; iz<vt+vnz+vb; iz++) {
        	    vexpmodel[iy][ix][iz] = vmodel[iy-vly][ix-vlx][vnz-1];
        	}
        }
    }
    

    // x direction
    for (iy=vly; iy<vly+vny;iy++) {
        for (ix=0; ix<vlx; ix++){
        	for (iz=0; iz<vt+vnz+vb; iz++) {
        	    vexpmodel[iy][ix][iz] = vexpmodel[iy][vlx][iz];
        	}
        }
        for (ix=vlx+vnx; ix<vlx+vnx+vrx; ix++) {
        	for (iz=0; iz<vt+vnz+vb; iz++) {
        	    vexpmodel[iy][ix][iz] = vexpmodel[iy][vlx+vnx-1][iz];
        	}
        }
    }
    
    
    
    // y direction
    for (iy=0; iy<vly;iy++) {
        for (ix=0; ix<vlx+vnx+vrx; ix++){
        	for (iz=0; iz<vt+vnz+vb; iz++) {
        	    vexpmodel[iy][ix][iz] = vexpmodel[vly][ix][iz];
        	}
        }
    }
    for (iy=vly+vny; iy<vly+vny+vry;iy++) {
        for (ix=0; ix<vlx+vnx+vrx; ix++){
        	for (iz=0; iz<vt+vnz+vb; iz++) {
        	    vexpmodel[iy][ix][iz] = vexpmodel[vly+vny-1][ix][iz];
        	}
        }
    }
    
    return 1;
 }


int main(int argc, char* argv[])
{
    int nx, nz, ny, leftx, rightx, lefty, righty, top, bottom;
    int lnx, lny, lnz;
    int ix,iy;
    float dx,dy,dz, ox,oy,oz;
    sf_file out, in;
    
    float ***exmodel, ***origmodel;
    
    sf_init(argc, argv);
    out = sf_output("out");
    in  = sf_input("in");
    
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input!");
    if (!sf_histint(in, "n1", &nz)) sf_error("No n1 in input!");
    if (!sf_histint(in, "n2", &nx)) sf_error("No n2 in input!");
    if (!sf_histint(in, "n3", &ny)) sf_error("No n3 in input!");

    if (!sf_histfloat(in,"d1",&dz)) dz = 10.0;
    if (!sf_histfloat(in,"d2",&dx)) dx = 10.0;
    if (!sf_histfloat(in,"d3",&dy)) dy = 10.0;
    
    if (!sf_histfloat(in,"o1",&oz)) oz = 0.0;
    if (!sf_histfloat(in,"o2",&ox)) ox = 0.0;
    if (!sf_histfloat(in,"o3",&oy)) oy = 0.0;

    if (!sf_getint("leftx", &leftx)) leftx = 0.5*nx;
    if (!sf_getint("rightx", &rightx)) rightx = 0.5*nx;
    if (!sf_getint("lefty", &lefty)) lefty = 0.5*ny;
    if (!sf_getint("righty", &righty)) righty = 0.5*ny;
    if (!sf_getint("top", &top)) top = 0;
    if (!sf_getint("bottom", &bottom)) bottom = 0;

    lnx = nx+leftx+rightx;
    lny = ny+lefty+righty;
    lnz = nz+top+bottom;

    sf_putint(out, "n1", lnz);
    sf_putint(out, "n2", lnx);
    sf_putint(out, "n3", lny);
    sf_putfloat(out,"d1",dz);
    sf_putfloat(out, "d2", dx);
    sf_putfloat(out, "d3", dy);
    sf_putfloat(out,"o1",oz-top*dz);
    sf_putfloat(out, "o2", ox-leftx*dx);
    sf_putfloat(out, "o3", oy-lefty*dy);


    origmodel = sf_floatalloc3(nz, nx, ny);
    for (iy=0; iy<ny; iy++) {
    for (ix=0; ix<nx; ix++) {
	sf_floatread(origmodel[iy][ix], nz, in);
    }
    }

    exmodel = sf_floatalloc3(lnz,lnx,lny);    
    expandb(exmodel, origmodel, nx, ny, nz, leftx, rightx, lefty, righty, top, bottom);

    sf_warning(">>>nx=%d ny=%d  nz=%d leftx=%d rightx=%d lefty=%d righty=%d top=%d bottom=%d lnx=%d lny=%d lnz=%d", nx, ny, nz, leftx, rightx, lefty, righty, top, bottom, lnx, lny, lnz);
    
    for (iy=0; iy<lny; iy++) {
        for(ix=0; ix<lnx; ix++){
        	sf_floatwrite(exmodel[iy][ix],lnz,out);
        }
    }
    
    free(**exmodel);
    free(*exmodel);
    free(exmodel);
    free(**origmodel);    
    free(*origmodel);
    free(origmodel);
    
    exit(0);
    
}



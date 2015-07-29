/* Expand 2D data  */
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


int expandb(float **vexpmodel, float **vmodel, int vnx, int vnz/*orignal size*/, int vl, int vr, int vt, int vb)
/*<expand model>*/
{
    int ix, iz;
    
    if (vnx<0 || vnz<0 || vl<0 || vr<0 || vt<0 || vb<0) {
	sf_error("Cann't expand model!");	
    }
    
    if ( vexpmodel==NULL || vmodel==NULL) {
	sf_error("Need allocate memory first, expand model fault!");
    }

	
    for (ix=vl; ix<vl+vnx; ix++) {
	for (iz=vt; iz<vt+vnz; iz++) {
	    vexpmodel[ix][iz] = vmodel[ix-vl][iz-vt];
	}
	for (iz=0; iz<vt; iz++) {
	    vexpmodel[ix][iz] = vexpmodel[ix][vt];
	}
	for (iz=vt+vnz; iz<vt+vnz+vb; iz++) {
	    vexpmodel[ix][iz] = vexpmodel[ix][vnz+vt-1];
	}
    }
   
    for (ix=0; ix<vl; ix++){
	for (iz=0; iz<vt+vnz+vb; iz++) {
	    vexpmodel[ix][iz] = vexpmodel[vl][iz];
	}
    }
    
    for (ix=vl+vnx; ix<vl+vnx+vr; ix++) {
	for (iz=0; iz<vt+vnz+vb; iz++) {
	    vexpmodel[ix][iz] = vexpmodel[vl+vnx-1][iz];
	}
    }
    return 1;
 }


int expandbc(float **vexpmodel, float **vmodel, int vnx, int vnz, int vbcout)
/*<Expand model for PML bandery condition>*/
{
    expandb(vexpmodel, vmodel, vnx, vnz, vbcout, vbcout, vbcout, vbcout);
    return 1;
}


int main(int argc, char* argv[])
{
    int nx, nz, left, right, top, bottom;
    int lnx, lnz;
    int ix;
    float dx,dz;
    sf_file out, in;
    
    float **exmodel, **origmodel;
    
    sf_init(argc, argv);
    out = sf_output("out");
    in  = sf_input("in");
    
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input!");
    if (!sf_histint(in, "n1", &nz)) sf_error("No n1 in input!");
    if (!sf_histint(in, "n2", &nx)) sf_error("No n2 in input!");
    if (!sf_histfloat(in,"d1",&dz)) dz = 10.0;
    if (!sf_histfloat(in,"d2",&dx)) dx = 10.0;

    if (!sf_getint("left", &left)) left = 0.5*nx;
    if (!sf_getint("right", &right)) right = 0.5*nx;
    if (!sf_getint("top", &top)) top = 0;
    if (!sf_getint("bottom", &bottom)) bottom = 0;

    lnx = nx+left+right;
    lnz = nz+top+bottom;

    sf_putint(out, "n1", lnz);
    sf_putint(out, "n2", lnx);
    sf_putfloat(out,"d1",dz);
    sf_putfloat(out, "d2", dx);

    origmodel = sf_floatalloc2(nz, nx);
    for (ix=0; ix<nx; ix++) {
	sf_floatread(origmodel[ix], nz, in);
    }

    exmodel = sf_floatalloc2(lnz,lnx);    
    expandb(exmodel, origmodel, nx, nz, left, right, top, bottom);

    sf_warning(">>>nx=%d nz=%d left=%d right=%d top=%d bottom=%d lnx=%d lnz=%d", nx, nz, left, right, top, bottom, lnx, lnz);
    
    for(ix=0; ix<lnx; ix++){
	sf_floatwrite(exmodel[ix],lnz,out);
    }

    sf_warning("nx=%d nz=%d left=%d right=%d top=%d bottom=%d lnx=%d lnz=%d", nx, nz, left, right, top, bottom, lnx, lnz);

    free(*exmodel);
    free(exmodel);
    free(*origmodel);
    free(origmodel);
    
    exit(0);
    
}



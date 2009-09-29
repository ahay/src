/* Maps a cube to a list, given a threshold (clip value). */
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

#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE
#endif

#include <rsf.h>

int main(int argc, char* argv[])
{
    bool verb;  /* verbosity flag */
    float clip; /* threshold (clip value) */
    sf_file Fc; /* cube file */
    sf_file Fl; /* list file */
    extern int fseeko(FILE *stream, off_t offset, int whence);

    sf_axis   ax,ay,az;
    sf_axis   aa;
    int   ix,iy,iz;
    int   nx,ny,nz,nj,na;
    int   nk=0,jk;
    
    float **cube;
    float dx,dy,dz;
    float x0,y0,z0;

    FILE* tfile;
    char* tname;
    float t2[3],t3[4];

    /*------------------------------------------------------------*/

    /* init RSF */
    sf_init(argc,argv);
    if(! sf_getbool("verb",&verb)) verb=false;
    if(! sf_getfloat("clip",&clip)) clip=0;

    Fc = sf_input ( "in"); /* input  cube */
    Fl = sf_output("out"); /* output list */
    
    /* read axes*/
    az=sf_iaxa(Fc,1); sf_setlabel(az,"z");
    ax=sf_iaxa(Fc,2); sf_setlabel(ax,"x");
    ay=sf_iaxa(Fc,3); sf_setlabel(ay,"y");

    nz=sf_n(az); z0=sf_o(az); dz=sf_d(az); 
    nx=sf_n(ax); x0=sf_o(ax); dx=sf_d(ax); 
    ny=sf_n(ay); y0=sf_o(ay); dy=sf_d(ay); 

    na=0; 
    if(ny>1) {
	if(verb) sf_warning("initiating 3D points");
	nj=4;
    } else {
	if(verb) sf_warning("initiating 2D points");
	nj=3;
    }
    /*------------------------------------------------------------*/

    cube = sf_floatalloc2(nz,nx);

    tfile = sf_tempfile(&(tname), "w+b");

    for (iy=0;iy<ny;iy++) {
/*	if(verb) sf_warning("iy=%d",iy);*/
	
	sf_floatread(cube[0],nz*nx,Fc);
	
	nk=0;
	for (ix=0; ix<nx; ix++) {
	    for (iz=0; iz<nz; iz++) {
		if( fabs(cube[ix][iz]) > clip) {
		    nk++;
		}
	    }
	}

	if(ny>1) {
	    jk=0;
	    for (ix=0; ix<nx; ix++) {
		for (iz=0; iz<nz; iz++) {
		    if( fabs(cube[ix][iz]) > clip) {
			t3[0] = x0 + ix * dx;
			t3[1] = y0 + iy * dy;
			t3[2] = z0 + iz * dz;
			t3[3] = cube[ix][iz];
			
			fseeko(tfile,jk*4*sizeof(float),SEEK_SET);
			fwrite(   t3,     sizeof(float),4,tfile);
			jk++;
		    }
		}
	    }
	} else {
	    jk=0;
	    for (ix=0; ix<nx; ix++) {
		for (iz=0; iz<nz; iz++) {
		    if( fabs(cube[ix][iz]) > clip) {
			t2[0] = x0 + ix * dx;
			t2[1] = z0 + iz * dz;
			t2[2] = cube[ix][iz];

			fseeko(tfile,jk*3*sizeof(float),SEEK_SET);
			fwrite(   t2,     sizeof(float),3,tfile);
			jk++;
		    }
		}
	    }
	} /* else ny=1 */

	na += nk;
    } /* iy */
    
    /* output axes */
    aa = sf_maxa(nj,0,1); sf_oaxa(Fl,aa,1); if(verb) sf_raxa(aa); free(aa);
    aa = sf_maxa(na,0,1); sf_oaxa(Fl,aa,2); if(verb) sf_raxa(aa); free(aa);

    if( ny>1) {
	for( jk=0; jk<nk; jk++) {
	    fseeko(tfile,jk*4*sizeof(float),SEEK_SET);
	    fread(    t3,     sizeof(float),4,tfile);
/*	    if(verb) sf_warning("%d, %g %g %g %g",jk,t3[0],t3[1],t3[2],t3[3]);*/
	    
	    sf_floatwrite(t3,4,Fl);
	}
    } else {
	for( jk=0; jk<nk; jk++) {
	    fseeko(tfile,jk*3*sizeof(float),SEEK_SET);
	    fread(    t2,     sizeof(float),3,tfile);
/*	    if(verb) sf_warning("%d, %g %g %g",jk,t2[0],t2[1],t2[2]);*/
	    
	    sf_floatwrite(t2,3,Fl);
	}
    }

    free(cube);
    unlink(tname);

    exit (0);
}

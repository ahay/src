/* 
 * Maps a cube to a list, given a threshold (clip value)
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

    axa   ax,ay,az; 
    axa   aa;       /* elements in the list */
    axa   aj;
    int   ix,iy,iz;
    int   nk,jk;
    
    float **cube;

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
    iaxa(Fc,&az,1); az.l="z";
    iaxa(Fc,&ax,2); ax.l="x";
    iaxa(Fc,&ay,3); ay.l="y";

    aa.n=0; aa.o=0; aa.d=1; aa.l=" ";
    if(ay.n>1) {
	if(verb) sf_warning("initiating 3D points");
	aj.n=4;aj.o=0;aj.d=1;aj.l=" ";
    } else {
	if(verb) sf_warning("initiating 2D points");
	aj.n=3;aj.o=0;aj.d=1;aj.l=" ";
    }
    /*------------------------------------------------------------*/

    cube = sf_floatalloc2(az.n,ax.n);

    tfile = sf_tempfile(&(tname), "w+b");

    for (iy=0;iy<ay.n;iy++) {
	if(verb) sf_warning("iy=%d",iy);
	
	sf_floatread(cube[0],az.n*ax.n,Fc);
	
	nk=0;
	for (ix=0; ix<ax.n; ix++) {
	    for (iz=0; iz<az.n; iz++) {
		if( abs(cube[ix][iz]) > clip) {
		    nk++;
		}
	    }
	}

	if(ay.n>1) {
	    jk=0;
	    for (ix=0; ix<ax.n; ix++) {
		for (iz=0; iz<az.n; iz++) {
		    if( abs(cube[ix][iz]) > clip) {
			t3[0] = ax.o + ix * ax.d;
			t3[1] = ay.o + iy * ay.d;
			t3[2] = az.o + iz * az.d;
			t3[3] = cube[ix][iz];
			
			fseeko(tfile,jk*4*SF_FLOAT,SEEK_SET);
			fwrite(   t3,   4*SF_FLOAT,1,tfile);
			jk++;
		    }
		}
	    }
	} else {
	    jk=0;
	    for (ix=0; ix<ax.n; ix++) {
		for (iz=0; iz<az.n; iz++) {
		    if( abs(cube[ix][iz]) > clip) {
			t2[0] = ax.o + ix * ax.d;
			t2[1] = az.o + iz * az.d;
			t2[2] = cube[ix][iz];
			
			fseeko(tfile,jk*3*SF_FLOAT,SEEK_SET);
			fwrite(   t2,   3*SF_FLOAT,1,tfile);
			
			jk++;
		    }
		}
	    }
	} /* else ay.n=1 */

	aa.n+=nk;
	if(verb) raxa(aa);
    } /* iy */
    
    /* output axes */
    oaxa(Fl,&aj,1); if(verb) raxa(aj);
    oaxa(Fl,&aa,2); if(verb) raxa(aa);

    if( ay.n>1) {
	for( jk=0; jk<nk; jk++) {
	    fseeko(tfile,jk*4*SF_FLOAT,SEEK_SET);
	    fread(    t3,   4*SF_FLOAT,1,tfile);
	    if(verb) sf_warning("%d, %g %g %g %g",jk,t3[0],t3[1],t3[2],t3[3]);
	    
	    sf_floatwrite(t3,4,Fl);
	}
    } else {
	for( jk=0; jk<nk; jk++) {
	    fseeko(tfile,jk*3*SF_FLOAT,SEEK_SET);
	    fread(    t2,   3*SF_FLOAT,1,tfile);
	    if(verb) sf_warning("%d, %g %g %g",jk,t2[0],t2[1],t2[2]);
	    
	    sf_floatwrite(t2,3,Fl);
	}
    }

    free(cube);
    unlink(tname);
    exit (0);
}

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
    float  x, y, z;
    int   ix,iy,iz;
    int   nk,ik;
    
    float **cube;
    pt2d  *list2, p2;
    pt3d  *list3, p3;

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
    oaxa(Fl,&aj,1); if(verb) raxa(aj);

    p2.x=0;         p2.z=0; p2.v=0;
    p3.x=0; p3.y=0; p3.z=0; p3.v=0;

    /*------------------------------------------------------------*/

    cube = sf_floatalloc2(az.n,ax.n);

/*
    for (iy=0;iy<ay.n;iy++) {
	y = ay.o + iy * ay.d;
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
	aa.n+=nk;

	if(ay.n>1) {
	    list3 = (pt3d*) sf_alloc(nk,sizeof(*list3)); 
	    for( ik=0;ik<nk;ik++) {
		list3[ik] = p3;
	    }
	    nk=0;
	    for (ix=0; ix<ax.n; ix++) {
		x = ax.o + ix * ax.d;
		for (iz=0; iz<az.n; iz++) {
		    z = az.o + iz * az.d;
		    if( abs(cube[ix][iz]) > clip) {
			list3[nk].x = x;
			list3[nk].y = y;
			list3[nk].z = z;
			list3[nk].v = cube[ix][iz];
			nk++;
		    }
		}
	    }
	    writept3d(Fl,list3,nk,4);
	    free(list3);
	} else {
	    list2 = (pt2d*) sf_alloc(nk,sizeof(*list2));
	    for( ik=0;ik<nk;ik++) {
		list2[ik] = p2;
	    }
	    nk=0;
	    for (ix=0; ix<ax.n; ix++) {
		x = ax.o + ix * ax.d;
		for (iz=0; iz<az.n; iz++) {
		    z = az.o + iz * az.d;
		    if( abs(cube[ix][iz]) > clip) {
			list2[nk].x = x;
			list2[nk].z = z;
			list2[nk].v = cube[ix][iz];
			nk++;
		    }
		}
	    }
	    writept2d(Fl,list2,nk,3);
	    free(list2);
       }	
    }
*/

    sf_floatread (cube[0],az.n*ax.n,Fc);

    /* output axes */
    aa.n=1;
    oaxa(Fl,&aa,2); if(verb) raxa(aa);

    sf_floatwrite(cube[0],az.n*ax.n,Fl);

    free(cube);
    exit (0);
}

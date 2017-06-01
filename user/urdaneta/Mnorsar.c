/* Traveltime and amplitude estimation using wavefront construction. */
/*
  Copyright (C) 1993 The Board of Trustees of Stanford University
  
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
#include <stdio.h>

#include <rsf.h>
#include <rsfplot.h>

/* Hector Urdaneta. 4/2/93 */

#include "norsar.h"
#include "caustics.h"
#include "gridding.h"
#include "interp.h"

#define XMIN    0.02
#define YMIN    .87
#define XMAX    13.66
#define YMAX    9.24
#define FAT     9
#define SKINNY  0
#define THIN    2

static void setgraphics (float ox, float oz, float lengthx, float lengthz);
static void draw_rays (struct heptagon *cube, int nr);
static void draw_wavefronts (struct heptagon *cube, int nr, float DSmax);

int main(int argc, char* argv[])
{
    int nang, N;
    float dx, dz;
    int nx, nz;
    float ox, oz;
    int gnx, gnz;
    float gdx, gdz, gox, goz;
    struct point length;
    int inter;
    float TETAMAX, alpha2;
    FILE *outfile;
    int ii;
    int rays, wfront, gap;
    int lomx, first;
    int nr, nrmax, nt;
    int prcube, pr, ste;
    int ns, nou;
    float DSmax, dt, T, freq;
    float *vel;
    float ds, os, goox;
    float xmin, xmax, zmin, zmax;
    float depth;
    struct point *pos;
    struct heptagon *cube;
    struct grid *out;
    sf_file inp, ampl, time,dirx,dirz;

    sf_init(argc,argv);
    inp = sf_input("in");

/* GET MODEL PARAMETERS	*/
    if (!sf_histint(inp,"n1",&nz)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histfloat(inp,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histfloat(inp,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(inp,"o1",&oz)) sf_error("No o1= in input");
    if (!sf_histfloat(inp,"o2",&ox)) sf_error("No o2= in input");

/* GET TRACING PARAMETERS */
    if(!sf_getint("nang",&nang))        nang = 10;    /* Number of take-off angles */
    if(!sf_getint("rays",&rays))        rays = 0;     /* If draw rays */
    if(!sf_getint("wfront",&wfront))    wfront = 0;   /* If draw wavefronts */ 
    if(!sf_getint("gap",&gap))        	gap = 1;      /* Draw wavefronts every gap intervals */
    if(!sf_getint("inter",&inter))      inter = 1;    /* If use linear interpolation */
    if(!sf_getfloat("DSmax",&DSmax))    DSmax = 5;    /* Maximum distance between contiguos points of a wavefront */
    if(!sf_getfloat("dt",&dt))          dt = 0.0005;  /* time step */
    if(!sf_getint("nt",&nt))            nt = 5;       /* Number of time steps between wavefronts */
    if(!sf_getint("nrmax",&nrmax))      nrmax = 2000; /* Maximum number of points that define a wavefront */
    if(!sf_getint("lomx",&lomx))        lomx = 1;     /* Use Lomax's waveray method */
    if(!sf_getint("first",&first))      first = 1;    /* Obtain first arrivals only */
    if(!sf_getint("nou",&nou))      	nou = 6;      

/* GET GRIDDING PARAMETERS */
    if(!sf_getint("gnx",&gnx))		gnx = nx;    /* Coordinates of output grid */
    if(!sf_getint("gnz",&gnz))		gnz = nz;
    if(!sf_getfloat("gdx",&gdx))	gdx = dx;
    if(!sf_getfloat("gdz",&gdz))	gdz = dz;
    if(!sf_getfloat("gox",&goox))	goox = ox;
    if(!sf_getfloat("goz",&goz))	goz = oz;

/* GET LOMAX SPECIFIC PARAMETERS */
    if(!sf_getint("N",&N))                N = 3;         /* Number of control points */
    if(!sf_getfloat("TETAMAX",&TETAMAX))  TETAMAX = 1.5; /* Truncation parameter */
    if(!sf_getfloat("alpha2",&alpha2))    alpha2 = 4.0;  /* Width of gaussian weighting function */
    if(!sf_getfloat("freq",&freq))	  freq = 100.;   /* Pseudo-frequency of waverays */

/* GET DEBUGGING INFO */
    if(!sf_getint("prcube",&prcube)) 	prcube=0;        /* For debugging porpouses */
    if(!sf_getint("pr",&pr)) 		pr=0;            /* For debugging porpouses */

/* GET SOURCE LOCATIONS */
    if(!sf_getint("ns",&ns) || ns==0)	  ns=1;          /* Number of source locations */
    if(!sf_getfloat("ds",&ds))		  ds=1.;         /* interval between sources */
    if(!sf_getfloat("os",&os))		  os=0.;         /* first source location */
    if(!sf_getfloat("depth",&depth))	  depth=dz;      /* Depth location of sources */

    pos = (struct point *) sf_alloc (ns,sizeof(struct point));
    for(ii=0;ii<ns;ii++) {
	pos[ii] = makepoint(ii*ds + os, depth);
    }

/* PREPARE OUTPUT */
    ampl = sf_output("ampl");
    sf_putint(ampl,"n1",gnz);
    sf_putint(ampl,"n2",gnx);
    sf_putint(ampl,"n3",ns);
    sf_putfloat(ampl,"d1",gdz);
    sf_putfloat(ampl,"d2",gdx);
    sf_putfloat(ampl,"d3",ds);
    sf_putfloat(ampl,"o1",goz);
    sf_putfloat(ampl,"o2",goox);
    sf_putfloat(ampl,"o3",os);

    time = sf_output("time");
    sf_putint(time,"n1",gnz);
    sf_putint(time,"n2",gnx);
    sf_putint(time,"n3",ns);
    sf_putfloat(time,"d1",gdz);
    sf_putfloat(time,"d2",gdx);
    sf_putfloat(time,"d3",ds);
    sf_putfloat(time,"o1",goz);
    sf_putfloat(time,"o2",goox);
    sf_putfloat(time,"o3",os);

    dirx = sf_output("dirx");
    sf_putint(dirx,"n1",gnz);
    sf_putint(dirx,"n2",gnx);
    sf_putint(dirx,"n3",ns);
    sf_putfloat(dirx,"d1",gdz);
    sf_putfloat(dirx,"d2",gdx);
    sf_putfloat(dirx,"d3",ds);
    sf_putfloat(dirx,"o1",goz);
    sf_putfloat(dirx,"o2",goox);
    sf_putfloat(dirx,"o3",os);

    dirz = sf_output("dirz");
    sf_putint(dirz,"n1",gnz);
    sf_putint(dirz,"n2",gnx);
    sf_putint(dirz,"n3",ns);
    sf_putfloat(dirz,"d1",gdz);
    sf_putfloat(dirz,"d2",gdx);
    sf_putfloat(dirz,"d3",ds);
    sf_putfloat(dirz,"o1",goz);
    sf_putfloat(dirz,"o2",goox);
    sf_putfloat(dirz,"o3",os);


/* READ VELOCITY MODEL */
    vel = sf_floatalloc(nx*nz+2);
    sf_floatread(vel,nx*nz,inp); 

/* ALLOCATE MEMORY FOR OUTPUT */
    out = (struct grid *) sf_alloc (1,sizeof(struct grid));

    out->time = sf_floatalloc (gnx*gnz);
    out->ampl = sf_floatalloc (gnx*gnz);
    out->dirx = sf_floatalloc (gnx*gnz);
    out->dirz = sf_floatalloc (gnx*gnz);
    out->flag = sf_intalloc (gnx*gnz);

    T = 1. / freq;

    length = makepoint((nx-1)*dx,(nz-1)*dz);

    cube = (struct heptagon *) sf_alloc (nrmax,sizeof(struct heptagon));

/* FOR DEBUGGING PORPOUSES, PRINT TO FILE */
    if(pr||prcube) {
	outfile = fopen("junk","w");
    } else {
	outfile = NULL;
    }

/*  SET DISPLAY IN ORDER TO SHOW RAYS ON SCREEN */
/*  NOTE: THIS PROGRAM USES DIRECT CALLS TO LIB_VPLOT
 *  TO DRAW THE RAYS AND THE WAVEFRONTS */
    if(rays || wfront) {
	setgraphics(ox, oz, length.x, length.z);
/*
	vp_color(BLUE);
	for(ii=0;ii<gnx;ii++)  {
	    vp_umove(ii*gdx+gox, goz); vp_udraw(ii*gdx+gox, (gnz-1)*gdz+goz); }
	for(ii=0;ii<gnz;ii++) {
	    vp_umove(gox, ii*gdz+goz); vp_udraw((gnx-1)*gdx+gox, ii*gdz+goz); } 
*/
    }

    norsar_init(gnx,gnz,
		TETAMAX,N,alpha2,inter,
		nx,nz,ox,oz,dx,dz,length);

/*  ALGORITHM: 
 *    For every source: */
    for(ii=0;ii<ns;ii++) {
	ste = 0;
	gox = goox + pos[ii].x;
	sf_warning("\nSource #%d\n", ii);

/*	1.- Construct the inital wavefront 			*/
	nr = nang;
    	initial (pos[ii], cube, vel, dt, nt, T, lomx, nr, out);

	gridding_init(gnx,gnz,
		      gdx,gdz,gox,goz,
		      outfile);

/*	run while the wavefront is not too small 		*/
	while (nr > 4) {
	    ste++;

/*	    2.- Propagate wavefront 				*/
	    wavefront (cube, nr, vel, dt, nt, T, lomx);
	    if(prcube || pr) {
		fprintf(outfile,"\n\nwavefront");
		printcube(cube, nr, outfile);
	    }

/*	    3.- Get rid of caustics				*/
            if(first) {
		if(ste%2==1) {			
		    caustics_up (cube, 0, nr);
		} else {
		    caustics_down (cube, nr-1, nr);
		}					
		if(prcube || pr) {	
		    fprintf(outfile,"\n\ncaustics");
		    printcube(cube, nr, outfile);
		}
	    }

/*          4.- Eliminate rays that cross boundaries, defined
		by xmin, xmax, zmin, zmax.
		Note that the computational grid is a little bigger 
		than the ouput grid. 				*/

	    xmin = gox-nou*gdx;	xmax = 2*pos[ii].x-gox+nou*gdx;
	    zmin = oz-nou*gdz;	zmax = length.z+oz+nou*gdz;
            mark_pts_outofbounds (cube, nr, xmin, xmax, zmin, zmax);
	    if(prcube) {
                fprintf(outfile, "\n\nboundaries");
                printcube(cube, nr, outfile);
            }

/*          5.- Rearrange cube                                  */
            makeup(cube, &nr);
	    if(nr<4) break;
	    if(prcube || pr) {	
                fprintf(outfile, "\n\nmakeup");
                printcube(cube, nr, outfile);
            }

/*          6.- Calculate amplitudes for new wavefront 		*/
	    amplitudes (cube, nr); 
	    if(prcube) {
                fprintf(outfile, "\n\namplitudes");
                printcube(cube, nr, outfile);
            }

/*	    7.- Draw rays 					*/
	    if(rays)
		draw_rays (cube, nr); 

/*          8.- Draw wavefront 					*/
	    if(wfront && (ste%gap==0)) 
		draw_wavefronts (cube, nr, DSmax); 

/*          9.- Parameter estimation at receivers 		*/
	    gridding (cube, nr, out, DSmax, (ste-1)*nt*dt, vel, first); /* pos[ii]); */

/*          10.- Interpolate new points of wavefront 		*/
            interpolation (cube, &nr, nrmax, DSmax); /* 0); */
	    if(prcube) {
                fprintf(outfile,"\n\ninterpolation");
                printcube(cube, nr, outfile);
            }

/*	    11.- Prepare to trace new wavefront 		*/
	    movwavf (cube, nr);
	    if(prcube) {
                fprintf(outfile,"\n\nmovwavf");
                printcube(cube, nr, outfile);
            }

	    if((wfront || rays) && (ste%gap==0)) vp_erase();
	}

/*	Finally interpolate amplitude and traveltime values to
        receivers that has not being covered.			*/

	TwoD_interp (out, gnx, gnz);

	sf_floatwrite(out->time, gnx * gnz, time); 
	sf_floatwrite(out->ampl, gnx * gnz, ampl);  
	sf_floatwrite(out->dirx, gnx * gnz, dirx);
	sf_floatwrite(out->dirz, gnx * gnz, dirz);
    }

    if(pr||prcube) fclose(outfile);

   exit(0);
}

static void setgraphics (float ox, float oz, float lengthx, float lengthz)
{
    float dash[1];
    float gap[1];

    dash[0] = 0.1;
    dash[0] = 0.1;

    vp_init();

    vp_erase ();
    vp_tjust (TH_CENTER, TV_TOP);
    vp_tfont (ROMANC, VP_STROKE, VP_NO_CHANGE);
    vp_purge ();

    vp_orig (XMIN + .165 * XMAX , YMIN + .875 * YMAX);
    vp_clip (-VP_MAX, -VP_MAX, VP_MAX, VP_MAX);
    vp_color (VP_RED);
    vp_fat (THIN);
    vp_setdash (dash, gap, 0);

    vp_uorig (ox , oz);
    vp_scale (0.75 * XMAX / (lengthx), -0.75 * YMAX / (lengthz));

    return;
}

/*----------------------------------------------------------------------------*/

static void draw_rays (struct heptagon *cube, int nr)
{
    int ii;
    float dp[1], gp[1];

    dp[0] = 0.01; gp[0] = 0.008;

    vp_setdash(dp, gp, 0);
    vp_fat(SKINNY);
    vp_color(VP_RED);
    vp_umove(cube[0].x0.x, cube[0].x0.z);
    vp_udraw(cube[0].x1.x, cube[0].x1.z);
    for(ii=1;ii<nr;ii++) {
	if(cube[ii].cf==END && cube[ii-1].cf==END) continue;
	vp_umove(cube[ii].x0.x, cube[ii].x0.z);
	vp_udraw(cube[ii].x1.x, cube[ii].x1.z);
    }
    return;
}

/*----------------------------------------------------------------------------*/

static void draw_wavefronts (struct heptagon *cube, int nr, float DSmax)
{
    int ii;
    float dp[1], gp[1];
    
    dp[0] = 0.01; gp[0] = 0.008;

    vp_fat(SKINNY);
    vp_color(VP_GREEN);        
    vp_setdash(dp, gp, 0);
    for(ii=0;ii<nr;ii++) {
	if(cube[ii].cf==END) continue;
        vp_umove(cube[ii].x1.x, cube[ii].x1.z);
        vp_udraw(cube[(ii+1+nr)%nr].x1.x, cube[(ii+1+nr)%nr].x1.z);
    }   
    return;
}


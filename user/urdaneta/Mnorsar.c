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

/* Hector Urdaneta. 4/2/93 */

#include "norsar.h"

int main(int argc, char* argv[])
{
    int ii;
    int rays, wfront, gap;
    int lomx, first;
    int nr, nrmax, nt;
    int prcube, pr, ste;
    int ns, n4, nou;
    char label1[25], label2[25];
    float DSmax, dt, T, freq;
    float *vel;
    float d4, o4, ds, os, goox;
    float xmin, xmax, zmin, zmax;
    float offset, depth;
    struct point *pos;
    struct heptagon *cube;
    struct grid *out;
    sf_file inp;

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

    pos = (struct point *) alloc (sizeof(struct point) * ns);
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

/* READ VELOCITY MODEL */
    vel = sf_floatalloc(nx*nz+2);
    sf_floatread(vel,nx*nz,inp); 

/* ALLOCATE MEMORY FOR OUTPUT */
    out = (struct grid *) sf_alloc (1,sizeof(struct grid));

    out->time = sf_floatalloc (gnx*gnz);
    out->ampl = sf_floatalloc (gnx*gnz);
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
	    gridding (cube, nr, out, DSmax, (ste-1)*nt*dt, vel, first, pos[ii]);

/*          10.- Interpolate new points of wavefront 		*/
            interpolation (cube, &nr, nrmax, DSmax, 0);
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
	TwoD_interp (out);

	srite("time", out->time, sizeof(float) * gnx * gnz); 
	srite("ampl", out->ampl, sizeof(float) * gnx * gnz);  
    }

   if(pr||prcube)
	fclose(outfile);

   exit(0);
}



/*----------------------------------------------------------------------------*/

/*
*
* void mark_pts_outofbounds (location of points on wavefront,
*       number of points on wavefront, boundary) 
*
* "mark_pts_outofbounds" raises a flag on any point of the wavefront 
* that goes out of bounds. The boundary is define by points x0, x1,
* z0, z1.
*
*/

void mark_pts_outofbounds (cube, nr, x0, x1, z0, z1)
struct heptagon *cube;
float x0, x1, z0, z1;
int nr;
{
    register ii;
    struct point pt;
	
    for(ii=0;ii<nr;ii++) {
	pt = cube[ii].x1;
	if(smaller(pt.x, x0) || bigger(pt.x, x1) || smaller(pt.z, z0)                     || bigger(pt.z, z1)) 
	    cube[ii].cf=OUT_OF_BOUNDS;
    }
    return;
}

/*----------------------------------------------------------------------------*/

/*
*
* void makeup (location of points on wavefront,
*       number of points on wavefront)
*
* "makeup" takes off, from the wavefront, any point that goes
* out of boundaries or that belongs to a caustic.
*
*/

void makeup (cube, nr)
struct heptagon *cube;
int *nr;
{
    register ii;
    struct heptagon temp;
    int sp;

    sp = 0;
    for(ii=0;ii<*nr;ii++) {
	temp = cube[ii];
	switch (temp.cf) {
	    case OUT_OF_BOUNDS:
		if(cube[(ii-1+*nr)%*nr].cf != OUT_OF_BOUNDS) 		{
		    if(cube[(ii-sp-2+*nr)%*nr].cf!=END) 
			cube[(ii-sp-1+*nr)%*nr].cf = END;
		    else sp++;						}
	    case CAUSTIC:
		sp++;
		break;
	    case IN_BOUNDS:
	    default:
		cube[(ii-sp+*nr)%*nr] = temp;
		break;
	}
    }
    *nr -= sp;
    return;
}

/*----------------------------------------------------------------------------*/

/*
*
* amplitudes (location of points on wavefront,
*       number of points on wavefront)
*
* This subroutine obtains the amplitudes of the new wavefront
* by calculating the geometrical spreading factor and 
* multiplying it with the old amplitude value.
*
* Notice that in order to calculate the amplitude for ending
* points of the wavefront; i.e., ii==0 or ii==nr-1, the code
* "wraps" to the other end of the wavefront by using moduli
* nr. 
*
*
*/

void amplitudes (cube, nr)
struct heptagon *cube;
int nr;
{
    register ii;
    int cf1;
    struct point A0, A1, A2;
    struct point B0, B1, B2;
    float R1, R2, r1, r2;

    for(ii=0;ii<nr;ii++) {
	A1 = cube[ii].x0;
	B1 = cube[ii].x1;

	if(cube[(ii-1+nr)%nr].cf==END) {
	    A0 = A1; B0 = B1;
	} else {
	    A0 = cube[(ii-1+nr)%nr].x0;
	    B0 = cube[(ii-1+nr)%nr].x1;
	}

	if(cube[(ii+nr)%nr].cf==END) {
	    A2 = A1; B2 = B1;
	} else {
	    A2 = cube[(ii+1+nr)%nr].x0;
	    B2 = cube[(ii+1+nr)%nr].x1;
	}

	r1 = dist(A0, A1); r2 = dist(A1, A2);
	R1 = dist(B0, B1); R2 = dist(B1, B2);
	if(R1+R2 > 0)
	    cube[ii].ampl *= sqrt((r1+r2) / (R1+R2));
    }
    return;
}

/*----------------------------------------------------------------------------*/

/*
* void interpolation (location of points on wavefront, number of points
*	on wavefront, maximum number of points per wavefront, maximum
* 	allowed distance between two contiguos rays).
*
* interpolation checks the distance between two contiguos points in
* the wavefront, if its bigger than DSmax, then it calculates how
* many rays does it has to interpolate in between.
*
* interpolation calls a pair of subroutines which lay in ./interp.c
* Both subroutines use a third order polynome interpolation.
*
*/

void interpolation (cube, nr, nrmax, DSmax)
struct heptagon *cube;
int *nr, nrmax;
float DSmax;
{
    register ii, jj;
    struct point pt0, pt1, pt2, pt3;
    float A0, A1, A2, A3, s;
    float an0, an1, an2, an3;
    int cf2, sp, nnc;
    struct heptagon temp;

    for(ii=0;ii<*nr;ii++) {

	if(cube[ii].cf == END)  continue;	
	nnc = ROUND(dist(cube[ii].x1, cube[(ii+1+*nr)%*nr].x1) / DSmax);
	if(nnc==0) continue;

        if(*nr + 1 >= nrmax)
            seperr("\nThe wavefront has grown bigger than the amount of \nmemory you alocated for it. %d\n", *nr);

	pt1 = cube[ii].x1;
	pt2 = cube[(ii+1+*nr)%*nr].x1;
	A1 = cube[ii].ampl;
        A2 = cube[(ii+1+*nr)%*nr].ampl;
	an1 = cube[ii].angle;
	an2 = cube[(ii+1+*nr)%*nr].angle;
	if(ii==*nr-1) an2 += 2.*PI;
	cf2 = cube[(ii+1+*nr)%*nr].cf;

	sp=ii-1;
	while(cube[(sp+*nr)%*nr].cf==NEWC) sp--;
	sp = (sp+*nr)%*nr;
	if(cube[sp].cf==END) {
	    pt0 = pt1;
	    A0 = A1;
	    an0 = an1;
	} else {
	    pt0 = cube[sp].x1;
	    A0 = cube[sp].ampl;
	    an0 = cube[sp].angle;
	    if(ii==0) an0 -= 2.*PI;
	}

	sp=ii+2;
	while(cube[(sp+*nr)%*nr].cf==NEWC) sp++;
	sp = (sp+*nr)%*nr;
	if(cube[sp].cf==END) {
	    pt3 = pt2;
	    A3 = A2;	
	    an3 = an2;
	} else {
	    pt3 = cube[sp].x1;
	    A3 = cube[sp].ampl;
	    an3 = cube[sp].angle;
	    if(ii>=*nr-2) an3 += 2.*PI;
	}

	for(jj=0;jj<nnc;jj++) {	
	    s = (float) (jj + 1) / (nnc + 1);
	    temp.x0 = cube[ii].x0;
	    temp.cf = NEWC;
/* Call interpolating subroutines to obtain ray angle, amplitude, and
   position on the wavefront, for the new points.		*/
	    temp.angle = realinterp (pt0, pt1, pt2, pt3, an0, an1, an2, an3, s);
	    temp.ampl = realinterp (pt0, pt1, pt2, pt3, A0, A1, A2, A3, s);
	    temp.x1 = ptinterp (pt0, pt1, pt2, pt3, s);
/* Push new wavefront point into cube				*/
	    push(temp, cube, nr, ++ii);
	}
    }
    return;
}

/*----------------------------------------------------------------------------*/

void draw_rays (cube, nr)
struct heptagon *cube;
int nr;
{
    register ii;
    float dp[1], gp[1];

    dp[0] = 0.01; gp[0] = 0.008;

    vp_setdash(dp, gp, 0);
    vp_fat(SKINNY);
    vp_color(RED);
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

void draw_wavefronts (cube, nr, DSmax)
struct heptagon *cube;
int nr;
float DSmax;
{
    register ii;
    float dp[1], gp[1];
    
    dp[0] = 0.01; gp[0] = 0.008;

    vp_fat(SKINNY);
    vp_color(GREEN);        
    vp_setdash(dp, gp, 0);
    for(ii=0;ii<nr;ii++) {
	if(cube[ii].cf==END) continue;
        vp_umove(cube[ii].x1.x, cube[ii].x1.z);
        vp_udraw(cube[(ii+1+nr)%nr].x1.x, cube[(ii+1+nr)%nr].x1.z);
    }   
    return;
}

/*----------------------------------------------------------------------------*/

void movwavf (cube, nr)
struct heptagon *cube;
int nr;
{
    register ii;

    for(ii=0;ii<nr;ii++) {
	cube[ii].x0 = cube[ii].x1;
	if(cube[ii].cf==NEWC) cube[ii].cf = IN_BOUNDS;
    }
    return;
}

/******************************** END *****************************************/

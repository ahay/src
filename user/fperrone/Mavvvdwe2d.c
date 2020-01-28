/* 2D acoustic variable-velocity variable-density time-domain FD modeling 

The code uses a standard second-order stencil in time.
The coefficients of the spatial stencil are computed 
by matching the transfer function of the discretized 
first-derivative operator to the ideal response. 
The optimized coefficients minimize dispersion 
given that particular size of the stencil.

The term 
	ro div (1/ro grad (u))
where
	ro   : density
	div  : divergence op
	grad : gradient  op
	u    : wavefield

is implemented in order to obtain a positive semi-definite operator.

The code implements both the forward (adj=n) and adjoint (adj=y) modeling operator.

============= FILE DESCRIPTIONS   ========================      

Fdat.rsf - An RSF file containing your data in the following format:
			axis 1 - Receiver location
			axis 2 - Time

Fwav.rsf - An RSF file containing your VOLUME DENSITY INJECTION RATE 
           wavelet information.  The sampling interval, origin time, 
           and number of time samples will be used as the defaults for the modeling code.
	       i.e. your wavelet needs to have the same length and parameters that you want to model with!
	       The first axis is the number of source locations.
	       The second axis is time.

Fvel.rsf - An N dimensional RSF file that contains the values for the velocity field at every point in the computational domain.

Fden.rsf - An N dimensional RSF file that contains the values for density at every point in the computational domain.

Frec.rsf - Coordinates of the receivers
		axis 1 - (x,z) of the receiver
		axis 2 - receiver index

Fsou.rsf - Coordinate of the sources
		axis 1 - (x,y,z) of the source
		axis 2 - source index

Fwfl.rsf - Output wavefield

verb=y/n - verbose flag

snap=y/n - wavefield snapshots flag

free=y/n - free surface on/off

dabc=y/n - absorbing boundary conditions on/off

jdata    - data sampling 

jsnap    - wavefield snapshots sampling


Author: Francesco Perrone
Date: November 2013
 */

/*
	Copyright (C) 2013 Colorado School of Mines

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
#ifdef _OPENMP
#include <omp.h>
#endif

#include <time.h>
/* check: dt<= 0.2 * min(dx,dz)/vmin */

#define NOP 3 /* derivative operator half-size */

/* LS coefficients */
#define C1 +1.1989919
#define C2 -0.08024696
#define C3 +0.00855954

enum adj_t{FWD,ADJ};


int main(int argc, char* argv[])
{
  // command line parameters
  // flags
  bool verb; // verbosity
  bool fsrf; // free surface
  bool dabc; // absorbing boundaries
  bool adj;  // adjoint flag
  // numerical arguments
  int  nb;

  /* I/O files */
  sf_file Fwav=NULL; /* wavelet   */
  sf_file Fsou=NULL; /* sources   */
  sf_file Frec=NULL; /* receivers */
  sf_file Fvel=NULL; /* velocity  */
  sf_file Fden=NULL; /* density   */
  sf_file Fdat=NULL; /* data      */
  sf_file Fwfl=NULL; /* wavefield */

  /* I/O arrays */
  float  *ww=NULL;  /* wavelet   */
  pt2d   *ss=NULL;  /* sources   */
  pt2d   *rr=NULL;  /* receivers */
  float  *dd=NULL;  /* data      */

  float *tt=NULL;
  float *ro=NULL;  /* density */
  float *uat=NULL; /* 1st derivative of wavefield */
  float *vp=NULL;  /* velocity */
  float *vt=NULL;  /* temporary vp*vp * dt*dt */

  /* for benchmarking */
  clock_t start_t, end_t;
  float total_t;

  /*------------------------------------------------------------*/
  /* init RSF */
  sf_init(argc,argv);

  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  // command line parameters
  if(! sf_getbool("verb",&verb)) verb=false; /* Verbosity flag */
  if(! sf_getbool("free",&fsrf)) fsrf=false; /* Free surface flag */
  if(! sf_getbool("dabc",&dabc)) dabc=false; /* Absorbing BC */
  if(! sf_getbool( "adj",&adj )) adj=false;  /* Adjoint flag*/
  /*------------------------------------------------------------*/

  /*------------------------------------------------------------*/
  /* I/O files */
  Fwav = sf_input ("in" );  /* wavelet   */
  Fvel = sf_input ("vel");  /* velocity  */
  Fden = sf_input ("den");  /* density   */
  Fsou = sf_input ("sou");  /* sources   */
  Frec = sf_input ("rec");  /* receivers */
  Fdat = sf_output("out");  /* data      */
  Fwfl = sf_output("wfl");  /* wavefield */

  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  /* absorbing boundary */
  if( !sf_getint("nb",&nb) || nb<NOP) nb=NOP;

  start_t=clock();
  switch (adj){
  case FWD:
    // fwdextrap2d();
    break;
  case ADJ:
    // adjextrap2d();
    break;
  }

  end_t = clock();
  if(verb) fprintf(stderr,"\n");

  if (verb){
    total_t = (float)(end_t - start_t) / CLOCKS_PER_SEC;
    fprintf(stderr,"Total time taken by CPU: %g\n", total_t  );
  }

  exit (0);
}

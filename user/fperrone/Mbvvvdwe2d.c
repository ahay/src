/* Born variable-density variable-velocity acoustic 2D time-domain FD modeling */
/*
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

The "reflectivity" that is used in the code is intended to be function of the 
change in VELOCITY. In particular, it is supposed to be equal to the product between 
the background and the perturbation in the velocity field, that is, the linear term in
the perturbation when you expand the square of the perturbed velocity
	
	v^2 = (vo + vp)^2 ~ vo^2 + 2*vo*vp
	
by assuming the perturbation is small compared to the background, the term vp^2 
can be neglected. The factor 2 is included in the source injection term in the code.

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

Fref.rsf - Reflectivity (same dimensions of the velocity model)

Frec.rsf - Coordinates of the receivers
		axis 1 - (x,z) of the receiver
		axis 2 - receiver index

Fsou.rsf - Coordinate of the sources
		axis 1 - (x,y,z) of the source
		axis 2 - source index

Fwfl.rsf - Output wavefield

Fliw.rsf - linearized scattered wavefield

Flid.rsf - linearized scattered data

verb=y/n - verbose flag

snap=y/n - wavefield snapshots flag

free=y/n - free surface on/off



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
#include "prep_utils.h"
#include "kernels.h"

int main(int argc, char* argv[])
{
  /* I/O files */
  sf_file Fwav=NULL; /* wavelet   */
  sf_file Fvpert=NULL; /* velocity perturbation */
  sf_file Frpert=NULL; /* density perturbation */
  sf_file Fsou=NULL; /* sources   */
  sf_file Frec=NULL; /* receivers */
  sf_file Fvel=NULL; /* velocity  */
  sf_file Fden=NULL; /* density   */
  sf_file Fdat=NULL; /* data      */
  sf_file Fwfl=NULL; /* wavefield */

  /* Other parameters */
  bool verb;
  bool fsrf;
  bool adj;
  int nb;

  /* for benchmarking */
  clock_t start_t, end_t;
  float total_t;

  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  /*                   RSF INITIALISATION                       */
  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  sf_init(argc,argv);

  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  /*                   PARSE THE PARAMETERS                     */
  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  if (!sf_getbool("verb",&verb)) verb=true;   /* verbosity flag*/
  if (!sf_getbool("free",&fsrf)) fsrf=false;  /* free surface */
  if (!sf_getbool("adj",&adj))    adj=false;  /* adjoint flag */

  if( !sf_getint("nb",&nb)) nb=NOP; /* thickness of the absorbing boundary: NOP is the width of the FD stencil*/
  if (nb<NOP) nb=NOP;

  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  /*                       OPEN FILES                           */
  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  Fwav = sf_input ("in" );  /* wavelet   */
  Fvel = sf_input ("vel");  /* velocity  */
  Fden = sf_input ("den");  /* density   */
  Fsou = sf_input ("sou");  /* sources   */
  Frec = sf_input ("rec");  /* receivers */
  Fdat = sf_output("out");  /* data      */

  /* -------------------------------------------------------------*/
  /* -------------------------------------------------------------*/
  /*                   CLOSE FILES AND EXIT                       */
  /* -------------------------------------------------------------*/
  /* -------------------------------------------------------------*/
  if (Fwav!=NULL) sf_fileclose(Fwav);
  if (Fvpert!=NULL) sf_fileclose(Fvpert);
  if (Frpert!=NULL) sf_fileclose(Frpert);
  if (Fsou!=NULL) sf_fileclose(Fsou);
  if (Frec!=NULL) sf_fileclose(Frec);
  if (Fvel!=NULL) sf_fileclose(Fvel);
  if (Fden!=NULL) sf_fileclose(Fden);
  if (Fdat!=NULL) sf_fileclose(Fdat);
  if (Fwfl!=NULL) sf_fileclose(Fwfl);

  exit (0);
}

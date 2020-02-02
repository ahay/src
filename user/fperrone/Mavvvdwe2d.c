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

#include "prep_utils.h"
#include "kernels.h"
#include <time.h>
/* check: dt<= 0.2 * min(dx,dz)/vmin */

static void print_param(in_para_struct_t in_para){
  sf_warning("verbosity          = %s",((in_para.verb==false)?"no":"yes"));
  sf_warning("free surface       = %s",((in_para.fsrf==false)?"no":"yes"));
  sf_warning("absorbing boundary = %s",((in_para.dabc==false)?"no":"yes"));
  if (in_para.dabc) sf_warning("- sponge thickness = %d",in_para.nb);
}

int main(int argc, char* argv[])
{
  // command line parameters
  in_para_struct_t in_para;

  /* I/O files */
  sf_file Fwav=NULL; /* wavelet   */
  sf_file Fsou=NULL; /* sources   */
  sf_file Frec=NULL; /* receivers */
  sf_file Fvel=NULL; /* velocity  */
  sf_file Fden=NULL; /* density   */
  sf_file Fdat=NULL; /* data      */
  sf_file Fwfl=NULL; /* wavefield */

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
  /*                COMMAND LINE PARAMETERS                     */
  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  if(! sf_getbool("verb",&(in_para.verb))) in_para.verb=false; /* Verbosity    */
  if(! sf_getbool("free",&(in_para.fsrf))) in_para.fsrf=false; /* Free surface */
  if(! sf_getbool("dabc",&(in_para.dabc))) in_para.dabc=false; /* Absorbing BC */
  if(! sf_getbool( "adj",&(in_para.adj) )) in_para.adj=false;  /* Adjointness  */

  if( !sf_getint("nb",&(in_para.nb)) || in_para.nb<NOP) in_para.nb=NOP;

  if (in_para.verb)
    print_param(in_para);

  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  /*                       OPEN FILES                           */
  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  Fwav = sf_input ("wav" ); /* wavelet   */
  Fvel = sf_input ("vel");  /* velocity  */
  Fden = sf_input ("den");  /* density   */
  Fsou = sf_input ("sou");  /* sources   */
  Frec = sf_input ("rec");  /* receivers */
  Fdat = sf_output("data"); /* data      */
  Fwfl = sf_output("wfl");  /* wavefield */

  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  /*                      READ AXES                             */
  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  sf_warning("WAVELET axes..");
  sf_axis axWav[2];
  axWav[0] = sf_iaxa(Fwav,1);
  sf_setlabel(axWav[0],"shot");
  sf_setunit(axWav[0],"1");
  if(in_para.verb) sf_raxa(axWav[0]); /* shot */

  axWav[1] = sf_iaxa(Fwav,2);
  sf_setlabel(axWav[1],"time");
  sf_setunit(axWav[1],"s");
  if(in_para.verb) sf_raxa(axWav[1]); /* time */

  sf_warning("VELOCITY model axes..");
  sf_axis axVel[2];
  axVel[0] = sf_iaxa(Fvel,1);
  sf_setlabel(axVel[0],"z");
  sf_setunit(axVel[0],"m");
  if(in_para.verb) sf_raxa(axVel[0]); /* depth */

  axVel[1] = sf_iaxa(Fvel,2);
  sf_setlabel(axVel[1],"x");
  sf_setunit(axVel[1],"m");
  if(in_para.verb) sf_raxa(axVel[1]); /* lateral */

  sf_warning("DENSITY model axes..");
  sf_axis axDen[2];
  axDen[0] = sf_iaxa(Fden,1);
  sf_setlabel(axDen[0],"z");
  sf_setunit(axDen[0],"m");
  if(in_para.verb) sf_raxa(axDen[0]); /* depth */

  axDen[1] = sf_iaxa(Fden,2);
  sf_setlabel(axDen[1],"x");
  sf_setunit(axDen[1],"m");
  if(in_para.verb) sf_raxa(axDen[1]); /* lateral */

  sf_warning("SHOT COORDINATES axes..");
  sf_axis axSou[2];
  axSou[0] = sf_iaxa(Fsou,1);
  sf_setlabel(axSou[0],"shot");
  sf_setunit(axSou[0],"1");
  if(in_para.verb) sf_raxa(axSou[0]); /* shot */

  axSou[1] = sf_iaxa(Fsou,2);
  sf_setlabel(axSou[1],"coords");
  sf_setunit(axSou[1],"1");
  if(in_para.verb) sf_raxa(axSou[1]); /* coords */

  sf_warning("RECEIVER COORDINATES axes..");
  sf_axis axRec[2];
  axRec[0] = sf_iaxa(Frec,1);
  sf_setlabel(axRec[0],"s");
  sf_setunit(axRec[0],"1");
  if(in_para.verb) sf_raxa(axRec[0]); /* shot */

  axRec[1] = sf_iaxa(Frec,2);
  sf_setlabel(axRec[1],"coords");
  sf_setunit(axRec[1],"1");
  if(in_para.verb) sf_raxa(axRec[1]); /* coords */

  sf_warning("CHECK MODEL DIMENSIONS..");
  if ((sf_n(axDen[0])!=sf_n(axVel[0])) ||
      (sf_n(axDen[1])!=sf_n(axVel[1]))){
    sf_error("Inconsistent model dimensions!");
    return -1;
  }
  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  /*                       PREPARE STRUCTURES                   */
  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  if (in_para.verb) sf_warning("Allocate structures..");
  wfl_struct_t *wfl = calloc(1,sizeof(wfl_struct_t));
  acq_struct_t *acq = calloc(1,sizeof(acq_struct_t));
  mod_struct_t *mod = calloc(1,sizeof(mod_struct_t));

  // PREPARE THE ACQUISITION STRUCTURE
  if (in_para.verb) sf_warning("Prepare the acquisition geometry structure..");
  prepare_acquisition_2d(acq, axSou, axRec, axWav, Fsou, Frec,Fwav);

  // PREPARE THE MODEL PARAMETERS CUBES
  if (in_para.verb) sf_warning("Read parameter cubes..");
  prepare_model_2d(mod,in_para,axVel,axDen,Fvel,Fden);

  // PREPARATION OF THE WAVEFIELD STRUCTURE
  if (in_para.verb) sf_warning("Prepare the wavefields for modeling..");
  prepare_wfl_2d(wfl,mod,Fdat,Fwfl,in_para);
  if (in_para.verb) sf_warning("Prepare the absorbing boundary..");
  setupABC(wfl);

  // WAVEFIELD HEADERS
  sf_oaxa(Fwfl,axVel[0],1);
  sf_oaxa(Fwfl,axVel[1],2);
  sf_oaxa(Fwfl,axWav[1],3);

  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  /*                  EXTRAPOLATION KERNEL                      */
  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  if (in_para.verb) sf_warning("Start Extrapolation..");
  start_t=clock();

  if (!in_para.adj)
    fwdextrap2d(wfl,acq,mod);
  else
    adjextrap2d(wfl,acq,mod);

  end_t = clock();

  if (in_para.verb){
    total_t = (float)(end_t - start_t) / CLOCKS_PER_SEC;
    sf_warning("Total time taken by CPU: %g", total_t );
  }

  /* -------------------------------------------------------------*/
  /* -------------------------------------------------------------*/
  /*                            FREE MEMORY                       */
  /* -------------------------------------------------------------*/
  /* -------------------------------------------------------------*/
  if (in_para.verb) sf_warning("Free memory..");
  clear_wfl_2d(wfl);
  free(wfl);

  clear_acq_2d(acq);
  free(acq);

  clear_model_2d(mod);
  free(mod);

  /* -------------------------------------------------------------*/
  /* -------------------------------------------------------------*/
  /*                   CLOSE FILES AND EXIT                       */
  /* -------------------------------------------------------------*/
  /* -------------------------------------------------------------*/
  if (in_para.verb) sf_warning("Close files..");
  if (Fwav!=NULL) sf_fileclose(Fwav);
  if (Fsou!=NULL) sf_fileclose(Fsou);
  if (Frec!=NULL) sf_fileclose(Frec);
  if (Fvel!=NULL) sf_fileclose(Fvel);
  if (Fden!=NULL) sf_fileclose(Fden);
  if (Fdat!=NULL) sf_fileclose(Fdat);
  if (Fwfl!=NULL) sf_fileclose(Fwfl);

  if (in_para.verb) sf_warning("ALL DONE!");
  exit (0);
}

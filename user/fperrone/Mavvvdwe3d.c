/* 3D acoustic variable-velocity variable-density time-domain FD modeling 

The code uses a standard second-order stencil in time.
The coefficients of the spatial stencil are computed
by matching the transfer function of the 6-point discretized
first-derivative operator to the ideal response.

The code implements the linearized operator obtained from the
system of first-order PDEs parametrized in incompressibility and density

dv/dt = - 1./rho * grad(p)
dp/dt = - K * div(v)

where
  rho  : density
  K    : incompressibility
  div  : divergence operator
  grad : gradient  operator
  p,v    : pressure and particle velocity wavefields

The models supplied by the user are wave speed and density, the code performs
the conversion internally to buoyancy (inverse density) and incompressibility.

Author: Francesco Perrone
Date: November 2020
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
#include <time.h>

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

  /* Other parameters */
  bool verb;
  bool fsrf;
  bool adj;
  bool dptf;
  bool snap;
  bool dabc;
  int nb;
  int jsnap;


  /* for benchmarking */
  clock_t start_fwd_t, end_fwd_t;
  float total_fwd_t=0;
  clock_t start_adj_t, end_adj_t;
  float total_adj_t=0;

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
  if(! sf_getbool("verb",&verb)) verb=false; /* Verbosity    */
  if(! sf_getbool("free",&fsrf)) fsrf=false; /* Free surface */
  if(! sf_getbool("dabc",&dabc)) dabc=false; /* Absorbing BC */
  if(! sf_getbool( "adj",&adj)) adj=false;  /* Adjointness  */
  if(! sf_getbool("snap",&snap)) snap=true; /* wavefield snapshots */

  if( !sf_getint("nb",&nb)) nb=NOP; /* thickness of the absorbing boundary: NOP is the width of the FD stencil */
  if (nb<NOP) nb=NOP;
  jsnap=1;
  if (snap){
    if (!sf_getint("jsnap",&jsnap)) jsnap=1;  /* undersampling factor for the wavefields*/
    if (jsnap<1) jsnap=1;
  }

  if (!sf_getbool( "dpt",&dptf)) dptf=false;  /* run dot product test */

  // fill the structure;
  in_para.verb=verb;
  in_para.fsrf=fsrf;
  in_para.dabc=dabc;
  in_para.adj=adj;
  in_para.snap=snap;
  in_para.nb=nb;
  in_para.jsnap=jsnap;
  in_para.dpt=dptf;

  if (in_para.verb)
    print_param(in_para);

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

  if (in_para.snap)
    Fwfl = sf_output("wfl");  /* wavefield */

  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  /*                      READ AXES                             */
  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  if(in_para.verb) sf_warning("WAVELET axes..");
  sf_axis axWav[2];
  axWav[0] = sf_iaxa(Fwav,1);
  sf_setlabel(axWav[0],"shot");
  sf_setunit(axWav[0],"1");
  if(in_para.verb) sf_raxa(axWav[0]); /* shot */

  axWav[1] = sf_iaxa(Fwav,2);
  sf_setlabel(axWav[1],"time");
  sf_setunit(axWav[1],"s");
  if(in_para.verb) sf_raxa(axWav[1]); /* time */

  if(in_para.verb) sf_warning("VELOCITY model axes..");
  sf_axis axVel[3];
  axVel[0] = sf_iaxa(Fvel,1);
  sf_setlabel(axVel[0],"z");
  sf_setunit(axVel[0],"m");
  if(in_para.verb) sf_raxa(axVel[0]); /* depth */

  axVel[1] = sf_iaxa(Fvel,2);
  sf_setlabel(axVel[1],"x");
  sf_setunit(axVel[1],"m");
  if(in_para.verb) sf_raxa(axVel[1]); /* lateral */

  axVel[2] = sf_iaxa(Fvel,3);
  sf_setlabel(axVel[2],"y");
  sf_setunit(axVel[2],"m");
  if(in_para.verb) sf_raxa(axVel[2]); /* lateral */

  if(in_para.verb) sf_warning("DENSITY model axes..");
  sf_axis axDen[3];
  axDen[0] = sf_iaxa(Fden,1);
  sf_setlabel(axDen[0],"z");
  sf_setunit(axDen[0],"m");
  if(in_para.verb) sf_raxa(axDen[0]); /* depth */

  axDen[1] = sf_iaxa(Fden,2);
  sf_setlabel(axDen[1],"x");
  sf_setunit(axDen[1],"m");
  if(in_para.verb) sf_raxa(axDen[1]); /* lateral */

  axDen[2] = sf_iaxa(Fden,3);
  sf_setlabel(axDen[2],"x");
  sf_setunit(axDen[2],"m");
  if(in_para.verb) sf_raxa(axDen[2]); /* lateral */

  if(in_para.verb) sf_warning("SHOT COORDINATES axes..");
  sf_axis axSou[2];
  axSou[0] = sf_iaxa(Fsou,1);
  sf_setlabel(axSou[0],"shot");
  sf_setunit(axSou[0],"1");
  if(in_para.verb) sf_raxa(axSou[0]); /* shot */

  axSou[1] = sf_iaxa(Fsou,2);
  sf_setlabel(axSou[1],"coords");
  sf_setunit(axSou[1],"1");
  if(in_para.verb) sf_raxa(axSou[1]); /* coords */

  if(in_para.verb) sf_warning("RECEIVER COORDINATES axes..");
  sf_axis axRec[2];
  axRec[0] = sf_iaxa(Frec,1);
  sf_setlabel(axRec[0],"s");
  sf_setunit(axRec[0],"1");
  if(in_para.verb) sf_raxa(axRec[0]); /* shot */

  axRec[1] = sf_iaxa(Frec,2);
  sf_setlabel(axRec[1],"coords");
  sf_setunit(axRec[1],"1");
  if(in_para.verb) sf_raxa(axRec[1]); /* coords */

  if(in_para.verb) sf_warning("CHECK MODEL DIMENSIONS..");
  if ((sf_n(axDen[0])!=sf_n(axVel[0])) ||
      (sf_n(axDen[1])!=sf_n(axVel[1])) ||
      (sf_n(axDen[2])!=sf_n(axVel[2])) ){
    sf_error("Inconsistent model dimensions!");

    exit(-1);
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

  // PREPARE THE MODEL PARAMETERS CUBES
  if (in_para.verb) sf_warning("Read parameter cubes..");
  prepare_model_3d(mod,in_para,axVel,axDen,Fvel,Fden);

  // PREPARE THE ACQUISITION STRUCTURE
  if (in_para.verb) sf_warning("Prepare the acquisition geometry structure..");
  prepare_acquisition_3d(acq, in_para, axSou, axRec, axWav, Fsou, Frec,Fwav);

  // PREPARATION OF THE WAVEFIELD STRUCTURE
  if (in_para.verb) sf_warning("Prepare the wavefields for modeling..");
  prepare_wfl_3d(wfl,mod,Fdat,Fwfl,in_para);


  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  /*                  EXTRAPOLATION KERNEL                      */
  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  if (in_para.verb) sf_warning("Start Extrapolation..");

  switch (in_para.adj){
  case FWD:
    start_fwd_t=clock();
//    fwdextrap3d(wfl,acq,mod);
    end_fwd_t=clock();
    total_fwd_t = (float)(end_fwd_t - start_fwd_t) / CLOCKS_PER_SEC;
    break;
  case ADJ:
    start_adj_t=clock();
//    adjextrap3d(wfl,acq,mod);
    end_adj_t=clock();
    total_adj_t = (float)(end_adj_t - start_adj_t) / CLOCKS_PER_SEC;
    break;
  }

  if (in_para.verb) sf_warning("Extrapolation completed..");
  /* -------------------------------------------------------------*/
  /* -------------------------------------------------------------*/
  /*                            FREE MEMORY                       */
  /* -------------------------------------------------------------*/
  /* -------------------------------------------------------------*/
  if (in_para.verb) sf_warning("Free memory..");
  clear_wfl_3d(wfl);
  free(wfl);

  clear_acq_3d(acq);
  free(acq);

  clear_model(mod);
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
  /*------------------------------------------------------------*/
  /*                    PROFILING INFO                          */
  /*------------------------------------------------------------*/
  sf_warning("=========================================================== ");
  sf_warning("PROFILING: [CPU time] ");
  sf_warning("=========================================================== ");
  if (in_para.adj==FWD){
    sf_warning("FWD operator                  :  %7.3g [s]", total_fwd_t );
  }
  else{
    sf_warning("ADJ operator                  : %7.3g [s]", total_adj_t );
  }
  sf_warning("=========================================================== ");
  sf_warning("=========================================================== ");

  exit (0);
}

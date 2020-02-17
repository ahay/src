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
#include<time.h>

int main(int argc, char* argv[])
{
  // command line parameters
  in_para_struct_t in_para;
  // running mode parameters
  born_setup_struct_t born_para;

  /* I/O files */
  // FWD and ADJ
  sf_file Fwav=NULL; /* wavelet   */
  sf_file Fsou=NULL; /* sources   */
  sf_file Frec=NULL; /* receivers */
  sf_file Fvel=NULL; /* velocity  */
  sf_file Fden=NULL; /* density   */

  // input for FWD, output for ADJ
  sf_file Fvpert=NULL; /* velocity perturbation */
  sf_file Frpert=NULL; /* density perturbation */

  // output FWD
  sf_file Fbdat=NULL; /* background data */
  // input ADJ
  sf_file Fsdat=NULL; /* scattered data */

  // auxiliary
  sf_file Fbwfl=NULL; /* background wavefield */
  sf_file Fswfl=NULL; /* scattered wavefield */

  /* Other parameters */
  bool verb;
  bool fsrf;
  bool adj;
  bool dabc;
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
  init_param(&in_para);
  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  /*                   PARSE THE PARAMETERS                     */
  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  if (!sf_getbool("verb",&verb)) verb=true;   /* verbosity flag*/
  if (!sf_getbool("free",&fsrf)) fsrf=false;  /* free surface */
  if (!sf_getbool("adj",&adj))    adj=false;  /* adjoint flag */
  if(! sf_getbool("dabc",&dabc)) dabc=false; /* Absorbing BC */

  if( !sf_getint("nb",&nb)) nb=NOP; /* thickness of the absorbing boundary: NOP is the width of the FD stencil*/
  if (nb<NOP) nb=NOP;

  // fill the structure;
  in_para.verb=verb;
  in_para.fsrf=fsrf;
  in_para.dabc=dabc;
  in_para.adj=adj;
  in_para.snap=true;  // we need the source wavefield
  in_para.nb=nb;

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

  bool vpert=false;
  bool rpert=false;

  if (in_para.adj){
    // ADJ BORN OPERATOR

    Fsdat = sf_input("pdata"); /* input pressure data to backproject (REQUIRED)*/

    born_para.outputVelPertImage=true;
    Fvpert= sf_output ("out");  /* default output: velocity perturbation image */

    born_para.outputDenPertImage=false;
    if (sf_getstring("dpert")){
      Frpert= sf_input ("rpert");  /* density perturbation */
      born_para.outputDenPertImage=true;
    }

    born_para.outputScatteredWfl=false;
    if (sf_getstring("swfl")){
      Fswfl = sf_output("swfl");  /* scattered wavefield*/
      born_para.outputScatteredWfl=true;
    }

  }
  else{
    // FWD BORN OPERATOR
    born_para.inputVelPerturbation=true;
    Fvpert= sf_input ("vpert");  /* velocity perturbation */
    vpert=true;

    born_para.inputDenPerturbation=false;
    if (sf_getstring("dpert")){
      Frpert= sf_input ("rpert");  /* density perturbation */
      born_para.inputDenPerturbation=true;
      rpert=true;
    }

    // these are aux output of the born forward modeling
    born_para.outputBackgroundWfl=false;
    if (sf_getstring("bwfl")){
      Fbwfl = sf_output("bwfl");  /* background wavefield*/
      born_para.outputBackgroundWfl=true;
    }

    born_para.outputBackgroundData=false;
    if (sf_getstring("bdat")){
      Fbdat = sf_output("bdat");  /* background data*/
      born_para.outputBackgroundData=true;
    }

    born_para.outputScatteredWfl=false;
    if (sf_getstring("swfl")){
      Fswfl = sf_output("swfl");  /* scattered wavefield*/
      born_para.outputScatteredWfl=true;
    }

    // this is the output of the born forward modeling
    Fsdat = sf_output("out");  /* scattered data*/
  }

  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  /*                      READ AXES                             */
  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  // from the input common to FWD and ADJ
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

  // ====================================================
  // CHECK
  // ====================================================
  sf_warning("CHECK MODEL DIMENSIONS..");
  if ((sf_n(axDen[0])!=sf_n(axVel[0])) ||
      (sf_n(axDen[1])!=sf_n(axVel[1]))){
    sf_error("Inconsistent model dimensions!");

    exit(-1);
  }
  // ====================================================
  // ====================================================

  sf_warning("WAVELET axes..");
  sf_axis axWav[2];
  axWav[0] = sf_iaxa(Fwav,1);
  axWav[1] = sf_iaxa(Fwav,2);
  sf_setlabel(axWav[0],"shot");
  sf_setlabel(axWav[1],"time");
  sf_setunit(axWav[0],"1");
  sf_setunit(axWav[1],"s");
  if(in_para.verb){
    sf_raxa(axWav[0]); /* shot */
    sf_raxa(axWav[1]); /* time */
  }

  sf_warning("SHOT COORDINATES axes..");
  sf_axis axSou[2];
  axSou[0] = sf_iaxa(Fsou,1);
  axSou[1] = sf_iaxa(Fsou,2);
  sf_setlabel(axSou[0],"coords");
  sf_setlabel(axSou[1],"shot");
  sf_setunit(axSou[0],"1");
  sf_setunit(axSou[1],"1");
  if(in_para.verb) {
    sf_raxa(axSou[0]); /* shot */
    sf_raxa(axSou[1]); /* coords */
  }

  sf_warning("RECEIVER COORDINATES axes..");
  sf_axis axRec[2];
  axRec[0] = sf_iaxa(Frec,1);
  axRec[1] = sf_iaxa(Frec,2);
  sf_setlabel(axRec[0],"coords");
  sf_setlabel(axRec[1],"receiver");
  sf_setunit(axRec[0],"1");
  sf_setunit(axRec[1],"1");
  if(in_para.verb) {
    sf_raxa(axRec[0]); /* coords */
    sf_raxa(axRec[1]); /* shots */
  }

  sf_axis axVelPert[2];
  sf_axis axDenPert[2];
  sf_axis axScData[2];
  if (!in_para.adj){
    // FWD BORN OPERATOR

    if (vpert){
      sf_warning("VELOCITY PERTURBATION model axes..");
      axVelPert[0] = sf_iaxa(Fvpert,1);
      axVelPert[1] = sf_iaxa(Fvpert,2);
      if(in_para.verb){
        sf_raxa(axVelPert[0]); /* depth */
        sf_raxa(axVelPert[1]); /* lateral */
      }
      sf_warning("CHECK MODEL DIMENSIONS..");
      if ((sf_n(axVel[0])!=sf_n(axVelPert[0]))||
          (sf_n(axVel[1])!=sf_n(axVelPert[1]))){
        sf_error("Inconsistent model dimensions!");

        exit(-1);
      }
    }

    if (rpert){
      sf_warning("DENSITY PERTURBATION model axes..");
      axDenPert[0] = sf_iaxa(Frpert,1);
      axDenPert[1] = sf_iaxa(Frpert,2);
      if(in_para.verb){
        sf_raxa(axDenPert[0]); /* depth */
        sf_raxa(axDenPert[1]); /* lateral */
      }

      sf_warning("CHECK MODEL DIMENSIONS..");
      if ((sf_n(axDen[0])!=sf_n(axDenPert[0]))||
          (sf_n(axDen[1])!=sf_n(axDenPert[1]))){
        sf_error("Inconsistent model dimensions!");

        exit(-1);
      }
    }

  }
  else{
    // ADJ BORN OPERATOR
    sf_warning("SCATTERED DATA axes..");
    axScData[0] = sf_iaxa(Fsdat,1);
    axScData[1] = sf_iaxa(Fsdat,2);

    if ((sf_n(axScData[0])!=sf_n(axRec[1]))){
      sf_error("Inconsistent receiver dimensions!");

      exit(-1);
    }

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
  prepare_model_2d(mod,in_para,axVel,axDen,Fvel,Fden);
  if (!in_para.adj)
    prepare_born_model_2d(mod,axVel,Fvpert,Frpert);


  // PREPARE THE ACQUISITION STRUCTURE
  if (in_para.verb) sf_warning("Prepare the acquisition geometry structure..");
  prepare_acquisition_2d(acq, in_para, axSou, axRec, axWav, Fsou, Frec,Fwav);
  if (in_para.adj)
    prepare_scatt_data_2d(acq,Fsdat);

  // PREPARATION OF THE WAVEFIELD STRUCTURE
  if (in_para.verb) sf_warning("Prepare the wavefields for modeling..");
  prepare_wfl_2d(wfl,mod,Fbdat,Fbwfl,in_para);
  prepare_born_wfl_2d(wfl,Fsdat,Fswfl);

  if (in_para.verb) sf_warning("Prepare the absorbing boundary..");
  setupABC(wfl);
  if (in_para.verb) sf_warning("Prepare the interpolation coefficients for source and receivers..");
  set_sr_interpolation_coeffs(acq,wfl);

  // WAVEFIELD HEADERS
  sf_axis axTimeWfl = sf_maxa(acq->ntsnap,
                              acq->ot,
                              acq->dt*in_para.jsnap);
  sf_setlabel(axTimeWfl,"time");
  sf_setunit(axTimeWfl,"s");

  sf_oaxa(Fbwfl,axVel[0],1);
  sf_oaxa(Fbwfl,axVel[1],2);
  sf_oaxa(Fbwfl,axTimeWfl,3);

  sf_oaxa(Fswfl,axVel[0],1);
  sf_oaxa(Fswfl,axVel[1],2);
  sf_oaxa(Fswfl,axTimeWfl,3);

  sf_oaxa(Fbdat,axRec[1],1);
  sf_axis axTimeData = sf_maxa( acq->ntdat,
                                acq->ot,
                                acq->dt);
  sf_oaxa(Fbdat,axTimeData,2);

  // HEADERS
  if (!in_para.adj){
    // DATA HEADERS
    sf_oaxa(Fsdat,axRec[1],1);
    sf_oaxa(Fsdat,axTimeData,2);
  }
  else{
    if (Fvpert){
      sf_axis axVpImg[2];
      axVpImg[0] = sf_maxa(sf_n(axVel[0]),
                           sf_o(axVel[0]),
                           sf_d(axVel[0]));
      axVpImg[1] = sf_maxa(sf_n(axVel[1]),
                           sf_o(axVel[1]),
                           sf_d(axVel[1]));

      sf_oaxa(Fvpert,axVpImg[0],1);
      sf_oaxa(Fvpert,axVpImg[1],2);
    }
    if (Frpert){
      sf_axis axRhImg[2];
      axRhImg[0] = sf_maxa(sf_n(axVel[0]),
                           sf_o(axVel[0]),
                           sf_d(axVel[0]));
      axRhImg[1] = sf_maxa(sf_n(axVel[1]),
                           sf_o(axVel[1]),
                           sf_d(axVel[1]));

      sf_oaxa(Frpert,axRhImg[0],1);
      sf_oaxa(Frpert,axRhImg[1],2);
    }


  }

  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  /*                  BACKGROUND WAVEFIELD                      */
  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  if (in_para.verb) sf_warning("Background wavefield extrapolation..");
  start_t=clock();
  fwdextrap2d(wfl,acq,mod);
  end_t = clock();

  if (in_para.verb){
    total_t = (float)(end_t - start_t) / CLOCKS_PER_SEC;
    sf_warning("Total time taken by CPU: %g", total_t );
  }

  //rewind the source wavefield
  sf_seek(wfl->Fwfl,0,SEEK_SET);

  // reset the wavefields
  reset_wfl(wfl);

  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  /*                  BORN OPERATOR                             */
  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  if (!in_para.adj){
    start_t=clock();
    // FWD BORN MODELING: model pert -> wfl
    if (in_para.verb) sf_warning("FWD Born operator..");

    // prepare the born sources
    make_born_sources_2d(wfl,mod,acq);

    // extrapolate secondary sources
    bornfwdextrap2d(wfl,acq,mod);
    end_t = clock();
  }
  else{
    start_t=clock();
    // ADJ BORN MODELING: wfl -> model pert
    if (in_para.verb) sf_warning("Adjoint Born operator..");

    // extrapolate data
    bornadjextrap2d(wfl,acq,mod);

    //rewind the scattered wavefield
    sf_seek(wfl->Fswfl,0,SEEK_SET);

    // prepare the born sources
    make_born_sources_2d(wfl,mod,acq);

    // stack wavefields
    stack_wfl_2d(Fvpert,Frpert,wfl,mod,acq);

    end_t = clock();
  }
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
  clear_born_wfl_2d(wfl);
  free(wfl);

  clear_acq_2d(acq);
  free(acq);

  clear_model_2d(mod);
  clear_born_model_2d(mod);
  free(mod);

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
  if (Fbdat!=NULL) sf_fileclose(Fbdat);
  if (Fsdat!=NULL) sf_fileclose(Fsdat);
  if (Fbwfl!=NULL) sf_fileclose(Fbwfl);
  if (Fswfl!=NULL) sf_fileclose(Fswfl);

  exit (0);
}

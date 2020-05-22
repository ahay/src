/* 1-D acoustic FWI with adjoint state method*/
/*
  Copyright (C) 2019 
  coded by Hongyu Zhou @ China University of Petroleum 
  (Email:carrot_cup@163.com)
  
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
#include"acoustic1D.h"
#include"adjoint_gradient.h"


int main(int argc, char *argv[])
{


    sf_init(argc,argv);bool verb;
    if(!sf_getbool("verb",&verb)) verb=false; /* verbosity */

    /* I/O files */
    sf_file Frec=sf_input("in");/*seismic record*/
    sf_file Fvel=sf_input("vel"); /*velocity model*/
    sf_file Fsrc = sf_input ("src"); /*source info*/
  
      
    int sx,rx,niter;
    if (!sf_getint("sx",&sx))sf_error("need sx");/*source position*/ 
    if (!sf_getint("rx",&rx)) sf_error("need rx"); /*reciever position*/
    if (!sf_getint("niter",&niter)) sf_error("need niter"); /*iteration number*/
    //sf_warning("/******************************************/");
    //sf_warning("sx=%d rx=%d,niter=%d",sx,rx,niter);

    //* cube axes */
    sf_axis at=sf_iaxa(Fsrc,1);int nt=sf_n(at);float dt=sf_d(at);
    sf_axis ax=sf_iaxa(Fvel,1);int nx=sf_n(ax);float dx=sf_d(ax);
   // sf_warning("nx=%d, dx=%g, nt=%d, dt=%g",nx,dx,nt,dt);
   // sf_warning("/******************************************/");

    float *vel_ini=sf_floatalloc(nx);sf_floatread(vel_ini,nx,Fvel);
    float *src=sf_floatalloc(nt);sf_floatread(src,nt,Fsrc);
    float *rec=sf_floatalloc(nt);sf_floatread(rec,nt,Frec);
    float *vel_ini_pert=sf_floatalloc(nx);
    float alpha=1.0;

     for (int iter = 0; iter < niter; ++iter)
    {

      // step 1. calc gradient with adjoint state method
      sf_warning("iter=%d   ",iter);
      float* gradient=adjoint_gradient(vel_ini,src,dx,nx,dt,nt,sx,rx,rec,1);
     //sf_error("TEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEST");
     // sf_warning("gradient[100]=%f",gradient[100]);
    //sf_file Fgradient = sf_output("gradient");/* seismic record @ rx*/
    //sf_oaxa(Fgradient,ax,1);
    //sf_floatwrite(gradient,nx,Fgradient);
      ///////////step 2. calc step length

      for (int ix = 0; ix< nx; ++ix) vel_ini_pert[ix]=vel_ini[ix]+alpha*gradient[ix];
      float* gradient_pert=adjoint_gradient(vel_ini_pert,src,dx,nx,dt,nt,sx,rx,rec,0);
     // sf_warning("gradient_pert[100]=%f",gradient_pert[100]);

      float num_sum=0.0,deno_sum=0.0;
      for (int ix = 0; ix < nx; ++ix)
      {
        num_sum=num_sum+gradient[ix]*gradient[ix];
        deno_sum=deno_sum+(gradient_pert[ix]-gradient[ix])*gradient[ix];
      } 
      alpha=-alpha*num_sum/deno_sum;//step length
      ///end of step 2
      for (int ix = 0; ix < nx; ++ix) vel_ini[ix]=vel_ini[ix]+alpha*gradient[ix]; //updating the velocity
      

      /* code */
    }


  
    sf_file Fvel_recover = sf_output("out");/* seismic record @ rx*/
    sf_oaxa(Fvel_recover,ax,1);
    sf_floatwrite(vel_ini,nx,Fvel_recover);
    sf_warning("ntttttttttttttttttttttttt_at_last=%d",nt);

    sf_close();
    exit(0);



}
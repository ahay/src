/* 1-D acoustic FWI with Perturbation method*/
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
#include <time.h>
#include"acoustic1D.h"

//Flow('vel_recover','rec','acoustic1D_FWI_ptb  vel=vel_ini.rsf src=src.rsf sx=10 rx=10' )


 int main(int argc, char  *argv[])
 {
    sf_init(argc,argv);bool verb;
    if(!sf_getbool("verb",&verb)) verb=false; /* verbosity */

    /* I/O files */
    sf_file Frec=sf_input("in");/*seismic record*/
    sf_file Fvel=sf_input("vel"); /*velocity model*/
    sf_file Fsrc = sf_input ("src"); /*source info*/
  
      
    int sx,rx,niter;float damp;
    if (!sf_getint("sx",&sx))sf_error("need sx");/*source position*/ 
    if (!sf_getint("rx",&rx)) sf_error("need rx"); /*reciever position*/
    if (!sf_getint("niter",&niter)) sf_error("need niter"); /*iteration number*/
    if (!sf_getfloat("damp",&damp)) sf_error("need damp"); /*damping factor*/
    sf_warning("/******************************************/");
    sf_warning("niter=%d,damp=%g, sx=%d rx=%d",niter,damp,sx,rx);
 
 
    //* cube axes */
    sf_axis at=sf_iaxa(Fsrc,1);int nt=sf_n(at);float dt=sf_d(at);
    sf_axis ax=sf_iaxa(Fvel,1);int nx=sf_n(ax);float dx=sf_d(ax);
    sf_warning("nx=%d, dx=%g, nt=%d, dt=%g",nx,dx,nt,dt);
    sf_warning("/******************************************/");


   



    float *vel_ini=sf_floatalloc(nx);sf_floatread(vel_ini,nx,Fvel);
    float *src=sf_floatalloc(nt);sf_floatread(src,nt,Fsrc);
    float *rec=sf_floatalloc(nt);sf_floatread(rec,nt,Frec);

	  float *rec_ini=sf_floatalloc(nt);float **rec_all_ini=sf_floatalloc2(nt,nx);
    float *vel_pert=sf_floatalloc(nx);float *rec_pert=sf_floatalloc(nt);
    float **Jac=sf_floatalloc2(nt,nx);float *gradient=sf_floatalloc(nx);float **Hessian_app=sf_floatalloc2(nx,nx);

    float pert=1.0; //float damping_factor=1e-9;

    for (int iter = 0; iter < niter; ++iter)
    {
          acoustic1D(vel_ini,src,dx,nx,dt,nt,sx,rx,rec_ini,rec_all_ini);

          for (int ix = 0; ix < nx; ++ix)
          {
            for (int ix = 0; ix < nx; ++ix) vel_pert[ix]=vel_ini[ix]; 
            vel_pert[ix]=vel_pert[ix]+pert;//adding pert
            acoustic1D(vel_pert,src,dx,nx,dt,nt,sx,rx,rec_pert,rec_all_ini);// pert response=>rec_pert
            for (int it = 0; it < nt; ++it) Jac[ix][it]=(rec_pert[it]-rec_ini[it])/pert;//obtain jacobian Matrix
          }

          for (int ix = 0; ix < nx; ++ix)
          {
          	gradient[ix]=0.0;
          	for (int it = 0; it < nt; ++it) gradient[ix]=gradient[ix]+Jac[ix][it]*(rec[it]-rec_ini[it]);// calc gradient
          } 
            
            

          for (int ix1 = 0; ix1 < nx; ++ix1)
            for (int ix2 = 0; ix2 < nx; ++ix2)
              {
              	 Hessian_app[ix1][ix2]=0.0;
              	 for (int it = 0; it < nt; ++it) Hessian_app[ix1][ix2]=Hessian_app[ix1][ix2]+Jac[ix1][it]*Jac[ix2][it];//calc Hessian app
                 if (ix1==ix2) Hessian_app[ix1][ix2]=Hessian_app[ix1][ix2]+damp;
              }

		  int NRHS=1,INFO;
		  int  *IPIV=sf_intalloc(nx);
		  sgesv_(&nx, &NRHS, Hessian_app[0], &nx, IPIV, gradient, &nx, &INFO);//gradient=Hessian_app\gradient
 

        for (int ix = 0; ix < nx; ++ix) vel_ini[ix]=vel_ini[ix]+gradient[ix]; //updating the velocity
        float error_sum=0;for (int it = 0; it < nt; ++it)  error_sum=error_sum+(rec_ini[it]-rec[it])*(rec_ini[it]-rec[it]); 
        sf_warning("iter=%d   ,error_sum=%f\n",iter,error_sum);

      /* code */
    }


    sf_file Fvel_recover = sf_output("out");/* seismic record @ rx*/
    sf_oaxa(Fvel_recover,ax,1);
    sf_floatwrite(vel_ini,nx,Fvel_recover);
    sf_warning("ntttttttttttttttttttttttt_at_last=%d",nt);

    //free(record);


  
    sf_close();
    exit(0);
//printf("hahahahha\n");



   return 0;
 }

/* 1-D acoustic wave propagation with CE absorbing boundary*/
/*
  Copyright (C) 2019 
  coded by carrot @ China University of Petroleum 
  (Email:1017014258@qq.com)
  
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

void acoustic1D(float *vel,float *src,float dx,int nx,float dt, int nt, int sx, int rx,float *record, float **record_all)
/*<1-D acoustic wave propagation with CE absorbing boundary>*/
{


      
    /*source position*/ 
    /*reciever position*/
    //sf_warning("sx=%d rx=%d",sx,rx);
   //sf_warning("nx=%d, dx=%g \n nt=%d, dt=%g",nx,dx,nt,dt);
//
//sf_warning("dt=%g",dt);
//printf("dt=%f\n",dt );
   //sf_error("test");

    //sf_error("test");

    float *prev=sf_floatalloc(nx);
    float *curr=sf_floatalloc(nx);
    float *next=sf_floatalloc(nx);

    //float **record_all=sf_floatalloc2(nt,nx);
    for (int ix = 0; ix < nx; ++ix)
    {
      prev[ix]=0.0;
      curr[ix]=0.0;
      next[ix]=0.0;
      /* code */
    }

    //float *record=sf_floatalloc(nt);
    for (int it = 0; it < nt; ++it)
    {
      record[it]=0.0;
    }


    for (int it = 0; it < nt; ++it)
    {
      for (int ix = 1; ix < nx-1; ++ix)
      {
        next[ix]=2*curr[ix]-prev[ix]+vel[ix]*vel[ix]*dt*dt*(curr[ix+1]+curr[ix-1]-2*curr[ix])/dx/dx;
      }
      curr[sx]= curr[sx]+src[it];//source term
      /* code */
      //CE boundary
      float rleft=vel[0]*dt/dx;    next[0]=(rleft*(next[1]-(prev[1]-prev[0]))-(prev[1]-2*curr[1]+next[1]+prev[0]-2*curr[0]))/(1+rleft);
      float rright=vel[nx-1]*dt/dx;next[nx-1]=(rright*(next[nx-2]+(prev[nx-1]-prev[nx-2]))-(prev[nx-2]-2*curr[nx-2]+next[nx-2]+prev[nx-1]-2*curr[nx-1]))/(1+rright);

     for (int ix = 0; ix < nx; ++ix)
      {
        prev[ix]=curr[ix];
        curr[ix]=next[ix];
        record_all[ix][it]=next[ix];
      }


      record[it]=next[rx];
      /*
      if (it%1000==0)
      {
        sf_warning("it=%d\n",it);
      
      }
      */
          

    }




    //sf_warning("nt=%d",nt);




}
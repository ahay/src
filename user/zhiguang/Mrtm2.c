/* 2-D exploding-reflector reverse-time migration/forward time modeling(exactly adjoint)
8th order in space, 2th order in time.
The adjoint of boundary condition was discussed with Pengliang Yang.
This code is written to implement least-squares exploding-reflector RTM, extra credit part in the sixth homework of "Seismic Imaging" taught by Sergey. */
/*
  Copyright (C) 2014 University of Texas at Austin
  
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

static int n1,n2,nb;
static float c[2][5],*abc,**vv;

void laplacian(float **u0,float **u1,float **u2,float **ud)
{ 
    int i1,i2,j;
    memset(ud[0],0,n1*n2*sizeof(float)); 

#ifdef _OPENMP
#pragma omp parallel for    \
   private(i1,i2,j)         \
   shared (u0,u1,u2,ud,n1,n2,vv,c)
#endif
    for(i2=4;i2<n2-4;i2++) 
	for(i1=4;i1<n1-4;i1++){
            for (j=0;j<=4;j++)
		ud[i2][i1]+=c[0][j]*(u1[i2][i1-j]+u1[i2][i1+j]);
		
	    for(j=0;j<=4;j++)
		ud[i2][i1]+=c[1][j]*(u1[i2-j][i1]+u1[i2+j][i1]);
		
	    u2[i2][i1]=2*u1[i2][i1]-u0[i2][i1]+ud[i2][i1]*vv[i2][i1];
		
	}
}

void adjoint(float **u1,float **ud)
{ 
    int i1,i2,j;
    memset(ud[0],0,n1*n2*sizeof(float));

#ifdef _OPENMP
#pragma omp parallel for    \
    private(i1,i2,j)        \
    shared(u1,ud,n1,n2,vv,c)
#endif
    for(i2=4;i2<n2-4;i2++) 
        for(i1=4;i1<n1-4;i1++){
	    for(j=0;j<=4;j++)
	        ud[i2][i1] += c[0][j]*(u1[i2][i1-j]*vv[i2][i1-j]+u1[i2][i1+j]*vv[i2][i1+j]);
	     
	    for(j=0;j<=4;j++)
	        ud[i2][i1] += c[1][j]*(u1[i2-j][i1]*vv[i2-j][i1]+u1[i2+j][i1]*vv[i2+j][i1]);
        }
}

void boundary(float **u0,float **u1)
{
    int i1,i2;

#ifdef _OPENMP
#pragma omp parallel for    \
    private(i1,i2)          \
    shared(u0,u1,abc,n1,n2,nb)
#endif
    for(i1=0;i1<n1;i1++){
	for(i2=0;i2<nb;i2++){			
	    u0[i2][i1]=abc[i2]*u0[i2][i1];
	    u1[i2][i1]=abc[i2]*u1[i2][i1];
	}	
	for(i2=n2-nb;i2<n2;i2++){			
	    u0[i2][i1]=abc[n2-i2-1]*u0[i2][i1];
	    u1[i2][i1]=abc[n2-i2-1]*u1[i2][i1];
	}	
    }
    
#ifdef _OPENMP
#pragma omp parallel for    \
    private(i1,i2)          \
    shared(u0,u1,abc,n1,n2,nb)
#endif
    for(i2=0;i2<n2;i2++){
	for(i1=0;i1<nb;i1++){			
	    u0[i2][i1]=abc[i1]*u0[i2][i1];
	    u1[i2][i1]=abc[i1]*u1[i2][i1];
	}
	for(i1=n1-nb;i1<n1;i1++){			
	    u0[i2][i1]=abc[n1-i1-1]*u0[i2][i1];
	    u1[i2][i1]=abc[n1-i1-1]*u1[i2][i1];
	}
    }
}

int main(int argc, char* argv[])
{
    bool adj; /* adjoint flag */
    int i1,i2,it,ib; /* index variables */
    int nt,nx,n0,n12;
    float dt,dx,dz,z0,t0,dt2,dz2,dx2,tmp;
    
    float **dd;
    float **u0,**u1,**u2,**ud; /* temporary arrays */      

    sf_file in,out,vel; /* I/O files */

    /* initialize Madagascar */
    sf_init(argc,argv);
    
    /* initialize OpenMP support */
    omp_init();

    if(!sf_getbool("adj",&adj)) adj=true;
    /* adjoint flag, 0: modeling, 1: migration */
    if(!sf_getint("nb",&nb)) nb=60;
    /* width of padding boundary */
    if(!sf_getint("n0",&n0)) n0=nb;
    /* surface */
	
    /* setup I/O files */
    in=sf_input("in");
    out=sf_output("out");
    vel = sf_input("velocity");
    /* velocity model */
    
    /* Dimensions */
    if(!sf_histint(vel,"n1",&n1)) sf_error("No n1 in velocity model");
    if(!sf_histint(vel,"n2",&n2)) sf_error("No n2 in velocity model");
    if(!sf_histfloat(vel,"d1",&dz)) sf_error("No d1 in velocity model");
    if(!sf_histfloat(vel,"d2",&dx)) sf_error("No d2 in velocity model");
    if(!sf_histfloat(vel,"o1",&z0)) sf_error("No o1 in velocity model");

    if(adj){/* migration */
        if(!sf_histint(in,"n1",&nt)) sf_error("No n1 in data file");
	if(!sf_histfloat(in,"d1",&dt)) sf_error("No d1 in data file");
	if(!sf_histfloat(in,"o1",&t0)) sf_error("No o1 in data file");
        if (!sf_histint(in,"n2",&nx) || nx!=n2)
        sf_error("Need n2=%d in data",n2);	    
	
	sf_putint(out,"n1",n1);
	sf_putfloat(out,"d1",dz);
	sf_putfloat(out,"o1",z0);
	sf_putstring(out,"label1","Depth");
	sf_putstring(out,"label2","Distance");
    }else{/* modeling */
	if(!sf_getint("nt",&nt)) sf_error("Need to input nt=");
	/* number of time steps */
	if(!sf_getfloat("dt",&dt)) sf_error("Need to input dt=");
	/* time sampling interval */
	if (!sf_getfloat("t0",&t0)) sf_error("Need to input t0=");
	/* time origin */
	
	sf_putint(out,"n1",nt);
	sf_putfloat(out,"d1",dt);
	sf_putfloat(out,"o1",t0);
	sf_putstring(out,"label1","Time");
	sf_putstring(out,"label2","Distance");
    }

    /* set Laplacian coefficients */
    dz2=1.0/(dz*dz);    
    c[0][1]=8.*dz2/5.;
    c[0][2]=-dz2/5.;
    c[0][3]=8.*dz2/315.;
    c[0][4]=-dz2/560;
    c[0][0]=-(c[0][1]+c[0][2]+c[0][3]+c[0][4]);
    
    dx2=1.0/(dx*dx);
    c[1][1]=8.*dx2/5.;
    c[1][2]=-dx2/5.;
    c[1][3]=8.*dx2/315.;
    c[1][4]=-dx2/560;
    c[1][0]=-(c[1][1]+c[1][2]+c[1][3]+c[1][4]);

    /* allocate array */
    abc=sf_floatalloc(nb);
	
    /* calculate boundary attenuation coefficients */
    for(ib=0;ib<nb;ib++){
	tmp=(nb-ib)*0.01;
	abc[ib]=expf(-tmp*tmp);
    }
        
    /* allocate arrays */
    vv=sf_floatalloc2(n1,n2);
    dd=sf_floatalloc2(nt,n2);
    
    n12=n1*n2;
    /* read velocity */
    sf_floatread(vv[0],n12,vel);
    
    /* read data or adjnull */
    if(adj) sf_floatread(dd[0],nt*n2,in);
    else memset(dd[0],0,4*nt*n2);
    
    /* allcoate temporary arrays */
    u0=sf_floatalloc2(n1,n2);
    u1=sf_floatalloc2(n1,n2);
    u2=sf_floatalloc2(n1,n2);
    ud=sf_floatalloc2(n1,n2);
    
    dt2=dt*dt;
    for(i2=0;i2<n2;i2++){
        for(i1=0;i1<n1;i1++){
            u0[i2][i1]=0.0;
            u1[i2][i1]=0.0;
            u2[i2][i1]=0.0;
            ud[i2][i1]=0.0;
            vv[i2][i1]*=vv[i2][i1]*dt2;
        }
    }    

    if(adj){/* migration */
        for(it=nt-1;it>=0;it--){
            sf_warning("Migration: %d;",it);

	    boundary(u0,u1);
	    adjoint(u1,ud);

#ifdef _OPENMP
#pragma omp parallel for    \
    private(i2,i1)       \
    shared(u0,u1,u2,ud,n1,n2,dd,n0,it)
#endif            
            for(i2=0;i2<n2;i2++){
                for(i1=0;i1<n1;i1++){
                    
                    u2[i2][i1]=2*u1[i2][i1]
                              -u0[i2][i1]
                              +ud[i2][i1];
                         
                    u0[i2][i1]=u1[i2][i1];
                    u1[i2][i1]=u2[i2][i1];
                }
                
                /* inject data */
                u1[i2][n0]+=dd[i2][it];            
            }            
        }            
        sf_warning(".");
        /* output image */
        sf_floatwrite(u1[0],n12,out);    	
    }else{/* modeling */
    	/* read reflector */
    	sf_floatread(u1[0],n12,in);
    	    
    	for(it=0;it<nt;it++){
    	    sf_warning("Modeling: %d;",it);
	  
            /* record data */
	    for(i2=0;i2<n2;i2++)
	        dd[i2][it]+=u1[i2][n0];
	
	    laplacian(u0,u1,u2,ud);
	    boundary(u1,u2);

#ifdef _OPENMP
#pragma omp parallel for    \
    private(i2,i1)       \
    shared(u0,u1,u2,n1,n2)
#endif
            for(i2=0;i2<n2;i2++){
	        for(i1=0;i1<n1;i1++){
	            u0[i2][i1]=u1[i2][i1];
	            u1[i2][i1]=u2[i2][i1];
	        }
	    }
    	}
    	sf_warning(".");
    	/* output data */
    	sf_floatwrite(dd[0],nt*n2,out);
    }
   
   free(*vv);free(vv);
   free(abc);    
    
   exit (0);
}

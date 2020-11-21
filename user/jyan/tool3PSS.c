#include <rsf.h>
#include <math.h>


#include "eigv.h"
#include "dsyevv3.h"
#include "dsyevc3.h"

#define  KMAP(i,n) (i<n/2.) ? SF_PI*i/(n/2.) : SF_PI*(-n+i)/(n/2.)
#define  KMAPK(i,n) (  ( i-n/2 )*SF_PI/(maxn/2.)       )    
#define TAPER(k) (k!=0 ?			\
		  8./5.  *sin(  k)/k -		\
		  2./5.  *sin(2*k)/k +		\
		  8./105.*sin(3*k)/k -		\
		  1./140.*sin(4*k)/k : 1)

#define TAPERG(kmag) exp(-(kmag*kmag)/1.5)/2.

#define  MAX2(a,b)   a>b? a: b
#define  MY_MAX(a,b,c)  MAX2(a,b) >c ? MAX2(a,b) :c
/*------------------------------------------------------------*/
/*2D*/
typedef struct wfs *wfs2d;
/*^*/

struct wfs{
    sf_complex ***temp;
    float        *eigv; /* eigenvalues */
    float        *polr; /* eigenvector */ 
    float       **ctfl; /* Christoffel matrix */
    float **UPz, **UPx;    
    sf_fft3d   ftz,   ftx;
    float **Cos, **Sin;
};/*^*/



wfs2d wfsep_init(sf_axis  ax,
		 sf_axis  az)
/*< init wavefield separator >*/
{
    wfs2d wfs;
    wfs = (wfs2d) sf_alloc(1,sizeof(*wfs));

    wfs->temp = sf_complexalloc3(sf_n(az),sf_n(ax),1);
    wfs->eigv = sf_floatalloc (2);           /* eigenvalues */
    wfs->polr = sf_floatalloc (2);           /* eigenvector */
    wfs->ctfl = sf_floatalloc2(2,2);         /* Christoffel matrix */

    wfs->UPz  = sf_floatalloc2(sf_n(az),sf_n(ax));
    wfs->UPx  = sf_floatalloc2(sf_n(az),sf_n(ax));

    wfs->ftz=sf_fft3a1_init(sf_n(az),sf_n(ax),1);
    wfs->ftx=sf_fft3a2_init(sf_n(az),sf_n(ax),1);

    wfs->Cos  = sf_floatalloc2(sf_n(az),sf_n(ax));
    wfs->Sin  = sf_floatalloc2(sf_n(az),sf_n(ax));

    return wfs;
}

/*------------------------------------------------------------*/
/*3D*/
typedef struct wfs3 *wfs3d;
/*^*/

struct wfs3{   
    sf_complex ***temp;
    float        *eigv; /* eigenvalues */
    float       **ctfl; /* Christoffel matrix */ 
    float          *up;/*eigenvectors */
    float       ***upx,***upy,***upz;
    sf_fft3d          ftx,fty,ftz;
};
/*^*/
    
wfs3d wfsep3_init(sf_axis  ax,
		  sf_axis  ay,
		  sf_axis  az)
/*< init wavefield separator >*/
{
    wfs3d wfs3;
    wfs3 = (wfs3d) sf_alloc(1,sizeof(*wfs3));

    wfs3->temp = sf_complexalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    wfs3->eigv = sf_floatalloc (3);           /* eigenvalues */
    wfs3->up   = sf_floatalloc (3);           /* P eigenvector */
    wfs3->ctfl = sf_floatalloc2(3,3);         /* Christoffel matrix */
    
    wfs3->upx  = sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    wfs3->upy  = sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    wfs3->upz  = sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    
    wfs3->ftz=sf_fft3a1_init(sf_n(az),sf_n(ax),sf_n(ay));
    wfs3->ftx=sf_fft3a2_init(sf_n(az),sf_n(ax),sf_n(ay));
    wfs3->fty=sf_fft3a3_init(sf_n(az),sf_n(ax),sf_n(ay));
    return wfs3;
}


void wfsep3_close(wfs3d wfs3)
/*< close wavefield separator >*/
{
    sf_fft3a1_close(wfs3->ftz);
    sf_fft3a2_close(wfs3->ftx);
    sf_fft3a3_close(wfs3->fty);
}




/*------------------------------------------------------------*/
/*Program  to calculate eigenvectors and eigenvalues*/
void wfsepK3(float ***pxdel, float ***pydel, float *** pzdel,
	     float ***vxdel, float ***vydel, float *** vzdel,
	     float ***hxdel, float ***hydel, float *** hzdel,	     
	     sf_axis  ax,
	     sf_axis  ay,
	     sf_axis  az,
	     float   c11,
	     float   c12,
	     float   c13,
	     float   c22,
	     float   c23,
	     float   c33,
	     float   c44,
	     float   c55,
	     float   c66,
	     wfs3d   wfs3
    )
/*<test>*/
{
    
    float a11,a22,a33,a12,a13,a23;
    float k=0.0;
    float upx,upy,upz;

    
    int    jx, jy, jz;
    int    nx, ny, nz;
    float  kx, ky, kz;
   
    float tmp;
    int   maxn;

    double A[3][3],w[3],Q[3][3];

    nx = sf_n(ax);
    ny = sf_n(ay);
    nz = sf_n(az);
    sf_warning("%d,%d,%d",nx,ny,nz);
    maxn=MY_MAX(nx,ny,nz);
    /*initialize here*/

    for(jy=0;jy<ny;jy++){
	for(jx=0;jx<nx;jx++){  
	    for(jz=0;jz<nz;jz++){ 
		kz = KMAPK(jz,nz);ky = KMAPK(jy,ny); kx = KMAPK(jx,nx); 
		wfs3->upx[jy][jx][jz]=kx;
		wfs3->upy[jy][jx][jz]=ky;
		wfs3->upz[jy][jx][jz]=kz;
	    }
	}
    }

    for(jy=0;jy<ny;jy++){
	for(jx=0;jx<nx;jx++){
	    for(jz=0;jz<nz;jz++){
		kz = KMAPK(jz,nz);ky = KMAPK(jy,ny); kx = KMAPK(jx,nx);
		if(kz!=0.0){
		    
		    k=sqrt(kx*kx+ky*ky+kz*kz);
		    a11=c11*kx*kx+
			c66*ky*ky +
			c55*kz*kz;
		    
		    a22=c66*kx*kx+
			c22*ky*ky +
			c44*kz*kz;
	    
		    a33=c55*kx*kx+
			c44*ky*ky +
			c33*kz*kz;
 
		    a12= (c12+c66)*kx*ky;
		    a13= (c13+c55)*kx*kz;
		    a23= (c23+c44)*ky*kz;
	    

		    wfs3->ctfl[0][0] = a11;
		    wfs3->ctfl[0][1] = a12;
		    wfs3->ctfl[0][2] = a13;
		    wfs3->ctfl[1][0] = a12;
		    wfs3->ctfl[1][1] = a22;
		    wfs3->ctfl[1][2] = a23;
		    wfs3->ctfl[2][0] = a13;
		    wfs3->ctfl[2][1] = a23;
		    wfs3->ctfl[2][2] = a33;

/* 		    eigval3(wfs3->ctfl,wfs3->eigv); */
/* 		    Leigvec3(wfs3->ctfl,wfs3->eigv,wfs3->up); */

		    A[0][0] = a11;
		    A[0][1] = a12;
		    A[0][2] = a13;
		    A[1][0] = a12;
		    A[1][1] = a22;
		    A[1][2] = a23;
		    A[2][0] = a13;
		    A[2][1] = a23;
		    A[2][2] = a33;
		    
		    dsyevc3(A,w);
		    dsyevv3(A, Q, w );



		}
		
		   
				
/* 		/\*normalized eigenvectors that correspond to largest eigenvalue*\/ */
/* 		upx=wfs3->up[0]; */
/* 		upy=wfs3->up[1];upz=wfs3->up[2]; */
		
/* 		/\* get the closest direction to k *\/ */
/* 		if(upx*kx + upy*ky+ upz*kz < 0.) { */
/* 		    upx=-wfs3->up[0]; */
/* 		    upy=-wfs3->up[1];upz=-wfs3->up[2]; */
/* 		} */
	
/* 		wfs3->upx[jy][jx][jz] = upx*k; */
/* 		wfs3->upy[jy][jx][jz] = upy*k;  */
/* 		wfs3->upz[jy][jx][jz] = upz*k; */
		

		upx=Q[0][0];
		upy=Q[1][0];
		upz=Q[2][0];



		/* get the closest direction to k */
		if(upx*kx + upy*ky+ upz*kz < 0.) {
		    upx=-Q[0][0];
		    upy=-Q[1][0];
		    upz=-Q[2][0];
		}
	
		wfs3->upx[jy][jx][jz] = upx*k;
		wfs3->upy[jy][jx][jz] = upy*k;
		wfs3->upz[jy][jx][jz] = upz*k;


		
		
	    }
	}
    } 
    
    for(jy=0;jy<ny;jy++){
	for(jx=0;jx<nx;jx++){  
	    for(jz=0;jz<nz;jz++){ 
		kz = KMAPK(jz,nz);ky = KMAPK(jy,ny); kx = KMAPK(jx,nx);
		tmp=sqrt(kx*kx+ky*ky+kz*kz);
		pxdel[jy][jx][jz]=wfs3->upx[jy][jx][jz]*TAPERG(tmp);
		pydel[jy][jx][jz]=wfs3->upy[jy][jx][jz]*TAPERG(tmp);
		pzdel[jy][jx][jz]=wfs3->upz[jy][jx][jz]*TAPERG(tmp);	

		hxdel[jy][jx][jz]=-ky*TAPERG(tmp);
		hydel[jy][jx][jz]= kx*TAPERG(tmp);
		hzdel[jy][jx][jz]= 0;

		vxdel[jy][jx][jz]=-kx*pzdel[jy][jx][jz];
		vydel[jy][jx][jz]=-ky*pzdel[jy][jx][jz];
		vzdel[jy][jx][jz]= kx*pxdel[jy][jx][jz]+ky*pydel[jy][jx][jz];

	    }
	}
    }

  
}

/*------------------------------------------------------------*/
void wfsep3(float ***pxdel, float ***pydel, float *** pzdel,
	    float ***vxdel, float ***vydel, float *** vzdel,
	    float ***hxdel, float ***hydel, float *** hzdel,	     
	    sf_axis  ax,
	    sf_axis  ay,
	    sf_axis  az,
	    float   c11,
	    float   c12,
	    float   c13,
	    float   c22,
	    float   c23,
	    float   c33,
	    float   c44,
	    float   c55,
	    float   c66,
	    wfs3d   wfs3
    )
/*<test>*/
{
    
    float a11,a22,a33,a12,a13,a23;
    float k=0.0;
    float upx,upy,upz;

    
    int    jx, jy, jz;
    int    nx, ny, nz;
    float  kx, ky, kz;
   
    float tmp=0.0;
/*    float maxn; */

    double A[3][3],w[3],Q[3][3];

    nx = sf_n(ax);
    ny = sf_n(ay);
    nz = sf_n(az);
    sf_warning("%d,%d,%d",nx,ny,nz);
    /*  maxn=MAX(nx,ny,nz); */
    /*initialize here*/
    
    for(jy=0;jy<ny;jy++){
	for(jx=0;jx<nx;jx++){  
	    for(jz=0;jz<nz;jz++){ 
		kz = KMAP(jz,nz);ky = KMAP(jy,ny); kx = KMAP(jx,nx); 
		wfs3->upx[jy][jx][jz]=kx;
		wfs3->upy[jy][jx][jz]=ky;
		wfs3->upz[jy][jx][jz]=kz;
	    }
	}
    }
 

    for(jy=0;jy<ny;jy++){
	for(jx=0;jx<nx;jx++){
	    for(jz=0;jz<nz;jz++){
    		kz = KMAP(jz,nz);ky = KMAP(jy,ny); kx = KMAP(jx,nx);
		if(kz!=0.0){
		    k=sqrt(kx*kx+ky*ky+kz*kz);
		    kx/=k;ky/=k;kz/=k;
					
		    a11=  c11*kx*kx+
			c66*ky*ky +
			c55*kz*kz;
	    
		    a22=  c66*kx*kx+
			c22*ky*ky +
			c44*kz*kz;
	    
		    a33=  c55*kx*kx+
			c44*ky*ky +
			c33*kz*kz;
 
		    a12= (c12+c66)*kx*ky;
		    a13= (c13+c55)*kx*kz;
		    a23= (c23+c44)*ky*kz;
	    


		    wfs3->ctfl[0][0] = a11;
		    wfs3->ctfl[0][1] = a12;
		    wfs3->ctfl[0][2] = a13;
		    wfs3->ctfl[1][0] = a12;
		    wfs3->ctfl[1][1] = a22;
		    wfs3->ctfl[1][2] = a23;
		    wfs3->ctfl[2][0] = a13;
		    wfs3->ctfl[2][1] = a23;
		    wfs3->ctfl[2][2] = a33;

/*		    eigval3(wfs3->ctfl,wfs3->eigv);*/
/*		    Leigvec3(wfs3->ctfl,wfs3->eigv,wfs3->up);*/


		    A[0][0] = a11;
		    A[0][1] = a12;
		    A[0][2] = a13;
		    A[1][0] = a12;
		    A[1][1] = a22;
		    A[1][2] = a23;
		    A[2][0] = a13;
		    A[2][1] = a23;
		    A[2][2] = a33;
		    
		    dsyevc3(A,w);
		    dsyevv3(A, Q, w );


		}

		/*normalized eigenvectors that correspond to largest eigenvalue*/
	/* 	upx=wfs3->up[0]; */
/* 		upy=wfs3->up[1];upz=wfs3->up[2]; */
		
/* 		/\* get the closest direction to k *\/ */
/* 		if(upx*kx + upy*ky+ upz*kz < 0.) { */
/* 		    upx=-wfs3->up[0]; */
/* 		    upy=-wfs3->up[1];upz=-wfs3->up[2]; */
/* 		} */
		upx=Q[0][0];
		upy=Q[1][0];
		upz=Q[2][0];



		/* get the closest direction to k */
		if(upx*kx + upy*ky+ upz*kz < 0.) {
		    upx=-Q[0][0];
		    upy=-Q[1][0];
		    upz=-Q[2][0];
		}
	
		wfs3->upx[jy][jx][jz] = upx*k;
		wfs3->upy[jy][jx][jz] = upy*k;
		wfs3->upz[jy][jx][jz] = upz*k;

	
	    }
	}
    }
   
/*------------------------------------------------------------*/
    /* 
     * Z derivative, P
     */
    for(jy=0;jy<ny;jy++){ky = KMAP(jy,ny);
	for    (jx=0;jx<nx;jx++){kx = KMAP(jx,nx);
	    for(jz=0;jz<nz;jz++){ kz = KMAP(jz,nz);
		tmp=sqrt(kx*kx+kz*kz+ky*ky);
		wfs3->temp[jy][jx][jz]=sf_cmplx( 0, wfs3->upz[jy][jx][jz] * TAPERG(tmp) );
	    }
	}   
    }

  
    sf_cnt3a3(wfs3->temp,wfs3->fty);
    sf_cnt3a2(wfs3->temp,wfs3->ftx);
    sf_cnt3a1(wfs3->temp,wfs3->ftz);
    sf_fft3a3(true,(kiss_fft_cpx***) wfs3->temp,wfs3->fty);sf_warning("stop a");
    sf_fft3a2(true,(kiss_fft_cpx***) wfs3->temp,wfs3->ftx);
    sf_fft3a1(true,(kiss_fft_cpx***) wfs3->temp,wfs3->ftz);    
    
    for(jy=0;jy<ny;jy++){
	for(    jx=0;jx<nx;jx++){
	    for(jz=0;jz<nz;jz++){
		pzdel[jy][jx][jz]=crealf(wfs3->temp[jy][jx][jz]);
	    }
	}
    } 

/* SV  */
    for(jy=0;jy<ny;jy++){ky = KMAP(jy,ny);
	for    (jx=0;jx<nx;jx++){kx = KMAP(jx,nx);
	    for(jz=0;jz<nz;jz++){ kz = KMAP(jz,nz);
		tmp=sqrt(kx*kx+kz*kz+ky*ky);
		wfs3->temp[jy][jx][jz]=sf_cmplx(  (wfs3->upx[jy][jx][jz]*kx +wfs3->upy[jy][jx][jz] *ky) * TAPERG(tmp),0 );
	    }
	}
    }
  
    sf_cnt3a3(wfs3->temp,wfs3->fty);
    sf_cnt3a2(wfs3->temp,wfs3->ftx);
    sf_cnt3a1(wfs3->temp,wfs3->ftz);
    sf_fft3a3(true,(kiss_fft_cpx***) wfs3->temp,wfs3->fty);sf_warning("stop a");
    sf_fft3a2(true,(kiss_fft_cpx***) wfs3->temp,wfs3->ftx);
    sf_fft3a1(true,(kiss_fft_cpx***) wfs3->temp,wfs3->ftz);
    
    for(jy=0;jy<ny;jy++){
	for(    jx=0;jx<nx;jx++){
	    for(jz=0;jz<nz;jz++){
		vzdel[jy][jx][jz]=crealf(wfs3->temp[jy][jx][jz]);
	    }
	}
    }
 /* Sh  */
  
    for(jy=0;jy<ny;jy++){
	for(jx=0;jx<nx;jx++){
	    for(jz=0;jz<nz;jz++){
		hzdel[jy][jx][jz]=0.;
	    }
	}
    }


    /*------------------------------------------------------------*/
  

    /*------------------------------------------------------------*/
    /* 
     * X derivative, P
     */
    for(jy=0;jy<ny;jy++){ ky = KMAP(jy,ny);
	for    (jx=0;jx<nx;jx++){ kx = KMAP(jx,nx);
	    for(jz=0;jz<nz;jz++){ kz = KMAP(jz,nz);
		tmp=sqrt(kx*kx+kz*kz+ky*ky);
		wfs3->temp[jy][jx][jz]=sf_cmplx( 0, wfs3->upx[jy][jx][jz] * TAPERG(tmp) );
		
	    }
	}
    }
 
    sf_cnt3a3(wfs3->temp,wfs3->fty);
    sf_cnt3a2(wfs3->temp,wfs3->ftx);
    sf_cnt3a1(wfs3->temp,wfs3->ftz);
    sf_fft3a3(true,(kiss_fft_cpx***) wfs3->temp,wfs3->fty);
    sf_fft3a2(true,(kiss_fft_cpx***) wfs3->temp,wfs3->ftx);
    sf_fft3a1(true,(kiss_fft_cpx***) wfs3->temp,wfs3->ftz);
    for(jy=0;jy<ny;jy++){
	for(    jx=0;jx<nx;jx++){
	    for(jz=0;jz<nz;jz++){
		pxdel[jy][jx][jz]=crealf(wfs3->temp[jy][jx][jz]);
	    }
	}
    }
/*SV   */
 for(jy=0;jy<ny;jy++){ ky = KMAP(jy,ny);
	for    (jx=0;jx<nx;jx++){ kx = KMAP(jx,nx);
	    for(jz=0;jz<nz;jz++){ kz = KMAP(jz,nz);
		tmp=sqrt(kx*kx+kz*kz+ky*ky);
		wfs3->temp[jy][jx][jz]=sf_cmplx(  -wfs3->upz[jy][jx][jz]*kx*TAPERG(tmp), 0);
		
	    }
	}
    }
 
    sf_cnt3a3(wfs3->temp,wfs3->fty);
    sf_cnt3a2(wfs3->temp,wfs3->ftx);
    sf_cnt3a1(wfs3->temp,wfs3->ftz);
    sf_fft3a3(true,(kiss_fft_cpx***) wfs3->temp,wfs3->fty);
    sf_fft3a2(true,(kiss_fft_cpx***) wfs3->temp,wfs3->ftx);
    sf_fft3a1(true,(kiss_fft_cpx***) wfs3->temp,wfs3->ftz);
    for(jy=0;jy<ny;jy++){
	for(    jx=0;jx<nx;jx++){
	    for(jz=0;jz<nz;jz++){
		vxdel[jy][jx][jz]=crealf(wfs3->temp[jy][jx][jz]);


	    }
	}
    }
/*SH*/
 for(jy=0;jy<ny;jy++){ ky = KMAP(jy,ny);
	for    (jx=0;jx<nx;jx++){ kx = KMAP(jx,nx);
	    for(jz=0;jz<nz;jz++){ kz = KMAP(jz,nz);
		tmp=sqrt(kx*kx+kz*kz+ky*ky);
		wfs3->temp[jy][jx][jz]=sf_cmplx( 0, -ky * TAPERG(tmp) );
		
	    }
	}
    }
 
    sf_cnt3a3(wfs3->temp,wfs3->fty);
    sf_cnt3a2(wfs3->temp,wfs3->ftx);
    sf_cnt3a1(wfs3->temp,wfs3->ftz);
    sf_fft3a3(true,(kiss_fft_cpx***) wfs3->temp,wfs3->fty);
    sf_fft3a2(true,(kiss_fft_cpx***) wfs3->temp,wfs3->ftx);
    sf_fft3a1(true,(kiss_fft_cpx***) wfs3->temp,wfs3->ftz);
    for(jy=0;jy<ny;jy++){
	for(    jx=0;jx<nx;jx++){
	    for(jz=0;jz<nz;jz++){
		hxdel[jy][jx][jz]=crealf(wfs3->temp[jy][jx][jz]);
	    }
	}
    }

    /*------------------------------------------------------------*/
    /*------------------------------------------------------------*/
    /* 
     * Y derivative,P
     */
    for(jy=0;jy<ny;jy++){ky = KMAP(jy,ny);
	for    (jx=0;jx<nx;jx++){ kx = KMAP(jx,nx);
	    for(jz=0;jz<nz;jz++){kz = KMAP(jz,nz);
		tmp=sqrt(kx*kx+kz*kz+ky*ky);
		wfs3->temp[jy][jx][jz]=sf_cmplx( 0, wfs3->upy[jy][jx][jz] * TAPERG(tmp)/tmp );
		
	    }
	}
    }

    sf_cnt3a3(wfs3->temp,wfs3->fty);
    sf_cnt3a2(wfs3->temp,wfs3->ftx);
    sf_cnt3a1(wfs3->temp,wfs3->ftz);
    sf_fft3a3(true,(kiss_fft_cpx***) wfs3->temp,wfs3->fty);
    sf_fft3a2(true,(kiss_fft_cpx***) wfs3->temp,wfs3->ftx);
    sf_fft3a1(true,(kiss_fft_cpx***) wfs3->temp,wfs3->ftz);
    for(jy=0;jy<ny;jy++){
	for(    jx=0;jx<nx;jx++){
	    for(jz=0;jz<nz;jz++){
		pydel[jy][jx][jz]=crealf(wfs3->temp[jy][jx][jz]);
	    }
	}
    }

/* v  */
for(jy=0;jy<ny;jy++){ky = KMAP(jy,ny);
	for    (jx=0;jx<nx;jx++){ kx = KMAP(jx,nx);
	    for(jz=0;jz<nz;jz++){kz = KMAP(jz,nz);
		tmp=sqrt(kx*kx+kz*kz+ky*ky);
		wfs3->temp[jy][jx][jz]=sf_cmplx(  -wfs3->upz[jy][jx][jz]*ky* TAPERG(tmp) ,0);
		
	    }
	}
    }

    sf_cnt3a3(wfs3->temp,wfs3->fty);
    sf_cnt3a2(wfs3->temp,wfs3->ftx);
    sf_cnt3a1(wfs3->temp,wfs3->ftz);
    sf_fft3a3(true,(kiss_fft_cpx***) wfs3->temp,wfs3->fty);
    sf_fft3a2(true,(kiss_fft_cpx***) wfs3->temp,wfs3->ftx);
    sf_fft3a1(true,(kiss_fft_cpx***) wfs3->temp,wfs3->ftz);
    for(jy=0;jy<ny;jy++){
	for(    jx=0;jx<nx;jx++){
	    for(jz=0;jz<nz;jz++){
		vydel[jy][jx][jz]=crealf(wfs3->temp[jy][jx][jz]);
	    }
	}
    }


/*h*/
   for(jy=0;jy<ny;jy++){ky = KMAP(jy,ny);
	for    (jx=0;jx<nx;jx++){ kx = KMAP(jx,nx);
	    for(jz=0;jz<nz;jz++){kz = KMAP(jz,nz);
		tmp=sqrt(kx*kx+kz*kz+ky*ky);
		wfs3->temp[jy][jx][jz]=sf_cmplx( 0,kx * TAPERG(tmp) );
	    }
	}
    }

    sf_cnt3a3(wfs3->temp,wfs3->fty);
    sf_cnt3a2(wfs3->temp,wfs3->ftx);
    sf_cnt3a1(wfs3->temp,wfs3->ftz);
    sf_fft3a3(true,(kiss_fft_cpx***) wfs3->temp,wfs3->fty);
    sf_fft3a2(true,(kiss_fft_cpx***) wfs3->temp,wfs3->ftx);
    sf_fft3a1(true,(kiss_fft_cpx***) wfs3->temp,wfs3->ftz);
    for(jy=0;jy<ny;jy++){
	for(    jx=0;jx<nx;jx++){
	    for(jz=0;jz<nz;jz++){
		hydel[jy][jx][jz]=crealf(wfs3->temp[jy][jx][jz]);
	    }
	}
    }
 

}

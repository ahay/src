#include <rsf.h>
#include <math.h>

#include "dsyevv3.h"
#include "dsyevc3.h"

#define  KMAP(i,n) (i<n/2.) ? SF_PI*i/(n/2.) : SF_PI*(-n+i)/(n/2.)
#define  KMAPK(i,n)  (  ( i-n/2 )*SF_PI/(n/2.)   ) 
#define  TAPER(k) (k!=0 ?			\
		  8./5.  *sin(  k)/k -		\
		  2./5.  *sin(2*k)/k +		\
		  8./105.*sin(3*k)/k -		\
		  1./140.*sin(4*k)/k : 1)

#define TAPERG(kmag,sigma) exp(-(kmag*kmag)/(1.5*sigma*sigma)  )/2.

#define  MAX2(a,b)   a>b? a: b
#define  MAX(a,b,c)  (MAX2(a,b)) >c ? (MAX2(a,b)) :c

#define  MIN2(a,b)   a<b? a: b
#define  MIN(a,b,c)  (MIN2(a,b)) <c ? (MIN2(a,b)) :c

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
};
/*^*/
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
    float        *up;/*eigenvectors */
    float      ***upx,***upy,***upz;
    sf_fft3d         ftx,fty,ftz;
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
void wfsepK3(float ***xdel, float ***ydel, float *** zdel,
	     sf_axis  ax,
	     sf_axis  ay,
	     sf_axis  az,
	     float c11, float c12, float c13, float c14, float c15, float  c16,
	     float c22, float c23, float c24, float c25, float c26,
	     float c33, float c34, float c35, float c36,
	     float c44, float c45, float c46,
	     float c55, float c56,
	     float c66, 
	     wfs3d wfs3)
/*< test >*/
{
    
    float a11,a22,a33,a12,a13,a23;
    float k=0.0;
    float upx,upy,upz;

    
    int    jx=0, jy=0, jz=0;
    int    nx, ny, nz;
    float  kx=0.0, ky=0.0, kz=0.0;
   
    float tmp,sigma=1.0;
    /* int   maxn; */

    double A[3][3],w[3],Q[3][3];

    nx = sf_n(ax);
    ny = sf_n(ay);
    nz = sf_n(az);
    sf_warning("%d,%d,%d",nx,ny,nz);
    /* maxn=MAX(nx,ny,nz); */
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
		if( kz!=0. ){
		    k=sqrt(kx*kx+ky*ky+kz*kz);
		    a11=c11*kx*kx+2*c16*kx*ky+
			c66*ky*ky+2*c15*kx*kz+
			c55*kz*kz+2*c56*ky*kz;
		    
		    a22=c66*kx*kx+    2*c26*kx*ky+
			c22*ky*ky+(c45+c46)*kx*kz+
			c44*kz*kz+    2*c24*ky*kz;
	    
		    a33=c55*kx*kx+ 2*c45*kx*ky+
			c44*ky*ky+ 2*c35*kx*kz+
			c33*kz*kz+ 2*c34*ky*kz;
 
		    a12=c16*kx*kx+(c12+c66)*kx*ky+
			c26*ky*ky+(c14+c56)*kx*kz+
			c45*kz*kz+(c25+c46)*ky*kz;
		    a13=c15*kx*kx+(c14+c56)*kx*ky+
			c46*ky*ky+(c13+c55)*kx*kz+
			c35*kz*kz+(c36+c45)*ky*kz;
		    a23=c56*kx*kx+(c25+c46)*kx*ky+
			c24*ky*ky+(c36+c45)*kx*kz+
			c34*kz*kz+(c23+c44)*ky*kz;


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
/*		upx=wfs3->up[0];*/
/*		upy=wfs3->up[1];upz=wfs3->up[2];*/
		
		/* get the closest direction to k */
/*		if(upx*kx + upy*ky+ upz*kz < 0.) {*/
/*		    upx=-wfs3->up[0];*/
/*		    upy=-wfs3->up[1];upz=-wfs3->up[2];*/
/*		}*/

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
    
  

/*    sigma=(float)(MIN(nx,ny,nz))/(float)(MAX(nx,ny,nz));*/

    for(jy=0;jy<ny;jy++){
	for(jx=0;jx<nx;jx++){  
	    for(jz=0;jz<nz;jz++){ 
		kz = KMAPK(jz,nz);ky = KMAPK(jy,ny); kx = KMAPK(jx,nx);
		tmp=sqrt(kx*kx+ky*ky+kz*kz);
		xdel[jy][jx][jz]=wfs3->upx[jy][jx][jz]*TAPERG(tmp,sigma);
		ydel[jy][jx][jz]=wfs3->upy[jy][jx][jz]*TAPERG(tmp,sigma);
		zdel[jy][jx][jz]=wfs3->upz[jy][jx][jz]*TAPERG(tmp,sigma);	
	    }
	}
    }

sf_warning("kx ky kz = %f %f %f",kx,ky,kz);  
kz = KMAPK(jz,nz);ky = KMAPK(jy,ny); kx = KMAPK(jx,nx);

sf_warning("kx ky kz = %f %f %f",kx,ky,kz);

sf_warning("test1 =%f",KMAPK(jz,nz));
sf_warning("test2 =%f",KMAPK(jx,nx));
sf_warning("test3 =%f",KMAPK(jy,ny));


}




/*------------------------------------------------------------*/
/*Program  to calculate eigenvectors and eigenvalues*/
void wfsep3(float ***xdel, float ***ydel, float *** zdel,
	     sf_axis  ax,
	     sf_axis  ay,
	     sf_axis  az,
	     float c11, float c12, float c13, float c14, float c15, float  c16,
	     float c22, float c23, float c24, float c25, float c26,
	     float c33, float c34, float c35, float c36,
	     float c44, float c45, float c46,
	     float c55, float c56,
	     float c66, 
	     wfs3d wfs3)
/*< test >*/
{
    
    float a11,a22,a33,a12,a13,a23;
    float k=0.0;
    float upx,upy,upz;

    
    int    jx, jy, jz;
    int    nx, ny, nz;
    float  kx, ky, kz;
   
    float tmp=0.0,sigma=1.0;
    /* int   maxn; */

    double A[3][3],w[3],Q[3][3];

    nx = sf_n(ax);
    ny = sf_n(ay);
    nz = sf_n(az);
    sf_warning("%d,%d,%d",nx,ny,nz);
    /* maxn=MAX(nx,ny,nz); */
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
		if( kz!=0. ){
		    k=sqrt(kx*kx+ky*ky+kz*kz);
/*		    kx/=k;ky/=k;kz/=k;*/
		    a11=c11*kx*kx+2*c16*kx*ky+
			c66*ky*ky+2*c15*kx*kz+
			c55*kz*kz+2*c56*ky*kz;
		    
		    a22=c66*kx*kx+    2*c26*kx*ky+
			c22*ky*ky+(c45+c46)*kx*kz+
			c44*kz*kz+    2*c24*ky*kz;
	    
		    a33=c55*kx*kx+ 2*c45*kx*ky+
			c44*ky*ky+ 2*c35*kx*kz+
			c33*kz*kz+ 2*c34*ky*kz;
 
		    a12=c16*kx*kx+(c12+c66)*kx*ky+
			c26*ky*ky+(c14+c56)*kx*kz+
			c45*kz*kz+(c25+c46)*ky*kz;
		    a13=c15*kx*kx+(c14+c56)*kx*ky+
			c46*ky*ky+(c13+c55)*kx*kz+
			c35*kz*kz+(c36+c45)*ky*kz;
		    a23=c56*kx*kx+(c25+c46)*kx*ky+
			c24*ky*ky+(c36+c45)*kx*kz+
			c34*kz*kz+(c23+c44)*ky*kz;


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
/*		upx=wfs3->up[0];*/
/*		upy=wfs3->up[1];upz=wfs3->up[2];*/
		
		/* get the closest direction to k */
/*		if(upx*kx + upy*ky+ upz*kz < 0.) {*/
/*		    upx=-wfs3->up[0];*/
/*		    upy=-wfs3->up[1];upz=-wfs3->up[2];*/
/*		}*/
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
  

    sigma=(MIN(nx,ny,nz)) / (MAX(nx,ny,nz));
   
/*------------------------------------------------------------*/
    /* 
     * Z derivative
     */
    for(jy=0;jy<ny;jy++){ky = KMAP(jy,ny);
	for    (jx=0;jx<nx;jx++){kx = KMAP(jx,nx);
	    for(jz=0;jz<nz;jz++){ kz = KMAP(jz,nz);
		tmp=sqrt(kx*kx+kz*kz+ky*ky);
		wfs3->temp[jy][jx][jz]=sf_cmplx( 0, wfs3->upz[jy][jx][jz] * TAPERG(tmp,sigma) );
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
		zdel[jy][jx][jz]=crealf(wfs3->temp[jy][jx][jz]);
	    }
	}
    }


    /*------------------------------------------------------------*/
  

    /*------------------------------------------------------------*/
    /* 
     * X derivative
     */
    for(jy=0;jy<ny;jy++){ ky = KMAP(jy,ny);
	for    (jx=0;jx<nx;jx++){ kx = KMAP(jx,nx);
	    for(jz=0;jz<nz;jz++){ kz = KMAP(jz,nz);
		tmp=sqrt(kx*kx+kz*kz+ky*ky);
		wfs3->temp[jy][jx][jz]=sf_cmplx( 0, wfs3->upx[jy][jx][jz] * TAPERG(tmp,sigma) );
		
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
		xdel[jy][jx][jz]=crealf(wfs3->temp[jy][jx][jz]);
	    }
	}
    }
   
    /*------------------------------------------------------------*/
    /*------------------------------------------------------------*/
    /* 
     * Y derivative
     */
    for(jy=0;jy<ny;jy++){ky = KMAP(jy,ny);
	for    (jx=0;jx<nx;jx++){ kx = KMAP(jx,nx);
	    for(jz=0;jz<nz;jz++){kz = KMAP(jz,nz);
		tmp=sqrt(kx*kx+kz*kz+ky*ky);
		wfs3->temp[jy][jx][jz]=sf_cmplx( 0, wfs3->upy[jy][jx][jz] * TAPERG(tmp,sigma) );
		
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
		ydel[jy][jx][jz]=crealf(wfs3->temp[jy][jx][jz]);
	    }
	}
    }

   
 
  
}

/* Anisotropy routines. */

/*
  Copyright (C) 2006 The Board of Trustees of Stanford University
  
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

#include "ani.h"
#include "mymath.h"
#include "arrayindex.h"

void coe(float* a,float phi,float ep,float dl,float w_vp,float f,float kx)
{
    float a0,a1,a2,a3,a4;
    float fm1,epdl,sin2phi,sin4phi,sinphi,cosphi,cos2phi,b1,b2;

    fm1=f-1;
    epdl=ep-dl;
    sin2phi=sinf(2.0*phi);
    cos2phi=cosf(2.0*phi);
    sin4phi=sinf(4.0*phi);
    sinphi=sinf(phi);
    cosphi=cosf(phi);

    a4=fm1+2.0*ep*fm1*powf(sinphi,2)-f/2.0*epdl*powf(sin2phi,2);
    a3=(f*epdl*sin4phi-2.0*ep*fm1*sin2phi)*kx;
    b2=f*epdl*powf(sin2phi,2)+2.0*fm1*(1+ep)-2.0*f*epdl*powf(cos2phi,2);
    a2=b2*kx*kx+w_vp*w_vp*(2.0+2.0*ep*powf(sinphi,2)-f);
    b1=2.0*ep*(-fm1)*sin2phi-f*epdl*sin4phi;
    a1=b1*powf(kx,3)-2.0*ep*sin2phi*w_vp*w_vp*kx;
    a0=(2.0+2.0*ep*powf(cosphi,2)-f)*w_vp*w_vp*kx*kx-powf(w_vp,4)-(-fm1*(1.0+2.0*ep*powf(cosphi,2))+f/2.0*epdl*powf(sin2phi,2))*powf(kx,4);
    a[0]=a0; a[1]=a1; a[2]=a2; a[3]=a3; a[4]=a4;
}

void coe3d(float *a, float phi,float ep,float dl,float w_vp,float f, float kx,float ky)
{
    float a0,a1,a2,a3,a4;
    float sinphi,cosphi,fm1,epdl,sin2phi,cos2phi,sin4phi;
    float A,B,C,D,E,F;
    sinphi=-sinf(phi); cosphi=cosf(phi);
    fm1=f-1; epdl=ep-dl; sin2phi=-sinf(2.0*phi); cos2phi=cosf(2.0*phi); sin4phi=-sinf(4.0*phi);
    //printf("sin=%f\n",sinphi);
    A=f-1;
    B=(f-1.0)*(1.0+2.0*ep);
    C=2.0*( (f-1.0)*(1.0+ep)-f*(ep-dl)    );
    D=(2.0+2.0*ep-f)*w_vp*w_vp;
    E=w_vp*w_vp*(2.0-f);
    F=-powf(w_vp,4);
   
    //printf("pow=%f,%f,%f,%f\n",powf(cosphi,4),A,B,C);   
 
    a4=A*powf(cosphi,4)+B*powf(sinphi,4)+C*powf(cosphi*sinphi,2);
    a3=(-A*4.0*powf(cosphi,3)*sinphi+4.0*B*cosphi*powf(sinphi,3)+2.0*C*(powf(cosphi,3)*sinphi-cosphi*powf(sinphi,3) ))*kx;
    a2=A*6.0*powf(cosphi*sinphi*kx,2)+B*( 6.0*powf(cosphi*sinphi*kx,2)+2.0*powf(sinphi*ky,2)  ) +
	C*( powf(cosphi,2)*(powf(cosphi*kx,2)+ky*ky)+ powf(sinphi,4)*kx*kx- 4.0*powf(cosphi*sinphi*kx,2))+D*powf(sinphi,2)+E*powf(cosphi,2);
    a1=A*(-4.0*cosphi*powf(sinphi*kx,3) )+B*4.0*(sinphi*powf(cosphi*kx,3)+cosphi*sinphi*kx*ky*ky) + 
	C*( 2.0*cosphi*powf(sinphi*kx,3) -2.0*cosphi*sinphi*kx*(powf(cosphi*kx,2)+ky*ky  )   ) +D*2.0*cosphi*sinphi*kx-2.0*E*cosphi*sinphi*kx;
    a0=A*powf(sinphi*kx,4)+B*( powf(cosphi*kx,4)+powf(ky,4)+2.0*powf(cosphi*kx*ky,2))+C*(powf(sinphi*kx,2)*(powf(cosphi*kx,2)+ky*ky  )   )+
	D*(powf(cosphi*kx,2)+ky*ky)+E*(powf(sinphi*kx,2)+ky*ky)+F;
    //printf("%f,%f,%f,%f,%f\n",a0,a1,a2,a3,a4);

   a4=fm1+2.0*ep*fm1*powf(sinphi,2)-f/2.0*epdl*powf(sin2phi,2);
     a3=(2.0*fm1*ep*sin2phi-f*epdl*sin4phi)*kx;
     a2=( 2*(f-1.0)*(1.0+ep)-f*epdl*(2.0*cos2phi*cos2phi-sin2phi*sin2phi ) )*kx*kx+ w_vp*w_vp*(2.0*ep*sinphi*sinphi+2.0-f)+
     2.0*ky*ky*( (f-1.0)*(1.0+ep)+(f-1.0)*ep*sinphi*sinphi-f*epdl*cosphi*cosphi);
     a1=powf(kx,3)*( 2.0*(f-1.0)*ep*sin2phi+f*epdl*sin4phi)+2.0*ep*powf(w_vp,2)*sin2phi*kx+ky*ky*2.0*sin2phi*kx*( (f-1.0)*ep+f*epdl);
     a0=powf(kx,4)*( fm1*(1.0+2.0*ep*cosphi*cosphi)-f*0.5*epdl*sin2phi*sin2phi  )+fm1*(1.0+2.0*ep)*powf(ky,4)+
     powf(kx*ky,2)*( fm1*(2.0*(1.0+ep)+2.0*ep*cosphi*cosphi)-2.0*f*epdl*sinphi*sinphi   )+
     powf(w_vp*kx,2)*(2.0-f+2.0*ep*cosphi*cosphi)-powf(w_vp,4)+powf(w_vp*ky,2)*(4.0+2.0*ep-2.0*f);

//   printf("fromaa%f,%f,%f,%f,%f\n",a0,a1,a2,a3,a4);

    a[0]=a0; a[1]=a1; a[2]=a2; a[3]=a3; a[4]=a4;

}

sf_complex kzani(float phi,float ep,float dl,float w_vp,float f,float kx)
{
    
    sf_complex x1,z[4];
    float  a[5];
    float kzsq;

    if (phi==0){
	coe(a,0.0,ep,dl,w_vp,1.0,kx);
	kzsq=-a[0]/a[2];
	if (kzsq >= 0 ) 
	    x1=sf_cmplx(sqrtf(kzsq),0.0);
	else
	    x1=sf_cmplx(0.0, -sqrtf(-kzsq));
    }
    else{
	coe(a,phi,ep,dl,w_vp,f,kx);
	//printf("2daaa=%f,%f,%f,%f,%f\n",a[4],a[3],a[2],a[1],a[0]);
	p4solver(a,z);
	x1=chooseproot(z);
    }
    return x1;
}

sf_complex kzani3d(float phi,float ep,float dl,float w_vp,float f,float kx,float ky)
/*< kz anisotropic >*/
{
    sf_complex x1,z[4];
    float  a[5];
    float kzsq;

    if (phi==0){
	coe3d(a,0.0,ep,dl,w_vp,1.0,kx,ky);
	kzsq=-a[0]/a[2];
	if (kzsq >= 0 ) 
	    x1=sf_cmplx(sqrtf(kzsq),0.0);
	else
	    x1=sf_cmplx(0.0, -sqrtf(-kzsq));
    }
    else{
	coe3d(a,phi,ep,dl,w_vp,f,kx,ky);
	//printf("3daaa=%f,%f,%f,%f,%f\n",a[0],a[1],a[2],a[3],a[4]);
	p4solver(a,z);
	//printf("3dz%f,%f,%f,%f\n",__real__ z[0],__real__ z[1],__real__ z[2],__real__ z[3]);
	//printf("3dz%f,%f,%f,%f\n",__imag__ z[0],__imag__ z[1],__imag__ z[2],__imag__ z[3]);
	//printf("2dz(%f,%f),(%f,%f),(%f,%f),(%f,%f)\n",__real__ z[0],__imag__ z[0], __real__ z[1],__imag__ z[1],__real__ z[2],__imag__ z[2],__real__ z[3],__imag__ z[3]);
	x1=chooseproot(z);
    }
    return x1;
}



sf_complex dkz(float phi,float ep,float dl,float w_vp,float f,float kx)
{
    sf_complex  kz1,dkzp;
    float kz,test;

    kz1=kzani(phi,ep,dl,w_vp,f,kx);
    //kkz1=kz1
    test=w_vp*w_vp-kx*kx;
    if (test >0.0)
	kz=sqrtf(test);
    else
	kz=0.0;

    if ( fabsf(cimagf(kz1)) <0.0000001*fabsf(crealf(kz1)) )
	dkzp=sf_cmplx(fabsf(crealf(kz1))-kz,0.0);    
    //dkzp=sf_cmplx(fabsf(crealf(kz1)),0.0);
    else
	dkzp=sf_cmplx(0.0,-fabsf(cimagf(kz1)));

    return dkzp;
}

sf_complex dkz3d(float phi,float ep,float dl,float w_vp,float f,float kx,float ky)
{
    sf_complex  kz1,dkzp;
    float kz,test;

    kz1=kzani3d(phi,ep,dl,w_vp,f,kx,ky);
    //kkz1=kz1
    test=w_vp*w_vp-(kx*kx+ky*ky);
    if (test >0.0){
	kz=sqrtf(test);
    }
    else{
	kz=0.0;
    }
    if ( fabsf(cimagf(kz1)) <0.0000001*fabsf(crealf(kz1)) )
	dkzp=sf_cmplx(fabsf(crealf(kz1))-kz,0.0);    
    //dkzp=sf_cmplx(fabsf(crealf(kz1)),0.0);
    else
	dkzp=sf_cmplx(0.0,-fabsf(cimagf(kz1)));
    return dkzp;
}


sf_complex dphase(float phi,float ep,float dl,float w_vp,float f,float kx,float dz)
{
    sf_complex dphp,dkzp;
    float pshift;

    dkzp=dkz(phi,ep,dl,w_vp,f,kx);
    if (fabsf(crealf(dkzp))*0.000001 > fabsf(cimagf(dkzp))){
	pshift=dz*crealf(dkzp);
	dphp=sf_cmplx(cosf(pshift),sinf(pshift));
    }
    else{
	pshift=dz*(-fabsf(cimagf(dkzp)));
	dphp=sf_cmplx(exp(pshift),0.0);
    }
    return dphp;
}

sf_complex dphase3d(float phi,float ep,float dl,float w_vp,float f,float kx,float ky,float dz)
{
    sf_complex dphp,dkzp;
    float pshift;

    dkzp=dkz3d(phi,ep,dl,w_vp,f,kx,ky);
    if (fabsf(crealf(dkzp))*0.000001 > fabsf(cimagf(dkzp))){
	pshift=dz*crealf(dkzp);
	dphp=sf_cmplx(cosf(pshift),sinf(pshift));
    }
    else{
	pshift=dz*(-fabsf(cimagf(dkzp)));
	dphp=sf_cmplx(exp(pshift),0.0);
    }
    return dphp;
}




float gauweit(int weitm,int im){
    float gau;
    if (im <weitm-5)
	gau=1.0;
    else{
    
	if (im<weitm && im>=weitm-5){
	    gau=(5.0-(im-weitm)-1.0)*0.1; 
	}
	else{
	    if (im <weitm+4)
		gau=(5.0-(im-weitm)-1.0)*0.02;
	    else
		gau=0.001;
	}
	gau=0.001;
    }
 
    return gau;
}

void dph(float phi,float ep,float dl,float w_vp,float f,float dkx,float dz,int m,sf_complex *dpp,int weitm)
{
    // dpp(m)
    float kx;
    sf_complex dkzp;
    float pshift;
    int im,nn;
    float cc,bb,dkzpweitm,dkzpweitm1;
    /* 
       for(im=0;im<m;im++){  //do im=1,m
       kx=im*dkx;
       dpp[im]=dphase(phi,ep,dl,w_vp,f,kx,dz);
       }
    */
  
    for (im=0;im<=weitm;im++){
	kx=im*dkx;
	dkzp=dkz(phi,ep,dl,w_vp,f,kx);
	pshift=dz*crealf(dkzp);
	dpp[im]=sf_cmplx(cosf(pshift),sinf(pshift));
    }
    dkzpweitm=crealf( dkz(phi,ep,dl,w_vp,f,dkx*weitm));
    dkzpweitm1=crealf( dkz(phi,ep,dl,w_vp,f,dkx*(weitm-1)));
    bb=(dkzpweitm- dkzpweitm1)/dkx;
    cc=fabsf(crealf(dkz(phi,ep,dl,w_vp,f,weitm*dkx)-dkz(phi,ep,dl,w_vp,f,0.0*dkx)));
    if (cc==0.0) cc=1.0;
    //printf("cc=%f,bb=%f\n",cc,bb); 

 
    for (im=weitm+1;im<m;im++){
	kx=weitm*dkx;
	kx=(im-weitm)*dkx;
	dkzp=cc*sinf(1.0/cc*bb*kx)+dkz(phi,ep,dl,w_vp,f,weitm*dkx);
	//dkzp=dkz(phi,ep,dl,w_vp,f,kx);
	pshift=dz*crealf(dkzp);
	dpp[im]=sf_cmplx(cosf(pshift),sinf(pshift));
	//if (dkx <0)
	//  printf("im=%d,(%f,%f)\n",im,__real__ dpp[im], __imag__ dpp[im]);
    }
    nn=weitm;
    nn=(m-weitm-2)/2;
    for(im=weitm;im<=weitm+2*nn+1;im++){
	if (im < m)
	    dpp[im]*=0.8+0.1*( 1.0+cosf((float)(im-weitm)*SF_PI/(float)(2*nn)  ));
    }
    //for(im=weitm+nn+nn+1;im<m;im++){
    //  dpp[im]*=0.9;
    //}
  
/*
  printf("ep=%f  dl=%f  w_vp=%f,f=%f, phi=%f\n",ep,dl,w_vp,f,phi);
  printf("dpppppppppppppp=%f\n",crealf(dpp[0]));
  printf("wwwwwwwwwwwwwwwwwwwwww=%d\n",weitm);
*/
    //weight(1:275)=1.0; weight(276:m)=0.001

    for (im=weitm+1; im <m; im++){ //do im=1,m  ?????????? wait
	//dpp[im]*=0.001;
	dpp[im]*=gauweit(weitm,im);
    }

}


void dph3d(float phi,float ep,float dl,float w_vp,float f,float dkx,float dky,float dz,int mx,int my,sf_complex *dpp,float *weit)
{
    // dpp(m)
    float kx,ky;
    float pshift;
    int imx,imy,weitmx,weitmy;
    float cc,bb=0,dkzpweitm,dkzpweitm1,nn,kx1,ky1,tmpr,scale,kx2,ky2,dkzp0;
    sf_complex phab;
    for(imy=0;imy<my;imy++){ 
	ky=imy*dky;
	for(imx=0;imx<mx;imx++){  //do im=1,m
	    kx=imx*dkx;
	    dpp[i2(imy,imx,mx)]=dphase3d(phi,ep,dl,w_vp,f,kx,ky,dz);
	}
    }
    for(imx=0;imx<mx;imx++)
	if (weit[i2(0,imx,mx)]<0.8) break;
    for(imy=0;imy<my;imy++)
	if (weit[i2(imy,0,mx)]<0.8) break;
    weitmx=imx-1; weitmy=imy-1; 
    phab=dphase3d(phi,ep,dl,w_vp,f,0.0,weitmy*dky,dz);
    //printf("%d,%dweitm",weitmx,weitmy);
    //printf("%f,%f,%f\n",weitmy*dky,__real__ phab,__imag__ phab);
    //phab=dkz3d(phi,ep,dl,w_vp,f,weitmx*dkx,0.0);
    //printf ("x(%f,%f)\n",__real__ phab,__imag__ phab);
    //phab=dkz3d(phi,ep,dl,w_vp,f,0.0,weitmy*dky);
    //printf("y(%f,%f)\n",__real__ phab,__imag__ phab);
  
    for(imy=0;imy<my;imy++){
	ky=imy*dky;
	for(imx=0;imx<mx;imx++){
	    kx=imx*dkx;
	    tmpr=(float)(imy*imy)/(float)(weitmy*weitmy)+(float)(imx*imx)/(float)(weitmx*weitmx);
	    nn=sqrtf((float)(my*my)/(float)(weitmy*weitmy))-1.0;
	    if (tmpr >=1){
		tmpr=sqrtf(tmpr);
		kx1=kx/tmpr; ky1=ky/tmpr;
		kx2=kx1*(1.0-dkx/kx1); ky2=ky1*(1.0-dkx/kx1);
		if (imy==0){
		    dkzpweitm=crealf(dkz3d(phi,ep,dl,w_vp,f,kx1,ky1));
		    dkzpweitm1=crealf(dkz3d(phi,ep,dl,w_vp,f,kx2,ky2));
		    dkzp0=crealf(dkz3d(phi,ep,dl,w_vp,f,0.0,0.0));
		    bb=(dkzpweitm-dkzpweitm1);
		}
		dkzp0=crealf(dkz3d(phi,ep,dl,w_vp,f,0.0,0.0)); 
		dkzpweitm=crealf(dkz3d(phi,ep,dl,w_vp,f,kx1,ky1));
		cc=fabsf(dkzpweitm-dkzp0);
		pshift=cc*sinf(bb/cc*weitmx*(tmpr-1))+crealf(dkz3d(phi,ep,dl,w_vp,f,kx1,ky1));
		pshift=pshift*dz;
		//pshift=-SF_PI/6.0;
		phab=sf_cmplx(cosf(pshift),sinf(pshift));
		dpp[i2(imy,imx,mx)]=phab;//dphase3d(phi,ep,dl,w_vp,f,kx1,ky1,dz)+phab;
		//if (imy==0) printf("%d,%f,%f\n",imx,bb,cc);
		//dpp[i2(imy,imx,mx)]=dphase3d(phi,ep,dl,w_vp,f,kx1,ky1,dz);
		scale=0.2+0.4*( 1.0+cosf(SF_PI*(tmpr-1.0)/nn));
		if ((tmpr-1.0)>nn ) scale=0.2;
		dpp[i2(imy,imx,mx)]*=scale;
	    }
	}
    }

    for(imy=0;imy<my;imy++)
	for(imx=0;imx<mx;imx++)
	    dpp[i2(imy,imx,mx)]*=weit[i2(imy,imx,mx)]; 
   

}




int weit(float phi,float ep,float dl,float w_vp,float f,float dkx,int m)
{
    int weitm=0,im;
    sf_complex  kz1,kzp1;
    float kx,weight;

    for (im=0; im <m; im++){
	kx=im*dkx;
	kz1=kzani(phi,ep,dl,w_vp,f,-kx);
	kzp1=kzani(phi,ep,dl,w_vp,f,kx);
	if ( fabsf(cimagf(kz1)) <0.0000001*fabsf(crealf(kz1)) &&  fabsf(cimagf(kzp1)) <0.0000001*fabsf(crealf(kzp1))){
	    if (crealf(kz1) >w_vp*0.0 && crealf(kzp1) >w_vp*0.0  )
		weight=1.00;
	    else
		weight=0.001;
	}
	else
	    weight=0.001;
	if (weight!=1.0){
	    weitm=im-2;
	    break;
	}
    }
    if (im==m) weitm=m-2;
    if (phi==0.0) weitm-=2;
    if ( phi<=15.0/180.0*SF_PI)
	weitm=weitm*(1.0-1.0/4.0);
    else
	weitm=weitm*(1-1.0/6.0);
  
    //printf("phi=%f,weitm=%d,m=%d\n",phi,weitm,m);
    if (phi<=15.0/180.0*SF_PI && weitm >1.5/3.0*m) weitm=1.5/3.0*m;
    return weitm;
}

void weit3d(float *weight,float phi,float ep,float dl,float w_vp,float f,float dkx,float dky,int mx,int my)
{
    int i_mx,i_my,im,ii_mx,ii_my;
    sf_complex  kz1,kzp1;
    float kx,ky;
    float tmpr;

    for (i_my=0;i_my<my;i_my++){
	ky=i_my*dky;
	for(i_mx=0;i_mx<mx;i_mx++){
	    kx=i_mx*dkx;
	    im=i_my*mx+i_mx;
	    kz1= kzani3d(phi,ep,dl,w_vp,f,-kx,ky);
	    kzp1=kzani3d(phi,ep,dl,w_vp,f,kx,ky);
	    if ( fabsf(cimagf(kz1)) <0.0000001*fabsf(crealf(kz1)) &&  fabsf(cimagf(kzp1)) <0.0000001*fabsf(crealf(kzp1))){
		if (crealf(kz1) >w_vp*0.0 && crealf(kzp1) >w_vp*0.0  )
		    weight[im]=1.00;
		else
		    weight[im]=0.001;
	    }
	    else
		weight[im]=0.001; 
	}
    }
    for (i_my=0;i_my<my;i_my++){
	for(i_mx=0;i_mx<mx;i_mx++){
	    im=i_my*mx+i_mx;
	    if( 1.0!=weight[im]) break;
	}
	if ( i_mx <2){
	    for (ii_mx=0;ii_mx<mx;ii_mx++)
		weight[i2(i_my,ii_mx,mx)]=0.001;
	    break;
	}
	else{
	    for (ii_mx=i_mx-2;ii_mx<mx;ii_mx++){
		weight[i2(i_my,ii_mx,mx)]=0.001;
	    }
	}
    }
    printf("%d dddd\n",i_my);
    if(i_my>=2)
	for(ii_my=i_my;ii_my<my;ii_my++)
	    for(ii_mx=0;ii_mx<mx;ii_mx++)
		weight[i2(ii_my,ii_mx,mx)]=0.001;
    for(i_mx=0;i_mx<mx;i_mx++){
	if (1.0!=weight[i_mx]) break;
    }
    printf("i_mx1=%d,i_my1=%d\n",i_mx,i_my);
    i_mx=i_mx*5.0/6.0; i_my=i_my*4.5/6.0;
    //i_my=i_mx;
    if (i_mx <2) i_mx=2; 
    if(i_my<2) i_my=2;
    for (ii_my=0;ii_my<my;ii_my++){
	ky=ii_my*dky;
	for(ii_mx=0;ii_mx<mx;ii_mx++){
	    im=ii_my*mx+ii_mx;
	    kx=ii_mx*dkx;
	    tmpr=(float)(ii_my*ii_my)/(float)(i_my*i_my)+(float)(ii_mx*ii_mx)/(float)(i_mx*i_mx);
	    if (tmpr<1) weight[im]=1.0;
	    else weight[im]=0.001;
	}
    }
    printf("i_mx=%d,i_my=%d\n",i_mx,i_my);
}

void phase_apro_matrix_cos(float *a,float dx,float dkx,int m,int n)
{
    // a(m,n)
    int i_m, i_n;
    float x,kx;

    for(i_m=0; i_m<m; i_m++){
	kx=i_m*dkx;
	for(i_n=0;i_n<n; i_n++){ 
	    x=i_n*dx;
	    a[i2(i_m,i_n,n)]=cosf(x*kx);
	}
    }
} 

void phase_apro_matrix_cos3d(float *a,float dx,float dy,float dkx,float dky,int mx,int my,int nx,int ny)
{
    // a(m,n)
    int i_mx,i_my, i_nx,i_ny,i_m,i_n,n;
    float x,kx,y,ky;
    n=nx*ny;
    for(i_my=0;i_my<my;i_my++){
	ky=i_my*dky;
	for(i_mx=0;i_mx<mx;i_mx++){
	    kx=i_mx*dkx;
	    i_m=i_my*mx+i_mx;
	    for(i_ny=0;i_ny<ny;i_ny++){
		y=i_ny*dy;
		for(i_nx=0;i_nx<nx;i_nx++){
		    x=i_nx*dx;
		    i_n=i_ny*nx+i_nx;
		    a[i2(i_m,i_n,n)]=cosf(y*ky)*cosf(x*kx);
		}
	    }
	}
    }
} 


void phase_apro_matrix_sin(float *a,float dx,float dkx,int m,int n)
{  // a(m-1,n-1)

    int i_m,i_n;
    float x,kx;

    for (i_m=0; i_m < m-1; i_m++){
	kx=(i_m+1)*dkx;
	for(i_n=0; i_n <n-1; i_n++){
	    x=(i_n+1)*dx;
	    a[i2(i_m,i_n,n-1)]=sinf(x*kx); 
	}
    }
}


void phase_apro_matrix_sin3d(float *a,float dx,float dy,float dkx,float dky,int mx,int my,int nx,int ny)
{  // a(m-1,n-1)
    int i_m,i_n,i_mx,i_my,i_nx,i_ny,n;
    float x,kx,y,ky;
    n=(nx-1)*ny;
    for(i_my=0;i_my<my;i_my++){
	ky=i_my*dky; 
	for(i_mx=0;i_mx<mx-1;i_mx++){
	    kx=(i_mx+1)*dkx;
	    i_m=i_my*(mx-1)+i_mx;
	    for(i_ny=0;i_ny<ny;i_ny++){
		y=i_ny*dy;
		for(i_nx=0;i_nx<nx-1;i_nx++){
		    x=(i_nx+1)*dx;
		    i_n=i_ny*(nx-1)+i_nx;
		    a[i2(i_m,i_n,n)]=cosf(y*ky)*sinf(x*kx);

		}
	    } 
	}
    }


}

void convlv_coe(int m,int n,sf_complex *conap,sf_complex *conam,sf_complex *dpc,sf_complex *dps, float *qc,float *rc,float *qs,float *rs){ // conap(n) conam(n) dpc(m) dps(m) qc(m,n) rc(n,n) qs(m-1,n-1) rs(n-1,n-1)
    sf_complex dpp1,dpp2;
    sf_complex *conas;
    int i;
    conas=(sf_complex *)malloc((n-1)*sizeof(sf_complex));
    for(i=0;i<m;i++){
	dpp1=(dpc[i]+dps[i])/2.0;
	dpp2=(dpc[i]-dps[i])/2.0;
	dpc[i]=dpp1;
	if (i > 0 ) dps[i-1]=dpp2;
    }

//for (i=0;i<m;i++)
//   printf("1dps[%d]=(%f,%f)\n",i+1,__real__ dps[i],__imag__ dps[i]);
//for(i=0;i<n;i++)
//    printf("rs(6,%d)=%f\n",i+1,rs[i2(5,i,n-1)]);

    qrcon(qc,rc,dpc,conap,m,n);
    qrcon(qs,rs,dps,conas,m-1,n-1);
    for(i=1;i<n;i++){
	conam[i]=0.5*conap[i]-0.5*(sf_cmplx(0.0,1.0))*conas[i-1];
	conap[i]=0.5*conap[i]+0.5*(sf_cmplx(0.0,1.0))*conas[i-1];
    }
//for (i=0;i<n;i++)
//   printf("1conap[%d]=(%f,%f)\n",i+1,__real__ conap[i],__imag__ conap[i]);
    free(conas);
}


void convlv_coe3d(int mx,int my,int nx,int ny,sf_complex *conap,sf_complex *conam,sf_complex *dpc,sf_complex *dps, float *qc,float *rc,float *qs,float *rs){ // conap(n) conam(n) dpc(m) dps(m) qc(m,n) rc(n,n) qs(m-1,n-1) rs(n-1,n-1)
    sf_complex dpp1,dpp2;
    sf_complex *conas;
    int im,m,n,mm1,nm1,imx,imy,imm,inx,iny,inm,in;
    m=mx*my; n=nx*ny; mm1=my*(mx-1); nm1=ny*(nx-1);
//conas=(sf_complex *)malloc((nm1)*sizeof(sf_complex));
    conas=sf_complexalloc(nm1);
    for(imy=0;imy<my;imy++){
	for(imx=0;imx<mx;imx++){
	    im=imy*mx+imx;
	    dpp1=(dpc[im]+dps[im])/2.0;
	    dpp2=(dpc[im]-dps[im])/2.0;
	    dpc[im]=dpp1;
	    if (imx>0){
		imm=imy*(mx-1)+imx-1;
		dps[imm]=dpp2;
	    }  
	}
    }
    qrcon(qc,rc,dpc,conap,m,n);
    qrcon(qs,rs,dps,conas,mm1,nm1);
//for(i=0;i<n;i++)
//      printf("i=%d,(%f,%f)   ",i,__real__ conas[i],__imag__ conas[i]);
//for(i=0;i<n;i++)
//        printf("i=%d,%f   ",i,rs[i]);
//printf("\n");
    for(iny=0;iny<ny;iny++){
	for(inx=1;inx<nx;inx++){
	    in=iny*nx+inx;
	    inm=iny*(nx-1)+inx-1;
	    conam[in]=0.5*conap[in]-0.5*(sf_cmplx(0.0,1.0))*conas[inm];
	    conap[in]=0.5*conap[in]+0.5*(sf_cmplx(0.0,1.0))*conas[inm];
	}
    }

    free(conas);
}





void ani_equation_con(float dx,float dkx,int m,int n,int weitm,float *qc,float *rc,float *qs,float *rs)
{  // qc(m,n),rc(n,n),qs(m-1,n-1),rs(n-1,n-1)
//real :: weight(:)
    float *ma,*q,*tmpd;
    int  sing=0,im,in;

    ma=sf_floatalloc(m*n);
    q=sf_floatalloc(m*m);
    tmpd=sf_floatalloc(n);
    //allocate(ma(m,n),q(m,m),tmpd(n))

    phase_apro_matrix_cos(ma,dx,dkx,m,n);
    for(im=weitm+1;im<m; im++)
	for(in=0;in<n;in++)
	    //ma[i2(im,in,n)]*=0.001;
	    ma[i2(im,in,n)]*=gauweit(weitm,im);
    //do im=1,m
    //   ma(im,:)=ma(im,:)*weight(im)
    //end do

    qrdcmp(ma,m,n,tmpd,sing);
    qrdpost(ma,m,n,tmpd,q,rc);
    for(im=0;im<m;im++)
	for(in=0;in<n;in++)
	    qc[i2(im,in,n)]=q[i2(im,in,m)];//qc=q(1:m,1:n)

    phase_apro_matrix_sin(ma,dx,dkx,m,n);
    for(im=weitm;im<m-1; im++)
	for(in=0;in<n-1;in++)
	    //ma[i2(im,in,n-1)]*=0.001; 
	    ma[i2(im,in,n-1)]*=gauweit(weitm,im);
    //do im=1,m-1
    //  ma(im,:)=ma(im,:)*weight(im+1) //! ????????????????????????????????
    //end do
    qrdcmp(ma,m-1,n-1,tmpd,sing);
    qrdpost(ma,m-1,n-1,tmpd,q,rs);
    for(im=0;im<m-1;im++)
	for(in=0;in<n-1;in++)
	    qs[i2(im,in,n-1)]=q[i2(im,in,m-1)];  //qs=q(1:m-1,1:n-1)
    free(ma);free(q);free(tmpd);
}



void ani_equation_con3d(float dx,float dy,float dkx,float dky,int mx,int my,int nx,int ny,float *weight,float *qc,float *rc,float *qs,float *rs)
{  // qc(m,n),rc(n,n),qs(m-1,n-1),rs(n-1,n-1)
//real :: weight(:)
    float *ma,*q,*tmpd;
    int  sing=0,im,in,m,n,mm1,nm1,imx,imy,iwm;
    m=mx*my; n=nx*ny;
    ma=sf_floatalloc(m*n);
    q=sf_floatalloc(m*m);
    tmpd=sf_floatalloc(n);
    //allocate(ma(m,n),q(m,m),tmpd(n))
    printf("1\n");
    phase_apro_matrix_cos3d(ma,dx,dy,dkx,dky,mx,my,nx,ny);
    for(im=0;im<m; im++)
	for(in=0;in<n;in++)
	    ma[i2(im,in,n)]*=weight[im];
    printf("2\n");
    qrdcmp(ma,m,n,tmpd,sing);
    printf("3\n");
    qrdpost(ma,m,n,tmpd,q,rc);
    printf("5\n");
    for(im=0;im<m;im++)
	for(in=0;in<n;in++)
	    qc[i2(im,in,n)]=q[i2(im,in,m)];//qc=q(1:m,1:n)
    phase_apro_matrix_sin3d(ma,dx,dy,dkx,dky,mx,my,nx,ny);
    printf("4\n");
    mm1=my*(mx-1); nm1=ny*(nx-1);
    for(imy=0;imy<my;imy++){
	for(imx=0;imx<mx-1;imx++){
	    im=imy*(mx-1)+imx;
	    iwm=imy*mx+imx+1;
	    for(in=0;in<nm1;in++)
		ma[i2(im,in,nm1)]*=weight[iwm];
	}
    }
    //do im=1,m-1
    //  ma(im,:)=ma(im,:)*weight(im+1) //! ????????????????????????????????
    //end do
    qrdcmp(ma,mm1,nm1,tmpd,sing);
    /*
      printf("\n****************\n");
      for(in=0;in<n;in++)
      printf("***%f",tmpd[in]);
      printf("\n*************\n");
    */
    qrdpost(ma,mm1,nm1,tmpd,q,rs);
    for(im=0;im<mm1;im++)
	for(in=0;in<nm1;in++)
	    qs[i2(im,in,nm1)]=q[i2(im,in,mm1)];  //qs=q(1:m-1,1:n-1)
    free(ma);free(q);free(tmpd);
}





void ani_init(float dx,float dkx,int m,int n,int nweit,int *weitm,float *qc_table, float *rc_table, float *qs_table,float *rs_table)
{  // weitm(nweit),qc_table(nweit,m,n),rc_table(nweit,n,n),qs_table(nweit,m-1,n-1),rs_table(nweit,n-1,n-1)
    int im,in,iweit;
    int nqc[2],nrc[2],nqs[2],nrs[2];
    float *qc,*rc,*qs,*rs;
    d3(m,n,nqc); d3(n,n,nrc); d3(m-1,n-1,nqs); d3(n-1,n-1,nrs); 
    qc=sf_floatalloc(m*n); rc=sf_floatalloc(n*n); qs=sf_floatalloc((m-1)*(n-1)); rs=sf_floatalloc((n-1)*(n-1));
    // qc(m,n),rc(n,n),qs(m-1,n-1),rs(n-1,n-1)

    for(iweit=0;iweit<nweit; iweit++){
	//printf("iweit=%d  nweit=%d   \n",iweit,nweit);
	ani_equation_con(dx,dkx,m,n,weitm[iweit],qc,rc,qs,rs);
	for(im=0;im<m;im++)
	    for(in=0;in<n;in++)
		qc_table[i3(iweit,im,in,nqc)]=qc[i2(im,in,n)]; //qc_table(iweit,:,:)=qc(:,:)
	for(im=0;im<n;im++)
	    for(in=0;in<n;in++)
		rc_table[i3(iweit,im,in,nrc)]=rc[i2(im,in,n)];  //rc_table(iweit,:,:)=rc(:,:)
	for(im=0;im<m-1;im++)
	    for(in=0;in<n-1;in++)
		qs_table[i3(iweit,im,in,nqs)]=qs[i2(im,in,n-1)];    //qs_table(iweit,:,:)=qs(:,:)
	for(im=0;im<n-1;im++)
	    for(in=0;in<n-1;in++)
		rs_table[i3(iweit,im,in,nrs)]=rs[i2(im,in,n-1)];     //rs_table(iweit,:,:)=rs(:,:)
    }
    free(qc);free(rc);free(qs);free(rs);
}

void explicit_table3d(sf_complex *conapp,sf_complex *conapm,
                      float phi,float f,float dx,float dy,float dkx,float dky,float dz,int mx,int my,int nx,int ny,
                      float ep,float dl,float w_vp,float *weight,float *qc,float *rc,float *qs,float *rs)
{
    int m;
    //float *qc,*rc,*qs,*rs,*weight;
    sf_complex *dppc,*dpps;
    m=my*mx;
    //qc=sf_floatalloc(m*n); rc=sf_floatalloc(n*n); qs=sf_floatalloc(mm1*nm1); rs=sf_floatalloc(nm1*nm1); weight=sf_floatalloc(m);
    dppc=sf_complexalloc(m); dpps=sf_complexalloc(m);
    printf("bbbbbbbbbbbbb\n");
    //weit3d(weight,phi,ep,dl,w_vp,f,dkx,dky,mx,my);
    dph3d(phi,ep,dl,w_vp,f,dkx,dky,dz,mx,my,dppc,weight);
    dph3d(phi,ep,dl,w_vp,f,-dkx,dky,dz,mx,my,dpps,weight);
    printf("aaaaaaaaaaa\n");
    //ani_equation_con3d(dx,dy,dkx,dky,mx,my,nx,ny,weight,qc,rc,qs,rs);
    printf("cccccccccccccc\n");
    convlv_coe3d(mx,my,nx,ny,conapp,conapm,dppc,dpps, qc,rc,qs,rs);
    printf("dddddddddddd \n");  

    //free(qc); free(rc); free(qs); free(rs); free(weight);
    free(dppc); free(dpps); 
  
}


void explicit_table2d(sf_complex *conapp,sf_complex *conapm,
                      float phi,float f,float dx,float dkx,float dz,int mx,int nx,
                      float ep,float dl,float w_vp)
{
    int m,n,mm1,nm1;
    float *qc,*rc,*qs,*rs,*weight;
    int my,ny; 
    float dy,dky;
    my=1; ny=1; dy=dx; dky=dkx;
    sf_complex *dppc,*dpps;
    m=my*mx; n=ny*nx; mm1=my*(mx-1); nm1=ny*(nx-1);
    qc=sf_floatalloc(m*n); rc=sf_floatalloc(n*n); qs=sf_floatalloc(mm1*nm1); rs=sf_floatalloc(nm1*nm1); weight=sf_floatalloc(m);
    dppc=sf_complexalloc(m); dpps=sf_complexalloc(m);
    printf("bbbbbbbbbbbbb\n");
    weit3d(weight,phi,ep,dl,w_vp,f,dkx,dky,mx,my);
    dph3d(phi,ep,dl,w_vp,f,dkx,dky,dz,mx,my,dppc,weight);
    dph3d(phi,ep,dl,w_vp,f,-dkx,dky,dz,mx,my,dpps,weight);
  
    printf("aaaaaaaaaaa\n");
    ani_equation_con3d(dx,dy,dkx,dky,mx,my,nx,ny,weight,qc,rc,qs,rs);
    printf("cccccccccccccc\n");
    convlv_coe3d(mx,my,nx,ny,conapp,conapm,dppc,dpps, qc,rc,qs,rs);
    printf("dddddddddddd \n");  
    free(qc); free(rc); free(qs); free(rs); free(dppc); free(dpps); free(weight);
  
}




void explicit_table(sf_complex *contablepp,sf_complex *contablepm,
                    float phi,float f,float dx,float dkx,float dz, int m,int n,
                    float o_ep,float d_ep,int n_ep,float o_dl,float d_dl,int n_dl,float o_w_vp,float d_w_vp,int n_w_vp)
/*< explicit table >*/
{  //contablepp(n_w_vp,n_ep,n_dl,n), contablepm(n_w_vp,n_ep,n_dl,n)
    sf_complex *conapp,*conapm,*dppc,*dpps;
    float *qc_table,*rc_table,*qs_table,*rs_table;
    float *qc,*rc,*qs,*rs;
    float w_vp,ep,dl;
    int *weitm,nweit,theweitm,iweit;
    int im,in,i_w_vp,i_ep,i_dl;
    int nqc[2],nrc[2],nqs[2],nrs[2],nconp[3];
    printf("begin 11111\n");
    d3(m,n,nqc); d3(n,n,nrc); d3(m-1,n-1,nqs); d3(n-1,n-1,nrs); d4(n_ep,n_dl,n,nconp);
    printf("begin 22222\n");
    weitm=sf_intalloc(m);
    printf("begin 33333\n");
    nweit=weitmlist(weitm,phi,dkx,f,m,o_ep,d_ep,n_ep,o_dl,d_dl,n_dl,o_w_vp,d_w_vp,n_w_vp);
    printf("begin 44444\n");
    printf("mmmm=%d,nweit=%d\n",m,nweit);
    conapp=sf_complexalloc(n); conapm=sf_complexalloc(n); dppc=sf_complexalloc(m); dpps=sf_complexalloc(m);
    qc_table=sf_floatalloc(nweit*m*n); 
    rc_table=sf_floatalloc(nweit*n*n); 
    qs_table=sf_floatalloc(nweit*(m-1)*(n-1) ); rs_table=sf_floatalloc(nweit*(n-1)*(n-1));
    qc=sf_floatalloc(m*n); rc=sf_floatalloc(n*n); qs=sf_floatalloc((m-1)*(n-1)); rs=sf_floatalloc((n-1)*(n-1));

//write(*,*) "begin ani_init"
    ani_init(dx,dkx,m,n,nweit,weitm,qc_table,rc_table,qs_table,rs_table);
//write(*,*) "end ani_init"
    printf("aaaaaaaaaaaa1\n");
//n_w_vp=1; n_ep=1; n_dl=1;
    for(i_w_vp=0;i_w_vp<n_w_vp; i_w_vp++){  //do i_w_vp=1,n_w_vp
	w_vp=o_w_vp+i_w_vp*d_w_vp;
	for(i_ep=0;i_ep<n_ep;i_ep++){//do i_ep=1,n_ep
	    ep=i_ep*d_ep+o_ep;
	    for(i_dl=0;i_dl<n_dl;i_dl++){ //do i_dl=1,n_dl
		dl=i_dl*d_dl+o_dl;
       
		theweitm=weit(phi,ep,dl,w_vp,f,dkx,m); 
		if (0>=theweitm) theweitm=1;
		dph(phi,ep,dl,w_vp,f,dkx,dz,m,dppc,theweitm);  //call dph(phi,ep,dl,w_vp,f,0.0,dkx,dz,m,dppc,dpmc)
		//if (i_ep==0 && i_dl==0 && i_w_vp==0 ){
		//   for(i=0;i<m;i++)
		//       printf("con[%d]=%f \n",i,__real__ dppc[i] );
		//}
       
     
		dph(phi,ep,dl,w_vp,f,-dkx,dz,m,dpps,theweitm); //call dph(phi,ep,dl,w_vp,f,0.0,-dkx,dz,m,dpps,dpms) 
		//if (i_ep==0 && i_dl==0 && i_w_vp==0 ){
		//   for(i=0;i<m;i++)
		//       printf("con[%d]=%f \n",i,__real__ dpps[i] );
		//}

		//!      call ani_equation_con(dx,dkx,m,n,weightp)
   
        

		for(iweit=0;iweit<nweit;iweit++){
		    if(theweitm==weitm[iweit]) break;
		}
      
		//printf("theweitm=%d,1\n",theweitm); 
		for(im=0;im<m;im++)
		    for(in=0;in<n;in++)
			qc[i2(im,in,n)]=qc_table[i3(iweit,im,in,nqc)];     //qc(:,:)=qc_table(iweit,:,:);
		//printf("2\n"); 
		for(im=0;im<n;im++)
		    for(in=0;in<n;in++)
			rc[i2(im,in,n)]=rc_table[i3(iweit,im,in,nrc)];    //  rc(:,:)=rc_table(iweit,:,:);
		//printf("3,iweit=%d,nweit=%d\n",iweit,nweit); 
		for(im=0;im<m-1;im++)
		    for(in=0;in<n-1;in++)
			qs[i2(im,in,n-1)]=qs_table[i3(iweit,im,in,nqs)];     //qs(:,:)=qs_table(iweit,:,:);
		//printf("4\n"); 
		for(im=0;im<n-1;im++)
		    for(in=0;in<n-1;in++)
			rs[i2(im,in,n-1)]=rs_table[i3(iweit,im,in,nrs)];   //rs(:,:)=rs_table(iweit,:,:);
      
		convlv_coe(m,n,conapp,conapm,dppc,dpps,qc,rc,qs,rs);
		//printf("i_dl=%d\n",i_dl); 
		//if (i_ep==0 && i_dl==0 && i_w_vp==0 ){
		//   for(i=0;i<n;i++)
		//     printf("con[%d]=(%f,%f) ",i,__real__ conapp[i],__imag__ conapp[i] ); 
		//   printf("\n");
		//}
      
		for(in=0;in<n;in++){
		    contablepp[i4(i_w_vp,i_ep,i_dl,in,nconp)]=conapp[in];  //contablepp(:,i_dl,i_ep,i_w_vp)=conapp(:)
		    contablepm[i4(i_w_vp,i_ep,i_dl,in,nconp)]=conapm[in];  //contablepm(:,i_dl,i_ep,i_w_vp)=conapm(:)
		}
      
		if (d_ep>fabsf(ep) && d_dl>fabsf(dl)){
		    for(in=0;in<n;in++){
			contablepp[i4(i_w_vp,i_ep,i_dl,in,nconp)]=sf_cmplx(0.0,0.0);
			contablepm[i4(i_w_vp,i_ep,i_dl,in,nconp)]=sf_cmplx(0.0,0.0);
		    }
		    contablepp[i4(i_w_vp,i_ep,i_dl,0,nconp)]=sf_cmplx(1.0,0.0);
		}
	    }
	}
    }
//write(*,*) "end concoe"
    free(conapp); free(conapm); 
    free(dppc); free(dpps);
    free(qc_table); free(rc_table); 
    free(qs_table); free(rs_table);
    free(qc);free(rc);free(qs);free(rs);
    free(weitm);
}

int weitmlist(int *weitm,float phi,float dkx,float f,int m,
	      float o_ep,float d_ep,int n_ep,
	      float o_dl,float d_dl,int n_dl,
	      float o_w_vp,float d_w_vp,int n_w_vp)
/*< WEI >*/
{
    int nweit,iweit;
    /*
      weitm_tmp=sf_intalloc(m);  
      printf("weit1\n");
      nweit=0;
      for(i_w_vp=0;i_w_vp<n_w_vp; i_w_vp++){
      w_vp=o_w_vp+i_w_vp*d_w_vp;
      for(i_ep=0;i_ep<n_ep; i_ep++){
      ep=o_ep+i_ep*d_ep;
      for(i_dl=0;i_dl<n_dl; i_dl++){
      dl=o_dl+d_dl*i_dl;
      theweitm=weit(phi,ep,dl,w_vp,f,dkx,m);
      for(iweit=0;iweit<nweit; iweit++){
      if (theweitm==weitm_tmp[iweit]) break;
      }
      if (iweit==nweit){
      weitm_tmp[nweit]=theweitm;
      nweit++;
      }
      }
      }
      }
      printf("weit2\n");
      for(iweit=0; iweit<nweit; iweit++)
      weitm[iweit]=weitm_tmp[iweit];
      free(weitm_tmp);
    */
    nweit=m*(1-1.0/6.0);
    for(iweit=0;iweit<nweit;iweit++)
	weitm[iweit]=iweit+1;
    return nweit;
}


float kycoe3d(float phi,float ep,float dl,float w_vp,float f, float kx,float kz){
    float sinphi,cosphi,fm1,epdl,sin2phi,cos2phi,sin4phi;
    float A,B,C,D,E,F;
    float b0,b2,b4;
    float b24ac,ab,ky=0,ky2;
    sinphi=-sinf(phi); cosphi=cosf(phi);
    fm1=f-1; epdl=ep-dl; sin2phi=-sinf(2.0*phi); cos2phi=cosf(2.0*phi); sin4phi=-sinf(4.0*phi);
    //printf("sin=%f\n",sinphi);
    A=f-1;
    B=(f-1.0)*(1.0+2.0*ep);
    C=2.0*( (f-1.0)*(1.0+ep)-f*(ep-dl)    );
    D=(2.0+2.0*ep-f)*w_vp*w_vp;
    E=w_vp*w_vp*(2.0-f);
    F=-powf(w_vp,4);
   
    //printf("pow=%f,%f,%f,%f\n",powf(cosphi,4),A,B,C);   
    b4=B; b2=(B*(2.0*powf(sinphi,2))+C*powf(cosphi,2))*powf(kz,2)+(B*4.0*(cosphi*sinphi*kx)-C*2.0*cosphi*sinphi*kx)*kz+
	      B*2.0*powf(cosphi*kx,2)+C*powf(sinphi*kx,2)+D+E;
    b0=(A*powf(cosphi,4)+B*powf(sinphi,4)+C*powf(cosphi*sinphi,2))*powf(kz,4)+((-A*4.0*powf(cosphi,3)*sinphi+4.0*B*cosphi*powf(sinphi,3)+2.0*C*(powf(cosphi,3)*sinphi-cosphi*powf(sinphi,3) ))*kx)*powf(kz,3)+(A*6.0*powf(cosphi*sinphi*kx,2)+
																								      B*( 6.0*powf(cosphi*sinphi*kx,2))+ 
																								      C*( powf(cosphi,2)*powf(cosphi*kx,2)+powf(sinphi,4)*kx*kx- 4.0*powf(cosphi*sinphi*kx,2))+D*powf(sinphi,2)+E*powf(cosphi,2))*powf(kz,2)+(A*(-4.0*cosphi*powf(sinphi*kx,3) )+B*4.0*(sinphi*powf(cosphi*kx,3))+
																																								       D*2.0*cosphi*sinphi*kx-2.0*E*cosphi*sinphi*kx+
																																								       C*( 2.0*cosphi*powf(sinphi*kx,3) -2.0*cosphi*sinphi*kx*(powf(cosphi*kx,2))))*kz+A*powf(sinphi*kx,4)+B*powf(cosphi*kx,4)+C*powf(sinphi*kx,2)*powf(cosphi*kx,2)+D*powf(cosphi*kx,2)+E*powf(sinphi*kx,2)+F;


    //printf("b0=%f,b2=%f,b4=%f\n",b0,b2,b4);

    b0=(fm1+2.0*ep*fm1*powf(sinphi,2)-f/2.0*epdl*powf(sin2phi,2))*powf(kz,4)+
	((2.0*fm1*ep*sin2phi-f*epdl*sin4phi)*kx)*powf(kz,3)+
	(( 2*(f-1.0)*(1.0+ep)-f*epdl*(2.0*cos2phi*cos2phi-sin2phi*sin2phi ) )*kx*kx+w_vp*w_vp*(2.0*ep*sinphi*sinphi+2.0-f))*powf(kz,2)+
	(powf(kx,3)*(2.0*(f-1.0)*ep*sin2phi+f*epdl*sin4phi)+2.0*ep*powf(w_vp,2)*sin2phi*kx)*kz+
	powf(kx,4)*( fm1*(1.0+2.0*ep*cosphi*cosphi)-f*0.5*epdl*sin2phi*sin2phi  )+powf(w_vp*kx,2)*(2.0-f+2.0*ep*cosphi*cosphi)-powf(w_vp,4);

    b4=fm1*(1.0+2.0*ep);
    b2=powf(kx,2)*( fm1*(2.0*(1.0+ep)+2.0*ep*cosphi*cosphi)-2.0*f*epdl*sinphi*sinphi)+powf(w_vp,2)*(4.0+2.0*ep-2.0*f)+
	2.0*sin2phi*kx*( (f-1.0)*ep+f*epdl)*kz+2.0*( (f-1.0)*(1.0+ep)+(f-1.0)*ep*sinphi*sinphi-f*epdl*cosphi*cosphi)*powf(kz,2);
    //printf("b0=%f,b2=%f,b4=%f,b^2-4ac=%f\n",b0,b2,b4,b2*b2-4*b4*b0);
    b24ac=b2*b2-4.0*b4*b0;
    ab=b4*b2;
    if (b24ac*10E9>0.0){
	if (ab*10E9<0.0){
          
	    ky=sqrtf((-b2-sqrtf(b24ac))/(2.0*b4));
	    ky2=sqrtf((-b2+sqrtf(b24ac))/(2.0*b4));
	    //printf("ky=%f,ky2=%f,%f,%f,%f\n",ky,ky2,b24ac,(-b2+sqrtf(b24ac))/(2.0*b4),ab);
	    if ((-b2+sqrtf(b24ac))/(2.0*b4)>0) ky=ky2;
	    else ky=0.0;
	    //if (ky2>0) ky=sqrtf(ky2);
	    //else ky=0.0;
	    //if (ky>ky2) ky=ky2;
	}
    }
    else ky=0.0;
    return ky;
}



//--------------------implicit method------------------------------------------------------

void kzvti(sf_complex *sz,float dsx,float ep,float dl,int nsx){

    float sx;
    int isx;
    for (isx=0;isx<nsx;isx++){
	sx=isx*dsx;
	sz[isx]=sqrtf((1-sx*sx*(1+2*ep))/(1-sx*sx*2.0*(ep-dl)))-1.0;

    }
}

void implicitleftmatrixisovtifd2(float **leftmatrix, sf_complex * sz,int nsx,float dsx)
{
    int isx;
    float sx,sx2,sx4;
    for (isx=0;isx<nsx;isx++){
	sx=dsx*isx;
	sx2=sx*sx; sx4=sx2*sx2;
	leftmatrix[isx][0]=sz[isx]*sx2;
	leftmatrix[isx][1]=sz[isx]*sx4;
	leftmatrix[isx][2]=-sx2;
	leftmatrix[isx][3]=-sx4;
    }
}


void implicitleftmatrixisovtifd1(float **leftmatrix, sf_complex * sz,int nsx,float dsx)
{
    int isx;
    float sx,sx2;
    for (isx=0;isx<nsx;isx++){
	sx=dsx*isx;
	sx2=sx*sx; 
	leftmatrix[isx][0]=sz[isx]*sx2;
	leftmatrix[isx][1]=-sx2;
    }
}


void qrsolver(float **leftmatrix,sf_complex *kz,sf_complex *conap,int m,int n)
{
    float *q,*qc,*rc,*tmpd;
    int sing=0;
    int im,in; 
    q=sf_floatalloc(m*m); qc=sf_floatalloc(m*n); rc=sf_floatalloc(n*n); tmpd=sf_floatalloc(n);
  
    qrdcmp(*leftmatrix,m,n,tmpd,sing);
    qrdpost(*leftmatrix,m,n,tmpd,q,rc);
    for(im=0;im<m;im++)
	for(in=0;in<n;in++)
	    qc[i2(im,in,n)]=q[i2(im,in,m)];
    for (im=0;im<m;im++) kz[im]=kz[im]*(-1);
    qrcon(qc,rc,kz,conap,m,n);
    free(q); free(qc); free(rc); free(tmpd);
}

void imfd1coesolver(float ep,float dl,int nsxall, float *alpha1,float *beta1)
{
    float **leftmatrixfd1;
    sf_complex *sz,*conapfd1;
    int nsx,nfd1=2;
    int isx,in;
    float dsx,sx,weight,maxsx;
    maxsx=sqrtf(1.0/(1.0+2.0*ep));
    dsx=maxsx/(nsxall-1);
    sz=sf_complexalloc(nsxall);  conapfd1=sf_complexalloc(nfd1);
    kzvti(sz,dsx,ep,dl,nsxall);
    for(isx=nsxall-1;isx>=0; isx--){
	sx=isx*dsx;
	if (sqrtf(sx*sx/(sx*sx+(sz[isx]+1.0)*(sz[isx]+1.0)))<=sinf(65.0/180.0*SF_PI)) break;
    }
    nsx=isx;
    leftmatrixfd1=sf_floatalloc2(nsx,nfd1);
    implicitleftmatrixisovtifd1(leftmatrixfd1,sz,nsx,dsx);
    for(isx=0;isx<nsx;isx++){
	weight=sqrtf(sx*sx+(sz[isx]+1.0)*(sz[isx]+1.0))/(sz[isx]+1.0);
	sz[isx]=sz[isx]*weight;
	for(in=0;in<nfd1;in++){
	    leftmatrixfd1[isx][in]=leftmatrixfd1[isx][in]*weight;
	}
    } 
    qrsolver(leftmatrixfd1,sz,conapfd1,nsx,nfd1);
    *beta1=crealf(conapfd1[0]);      
    *alpha1=crealf(conapfd1[1]);  
    free(sz); free(*leftmatrixfd1); free(leftmatrixfd1); free(conapfd1);
}


void imfd2coesolver(float ep,float dl,int nsxall, float *alpha1,float *alpha2,float *beta1,float *beta2)
{
    float **leftmatrixfd2;
    sf_complex *sz,*conapfd2;
    int nsx,nfd2=4;
    int isx,in;
    float dsx,sx,b2,b4,a2,a4,weight,maxsx;
    maxsx=sqrtf(1.0/(1.0+2.0*ep));
    dsx=maxsx/(nsxall-1);
    sz=sf_complexalloc(nsxall);  conapfd2=sf_complexalloc(nfd2);
    kzvti(sz,dsx,ep,dl,nsxall);
    for(isx=nsxall-1;isx>=0; isx--){
	sx=isx*dsx;
	if (sqrtf(sx*sx/(sx*sx+(sz[isx]+1.0)*(sz[isx]+1.0)))<=sinf(83.5/180.0*SF_PI)) break;
    }
    nsx=isx;
    leftmatrixfd2=sf_floatalloc2(nsx,nfd2);
    implicitleftmatrixisovtifd2(leftmatrixfd2,sz,nsx,dsx);
    for(isx=0;isx<nsx;isx++){
	weight=sqrtf(sx*sx+(sz[isx]+1.0)*(sz[isx]+1.0))/(sz[isx]+1.0);
	sz[isx]=sz[isx]*weight*weight*powf(weight,0.75);
	for(in=0;in<nfd2;in++){
	    leftmatrixfd2[isx][in]=leftmatrixfd2[isx][in]*weight*weight*powf(weight,0.75);
	}
    } 
    qrsolver(leftmatrixfd2,sz,conapfd2,nsx,nfd2);
    b2=crealf(conapfd2[0]); b4=crealf(conapfd2[1]); a2=crealf(conapfd2[2]); a4=crealf(conapfd2[3]);
    *beta1=0.5*(b2+sqrtf(b2*b2-4*b4));     *beta2=0.5*(b2-sqrtf(b2*b2-4*b4));
    *alpha1=(a4-a2*(*beta1))/((*beta2)-(*beta1));  *alpha2=(a4-a2*(*beta2))/((*beta1)-(*beta2));
    free(sz); free(*leftmatrixfd2); free(leftmatrixfd2); free(conapfd2);
}



//-----------------------------------------------------------------------------------------

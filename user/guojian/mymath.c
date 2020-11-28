/* Math functions. */

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
#include  <rsf.h>

#include "mymath.h"
#include "arrayindex.h"

void thryj(sf_complex* a, sf_complex* b, sf_complex* c, sf_complex* v, int h, int yJ)
/*< tridiagonal >*/
{
    int i;
    for(i=0; i<yJ; i++){
	b[i]=b[i]/a[i];
	v[i]=v[i]/a[i];
    }
    for(i=yJ; i<h; i++){
	a[i]=a[i]-b[i-yJ]*c[i-yJ];
	v[i]=(v[i]-c[i-yJ]*v[i-yJ] )/a[i];
	if(i<=h-yJ-1) b[i]=b[i]/a[i];
    }
    for(i=h-yJ-1; i>=0; i--) {
	v[i]=v[i]-b[i]*v[i+yJ];
	/*cout <<"test i="<< i << "value  "<< __real__ v[i]<<"  " << __imag__ v[i]<< "\n";*/
    }
}


void thrnew(sf_complex *a, sf_complex *b,sf_complex *c, 
	    sf_complex *v, int n)
{
    int i;
    float chir,chii;
    sf_complex cdlt;
    b[0]=b[0]/a[0];
    v[0]=v[0]/a[0];
    for (i=1;i<=n-2;i++){
	chir= crealf(a[i])-crealf(b[i-1])*crealf(c[i-1]) 
	    + cimagf(b[i-1])*cimagf(c[i-1]);
	chii= cimagf(a[i])-(crealf(b[i-1])*cimagf(c[i-1]) +
			    cimagf(b[i-1])*crealf(c[i-1]));
	cdlt=sf_cmplx(chir,-chii)/(chir*chir +chii*chii);
	b[i]   = b[i]*cdlt;
	v[i] = (v[i]-c[i-1]*v[i-1])*cdlt;
    }
    v[n-1]=(v[n-1]-c[n-2]*v[n-2])/(a[n-1]-c[n-2]*b[n-2]);
    for (i=n-2;i>=0;i--){
	v[i]=v[i]-b[i]*v[i+1];
    }
}



void p4solver(float *a, sf_complex* z) 
/*< solver >*/
{
    sf_complex y[3];
    float b3,b2,b1,b0,b[4];
    float a4,a3,a2,a1,a0;
    float ctmp1,ctmp2,yy;
    sf_complex  z1,z2,z3,z4,Dtmp,Rtmp,Etmp;
    int i;

    a0=a[0]; a1=a[1]; a2=a[2]; a3=a[3]; a4=a[4];
    //printf("solver %f,%f,%f,%f,%f\n",a0,a1,a2,a3,a4);
    a0=a0/a4; a1=a1/a4; a2=a2/a4; a3=a3/a4; a4=1.0;
    b3=1; b2=-a2; b1=a1*a3-4.0*a0; b0=4.0*a2*a0-pow(a1,2)-pow(a3,2)*a0;
    b[0]=b0; b[1]=b1; b[2]=b2; b[3]=b3;
    p3solver(b,y);
    //printf("solver y(%f,%f)--(%f,%f)--(%f,%f)",__real__ y[0],__imag__ y[0],__real__ y[1],__imag__ y[1],__real__ y[2],__imag__ y[2]);
    for (i=0; i<3; i++){
	if (fabs(cimagf(y[i]))< fabs(crealf(y[i]))*0.00001) {
	    yy=crealf(y[i]);
	    ctmp1=powf(a3,2)-4.0*a2+4.0*yy;
	    Rtmp=0.5*csqrtf( sf_cmplx(ctmp1,0.0));
	    ctmp2=powf(yy,2)-4.0*a0;
	    if (Rtmp==0.0) 
		Dtmp=csqrtf( 3.0/4.0*powf(a3,2) - 2.0*a2+2.0*csqrtf(  sf_cmplx(ctmp2,0.0)  ) );
	    else
		Dtmp=csqrtf( 3.0/4.0*powf(a3,2)-cpowf(Rtmp,2)-2.0*a2+0.25*(4.0*a3*a2-8.0*a1-powf(a3,3))/Rtmp);
	    if (Rtmp==0.0) 
		Etmp=csqrtf( 3.0/4.0*powf(a3,2) - 2.0*a2- 2.0*csqrtf(  sf_cmplx(ctmp2,0.0)  ) );
	    else
		Etmp=csqrtf( 3.0/4.0*powf(a3,2)-cpowf(Rtmp,2)-2.0*a2-0.25*(4.0*a3*a2-8.0*a1-powf(a3,3))/Rtmp);
	    z1=-0.25*a3+0.5*Rtmp+0.5*Dtmp;
	    z2=-0.25*a3+0.5*Rtmp-0.5*Dtmp;
	    z3=-0.25*a3-0.5*Rtmp+0.5*Etmp;
	    z4=-0.25*a3-0.5*Rtmp-0.5*Etmp;
	    z[0]=z1; z[1]=z2; z[2]=z3; z[3]=z4;
	    //printf("solver:z: (%f,%f),(%f,%f),(%f,%f),(%f,%f)\n",__real__ z1,__imag__
//z1,__real__ z2,__imag__ z2,__real__ z3,__imag__ z3,__real__ z4,__imag__ z4);
	    break;
	}
    }
}


void p3solver(float* a, sf_complex* x)
/*< solver >*/
{
    double a0,a1,a2,a3;
    double p,q,aa;
    sf_double_complex  x1,x2,x3,sqrtaa;
    sf_double_complex w3,w,w1,w2;

    a0=a[0]; a1=a[1]; a2=a[2]; a3=a[3];
    //printf("a33333=%f,%f,%f,%f\n",a3,a2,a1,a0);
    a2=a2/a3; a1=a1/a3; a0=a0/a3;
    p=(3.0*a1-a2*a2)/3.0;
    q=(9.0*a1*a2-27*a0-2.0*a2*a2*a2  )/27.0;
    //printf("pppppp=%f,%f\n",p,q);
    aa=0.25*q*q+1.0/27.0*p*p*p;
    //aa=10000.0*aa;
    //printf("aa=%f\n",aa);
    if (aa >=0)
	sqrtaa=sf_dcmplx(sqrt(aa),0.0);
    else
	sqrtaa=sf_dcmplx(0.0,sqrt(-aa));
    //sqrtaa=sqrtaa*0.01;
    w3=0.5*q+sqrtaa;   // csqrt(aa); 
    w=cpowf(w3,1.0/3.0);       //            w=(w3)**(1.0/3.0);  
    //printf("wwwww=(%f,%f)\n",__real__ w,__imag__ w);
    x1=w-p/(3.0*w)-1.0/3.0*a2;
    w1=w*sf_dcmplx( cos(SF_PI*2.0/3.0) , sin(SF_PI*2.0/3.0) );
    x2=w1-p/(3.0*w1)-1.0/3.0*a2;
    w2=w*sf_dcmplx(cos(SF_PI*4.0/3.0),sin(SF_PI*4.0/3.0));
    x3=w2-p/(3.0*w2)-1.0/3.0*a2;
    x[0]=x1;  x[1]=x2; x[2]=x3;
    //printf("solver3p y(%f,%f)--(%f,%f)--(%f,%f),(%f,%f)\n",__real__ x[0],__imag__ x[0],__real__ x[1],__imag__ x[1],__real__ x[2],__imag__ x[2],__real__ sqrtaa,__imag__ sqrtaa);
} 


void qrdcmp(float* a,int m,int n,float* d,int sing)
/*< QR >*/
{
    float scale,sigma,sum,tau;
    float* c;
    int i,j,k;

    c=sf_floatalloc(n); 
    scale=0.0;
    sing=0;

    for (k=0; k<n; k++){ 
	for(i=k; i<m; i++) 
	    if (fabs(a[i2(i,k,n)])>scale) scale=fabs(a[i2(i,k,n)]);
	if (scale==0.0){
	    sing=1;
	    c[k]=0.0; d[k]=0.0;
	}
	else{
	    for(i=k;i<m;i++)  
		a[i2(i,k,n)]=a[i2(i,k,n)]/scale;
	    sum=0.0;
	    for(i=k;i<m;i++)  
		sum=sum+a[i2(i,k,n)]*a[i2(i,k,n)];
	    if (a[i2(k,k,n)]>0) 
		sigma=sqrt(sum);
	    else
		sigma=-sqrt(sum);
	    a[i2(k,k,n)]=a[i2(k,k,n)]+sigma;
	    c[k]=sigma*a[i2(k,k,n)];
	    d[k]=-scale*sigma;
	    for(j=k+1;j<n;j++){  
		sum=0.0;
		for(i=k;i<m;i++) 
		    sum=sum+a[i2(i,k,n)]*a[i2(i,j,n)];
		tau=sum/c[k];
		for(i=k;i<m;i++) 
		    a[i2(i,j,n)]=a[i2(i,j,n)]-tau*a[i2(i,k,n)];
	    }
	}
    }
    free(c); 

}


void qrdpost(float* a,int m,int n,float* d,float* q,float* r)
/*< QR >*/
{
    float* qt; 
    float  *qtr,*atr;
    int k,l,i,j; 
    float con,akicon;   


    qt=sf_floatalloc(m*m);  
    qtr=sf_floatalloc(m*m); atr=sf_floatalloc(m*n);
    for(i=0;i<n;i++)
	for(j=0;j<m;j++)
	    atr[i2(i,j,m)]=a[i2(j,i,n)];
    for(k=0;k<n;k++){   
	for(l=0;l<n;l++){  
	    if (l>k){ 
		r[i2(k,l,n)]=a[i2(k,l,n)];
	    }
	    else{
		if (l<k) 
		    r[i2(k,l,n)]=0.0;
		else
		    r[i2(k,l,n)]=d[k];
	    }
	}
    }
//printf("post111111\n");
    for(k=0;k<m;k++){  
	for(l=0;l<m;l++){   
	    if (l==k) 
		q[i2(k,l,m)]=1.0;
	    else
		q[i2(k,l,m)]=0.0;
	}
    }
    for (i=0;i<m;i++)
	for(j=0;j<m;j++)
	    qtr[i2(i,j,m)]=q[i2(j,i,m)];
//printf("post222222\n");

    for(i=n-1; i>=0; i--){ 
	//printf("i=%d\n",i);
	con=0.0;
	for(k=i;k<m;k++)   
	    con=con+atr[i2(i,k,m)]*atr[i2(i,k,m)];
	con=con/2.0;

  
	for(k=i;k<m;k++){ 
	    akicon=a[i2(k,i,n)]/con; 
	    for(l=i;l<n;l++){ 
		qt[i2(k,l,m)]=0.0;
		for(j=i;j<m;j++)   
		    qt[i2(k,l,m)]+=qtr[i2(l,j,m)]*atr[i2(i,j,m)];
		qt[i2(k,l,m)]=qt[i2(k,l,m)]*akicon;
	    }
	}

	for(k=i;k<m;k++){ 
	    for(l=i;l<m;l++)  
		qtr[i2(l,k,m)]=qtr[i2(l,k,m)]-qt[i2(k,l,m)];
	}
    }
    for (i=0;i<m;i++)
	for(j=0;j<m;j++)
	    q[i2(i,j,m)]=qtr[i2(j,i,m)];
    free(qt); free(qtr); free(atr); 

}


void qrcon(float* q,float* r,
	   sf_complex* dpp,sf_complex* conap,int m,int n)
/*< QR >*/
{
    int i_n,i_m; 
    for(i_n=0; i_n<n; i_n++){ 
	conap[i_n]=sf_cmplx(0.0,0.0);
	for(i_m=0;i_m<m;i_m++) 
	    conap[i_n]=conap[i_n] +q[i2(i_m,i_n,n)]*dpp[i_m];
    }
    uptrisolver(r,conap,n);

}

void uptrisolver(float* r,sf_complex* conap,int n)
/*< up tri solver >*/
{
    int i_n,i_m;

    for(i_n=n-1;i_n>=0;i_n--){  
	for(i_m=i_n+1; i_m<n;i_m++)   
	    conap[i_n]=conap[i_n]-conap[i_m]*r[i2(i_n,i_m,n)];
	conap[i_n]=conap[i_n]/r[i2(i_n,i_n,n)];
    }  

}


sf_complex chooseproot(sf_complex *allroot)
/*< choose positive root >*/
{
    sf_complex  z1,z2,z3,z4,x1;
    sf_complex  z[4],xx1,xx2,tmpxx,zz[2];
    int i,j,k;
    float  absi,absx;


    z1=allroot[0]; z2=allroot[1]; z3=allroot[2]; z4=allroot[3];
    z[0]=z1; z[1]=z2; z[2]=z3; z[3]=z4;


    for(i=0,j=0;i<4;i++)  //do i=1,4
	if (fabsf(cimagf(z[i])) < 0.00001*fabsf(crealf(z[i])))
	    j++;  //     j=j+1;
   

    if (j==0 || j==2){ 
	// no real root
	absi=fabsf(cimagf(z[1]));
	for(i=0;i<4;i++)  //do i=1,4
	    if (absi <fabsf(cimagf(z[i])) ) 
		absi=fabsf(cimagf(z[i]));
	k=0;
	for(i=0;i<4;i++){   //do i=1,4
	    if (absi==fabsf(cimagf(z[i])) ){
		k++; //k=k+1
		zz[k-1]=z[i];
	    }
	}
    }

    if (j==4){ 
	// 4 real root
	absx=fabsf(crealf(z[0]));
	xx1=z[0];
	for(i=0;i<4;i++){   //do i=1,4
	    if (absx >fabsf(crealf(z[i])) ) {
		absx=fabsf(crealf(z[i]));
		xx1=z[i];
	    }
	}



	if ( xx1==z[0] ){
	    absx=fabsf(crealf(z[1]));
	    xx2=z[1];
	}
	else{
	    absx=fabsf(crealf(z[0]));
	    xx2=z[0];
	}
  
	for(i=0;i<4;i++){ //do i=1,4
	    if (absx > fabsf(crealf(z[i])) && xx1!=z[i] ){ 
		absx=fabsf(crealf(z[i]));
		xx2=z[i];
	    }
	}

	if (crealf(xx1) < 0 ){ 
	    tmpxx=xx1;
	    xx1=xx2;
	    xx2=tmpxx;
	}
	zz[0]=xx1;
	zz[1]=xx2;
	k=0; 
	for (i=0;i<4;i++){
	    if (crealf(z[i])*10000.0<0) k++;
	}
	if (k>2){ 
	    zz[0]=sf_cmplx(0.0,1.0);
	    zz[1]=sf_cmplx(0.0,1.0);
	}
	k=0;
	for (i=0;i<4;i++){
	    if (crealf(z[i])*10000.0>0) k++;
	}
	if (k>2){ 
	    zz[0]=sf_cmplx(0.0,1.0);
	    zz[1]=sf_cmplx(0.0,1.0);
	}

    }
  
    x1=zz[0];
    return x1;
}

void stoepd (int n, double r[], double g[], double f[], double a[])
{
    int i,j;
    double v,e,c,w,bot;
    if (r[0] == 0.0) return;
    a[0] = 1.0;
    v = r[0];
    f[0] = g[0]/r[0];
    for (j=1; j<n; j++) {
	/* solve Ra=v as in Claerbout, FGDP, p. 57 */
	a[j] = 0.0;
	f[j] = 0.0;
	for (i=0,e=0.0; i<j; i++)
	    e += a[i]*r[j-i];
	c = e/v;
	v -= c*e;
	for (i=0; i<=j/2; i++) {
	    bot = a[j-i]-c*a[i];
	    a[i] -= c*a[j-i];
	    a[j-i] = bot;
	}
	/* use a and v above to get f[i], i = 0,1,2,...,j */
	for (i=0,w=0.0; i<j; i++)
	    w += f[i]*r[j-i];
	c = (w-g[j])/v;
	for (i=0; i<=j; i++)
	    f[i] -= c*a[j-i];
    }
}



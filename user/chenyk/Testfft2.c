/* Test fft along second axis */
#include<rsf.h>
#include"fft.h"
int main(void)
{
	
    int i1,i2,nt,nx,nxfft;

 	kiss_fft_cpx **a,**b,**c;
    a = (kiss_fft_cpx**) sf_complexalloc2(4,4);
    b = (kiss_fft_cpx**) sf_complexalloc2(4,4);
    c = (kiss_fft_cpx**) sf_complexalloc2(4,4);

	a[0][0]=cmplx(1.0,0.0);
	a[0][1]=cmplx(1.5,0.0);
	a[0][2]=cmplx(2.0,0.0);
	a[0][3]=cmplx(3.0,0.0);
	a[1][0]=cmplx(1.0,0.0);
	a[1][1]=cmplx(1.5,0.0);
	a[1][2]=cmplx(2.0,0.0);
	a[1][3]=cmplx(3.0,0.0);
	a[2][0]=cmplx(1.0,0.0);
	a[2][1]=cmplx(1.5,0.0);
	a[2][2]=cmplx(2.0,0.0);
	a[2][3]=cmplx(3.0,0.0);
	a[3][0]=cmplx(1.0,0.0);
	a[3][1]=cmplx(1.5,0.0);
	a[3][2]=cmplx(2.0,0.0);
	a[3][3]=cmplx(3.0,0.0);

	/*for(i1=0;i1<4;i1++)
			kiss_fft_stride(forw,a[0]+i1,b[0]+i1,4); */

	/*for(i1=0;i1<4;i1++)
		{kiss_fft_stride(invs,b[0]+i1,temp,4);
		for(i2=0;i2<4;i2++)
				c[i2][i1]=temp[i2];} */

	nt=4; nx=4;
   	nxfft = 2*kiss_fft_next_fast_size((nx+1)/2);
	xfft(a,b,nt,nxfft);
	ixfft(b,c,nt,nxfft);

	
	for(i2=0;i2<4;i2++) 
		for(i1=0;i1<4;i1++)
			sf_warning("a[%d][%d]=%g+i %g",i2,i1,sf_crealf(a[i2][i1]),sf_cimagf(a[i2][i1]));
	for(i2=0;i2<4;i2++) 
		for(i1=0;i1<4;i1++)
			sf_warning("b[%d][%d]=%g+i %g",i2,i1,sf_crealf(b[i2][i1]),sf_cimagf(b[i2][i1]));
	for(i2=0;i2<4;i2++) 
		for(i1=0;i1<4;i1++) 
			sf_warning("c[%d][%d]=%g+i %g",i2,i1,sf_crealf(c[i2][i1]),sf_cimagf(c[i2][i1]));

    exit(0);
}

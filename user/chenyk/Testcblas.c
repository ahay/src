/* Test cblas */
#include<rsf.h>

int main(void)
{
	
     int i, N=4, inca=1, incb=1;
     float *a, *b, *bb, alpha=0.5;
     a=sf_floatalloc(N);
     b=sf_floatalloc(N);
     bb=sf_floatalloc(N);

     a[0]=1; a[1]=2; a[2]=3; a[3]=4; 
     b[0]=4; b[1]=3; b[2]=2; b[3]=1; 
     bb[0]=4; bb[1]=3; bb[2]=2; bb[3]=1; 
  
     for(i=0;i<N;i++)
     	sf_warning("a[%d]=%f",i,a[i]);
     for(i=0;i<N;i++)
     	sf_warning("b[%d]=%f",i,b[i]);

     cblas_saxpy(N, 0.5 , a, 1, b, 1);   

     for(i=0;i<N;i++)
     	sf_warning("(0.5*a+b)[%d]=%f",i,b[i]);

     saxpy_(&N, &alpha, a, &inca, bb, &incb);
 
     for(i=0;i<N;i++)
     	sf_warning("(0.5*a+bb)[%d]=%f",i,bb[i]);

    exit(0);
}

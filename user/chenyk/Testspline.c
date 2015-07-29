/* Test for 1-D B-spline interpolation. */
#include<rsf.h>

int main(void)
{
	
    int i, j, i1, nf=8, n1=101, n2=301, im;
    float *a, *a1, *a2, *w, f;   /* a is reference data, a1 is true data, a2 is interpolated data. */

    a=sf_floatalloc(n1);
    a1=sf_floatalloc(n2);
    a2=sf_floatalloc(n2);
    w=sf_floatalloc(nf);  
   
    for(i=0;i<n1;i++)
	a[i]=sinf(SF_PI/100*i)+cosf(SF_PI/100*i); /* Calculate original data point value. */
    
    for(i=0;i<n2;i++)
        a1[i]=sinf(SF_PI/100*(0+1.0/3*i))+cosf(SF_PI/100*(0+1.0/3*i)); /* Calculate true fine-grid data point value. */

    for(i1=0;i1<n2;i1++)
    {
	f=(1.0-0.5*nf)+(0+1.0/3*i1-0)/1;
	j=f;  /* j is grid location */
	f-=j; /* f is offset from grid location */
        sf_spline_int(f,nf,w);
	a2[i1]=0;
	for (i = SF_MAX(0,-j); i < SF_MIN(nf,n1-j); i++) { 
	    im = i+j;
	    a2[i1]+=a[im]*w[i];
	}
    }

    for(i=0;i<n2;i++)
	sf_warning("Interpolation result comparison: True[%d]=%f, Inter[%d]=%f",i,a1[i],i,a2[i]);

    exit(0);
}

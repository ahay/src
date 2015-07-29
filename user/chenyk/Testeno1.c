/* Test 1D eno interpolation, sf_eno_init, sf_eno_set, sf_eno_apply. */
#include<rsf.h>

int main(void)
{
	
    int i, j, i1, order=8, n1=11, n2=31;
    float *a, *a1, *a2, f, f1;   /* a is reference data, a1 is true data, a2 is interpolated data. */
    sf_eno map;

    a=sf_floatalloc(n1);
    a1=sf_floatalloc(n2);
    a2=sf_floatalloc(n2);

    for(i=0;i<n1;i++)
	a[i]=sinf(SF_PI/10*i)+cosf(SF_PI/10*i); /* Calculate original data point value. */
    
    for(i=0;i<n2;i++)
        a1[i]=sinf(SF_PI/10*(0+1.0/3*i))+cosf(SF_PI/10*(0+1.0/3*i)); /* Calculate true fine-grid data point value. */
	
    map = sf_eno_init(order,n1);
    sf_eno_set (map,a);

    for(i1=0;i1<n2;i1++)
    {
	f=(0+1.0/3*i1-0)/1;
	j=f;  /* j is grid location */
	f-=j; /* f is offset from grid location */
        sf_eno_apply(map, j, f, a2+i1, &f1, FUNC); /* Interpolate fine-grid data point. f1 is the derivative sequence. FUNC is a flag value for eno interpolation */
    }

    for(i=0;i<n2;i++)
	sf_warning("Interpolation result comparison: True[%d]=%f, Inter[%d]=%f",i,a1[i],i,a2[i]);

    exit(0);
}

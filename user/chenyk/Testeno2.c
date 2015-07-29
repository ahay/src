/* Test 2D eno interpolation, sf_eno2_init, sf_eno2_set, sf_eno2_apply. */
#include<rsf.h>

int main(void)
{
	
    int i, j, m, n, order=6, n1x=11, n1y=11, n2x=31, n2y=31; /*for 2D eno, order must be less than 7, or ng=2*(order-1) is too big */
    float **a, **a1, **a2, x, y, f1;   /* a is reference data, a1 is true data, a2 is interpolated data. */
    sf_eno2 map;

    a=sf_floatalloc2(n1x,n1y); /*n1x and n1y are original total grid points*/
    a1=sf_floatalloc2(n2x,n2y);/*n2x and n2y are interpolated total grid points*/
    a2=sf_floatalloc2(n2x,n2y);

    for(i=0;i<n1x;i++)
	for(j=0;j<n1y;j++)
		a[j][i]=sinf(SF_PI/100*i*j)+cosf(SF_PI/100*i*j); /* Calculate original data point value. */
    
    for(i=0;i<n2x;i++)
	for(j=0;j<n2y;j++)
		{a1[j][i]=sinf(SF_PI/100*(0+1.0/3*i)*(0+1.0/3*j))+cosf(SF_PI/100*(0+1.0/3*i)*(0+1.0/3*j));} /* Calculate true fine-grid data point value. */
	
    map = sf_eno2_init(order,n1x,n1y);
    sf_eno2_set (map,a);
    sf_warning("hello1");
    for(i=0;i<n2x;i++)
	for(j=0;j<n2y;j++)
    {
	x=(0+1.0/3*i-0)/1;
	y=(0+1.0/3*j-0)/1;
	m=x;  /* m is grid location in x direction */
	n=y;    /* n is grid location in y direction */
	x-=m; /* x is offset from grid location m */
	y-=n; /* y is offset from grid location n */
        sf_eno2_apply (map, m, n, x, y, &a2[j][i], &f1, FUNC); /* Interpolate fine-grid data point. f1 is the derivative sequence. FUNC is a flag value for eno interpolation */
    }
    sf_warning("hello2");
    for(i=0;i<n2x;i++)
	for(j=0;j<n2y;j++)
		{sf_warning("Interpolation comparison: True[%d][%d]=%f, Inter[%d][%d]=%f",j,i,a1[j][i],j,i,a2[j][i]);}

    exit(0);
}

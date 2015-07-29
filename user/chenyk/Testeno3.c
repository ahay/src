/* Test 3D eno interpolation, sf_eno3_init, sf_eno3_set, sf_eno3_apply. */
#include<rsf.h>

int main(void)
{
	
    int i, j, k, m, n, p, order=6, n1x=11, n1y=11, n1z=11, n2x=31, n2y=31, n2z=31; /*for 3D eno, order must be less than 7, or ng=2*(order-1) is too big */
    float ***a, ***a1, ***a2, x, y, z, f1;   /* a is reference data, a1 is true data, a2 is interpolated data. */
    sf_eno3 map;

    a=sf_floatalloc3(n1x,n1y,n1z); /*n1x, n1y and n1z are original total grid points*/
    a1=sf_floatalloc3(n2x,n2y,n2z);/*n2x, n2y and n2z are interpolated total grid points*/
    a2=sf_floatalloc3(n2x,n2y,n2z);

    for(i=0;i<n1x;i++)
	for(j=0;j<n1y;j++)
	     for(k=0;k<n1z;k++)
		  a[k][j][i]=sinf(SF_PI/1000*i*j*k)+cosf(SF_PI/1000*i*j*k); /* Calculate original data point value. */
    
    for(i=0;i<n2x;i++)
	for(j=0;j<n2y;j++)
	     for(k=0;k<n2z;k++)
		  a1[k][j][i]=sinf(SF_PI/1000*(0+1.0/3*i)*(0+1.0/3*j)*(0+1.0/3*k))+cosf(SF_PI/1000*(0+1.0/3*i)*(0+1.0/3*j)*(0+1.0/3*k)); /* Calculate true fine-grid data point value. */
	
    map = sf_eno3_init(order,n1x,n1y,n1z);
    sf_eno3_set (map,a);

    for(i=0;i<n2x;i++)
	for(j=0;j<n2y;j++)
	    for(k=0;k<n2z;k++)
    	    {
		x=(0+1.0/3*i-0)/1;
		y=(0+1.0/3*j-0)/1;
		z=(0+1.0/3*k-0)/1;
		m=x;  /* m is grid location in x direction */
		n=y;  /* n is grid location in y direction */
		p=z;  /* p is grid location in z direction */
		x-=m; /* x is offset from grid location m */
		y-=n; /* y is offset from grid location n */
		z-=p; /* z is offset from grid location p */
        	sf_eno3_apply (map, m, n, p, x, y, z, &a2[k][j][i], &f1, FUNC); /* Interpolate fine-grid data point. f1 is the derivative sequence. FUNC is a flag value for eno interpolation */ 
    	    }
    for(i=0;i<n2x;i++)
	for(j=0;j<n2y;j++)
	    for(k=0;k<n2z;k++)
		sf_warning("Comparison: True[%d][%d][%d]=%f, Inter[%d][%d][%d]=%f",k,j,i,a1[k][j][i],k,j,i,a2[k][j][i]);

    exit(0);
}


#include <rsf.h>
#include <math.h>

float fd1_1(float ***p, int i1, int i2, int i3, char id)
/*< first order difference >*/
{
	int j1, j2, j3;
	double a=0;
	

	switch(id)
	{
	case 'z':
		for (j3=-1; j3<=1; j3++)
		for (j2=-1; j2<=1; j2++)
			a += p[i3+j3][i2+j2][i1+1] - p[i3+j3][i2+j2][i1-1];
		break;
	case 'y':
		for (j2=-1; j2<=1; j2++)
		for (j1=-1; j1<=1; j1++)
			a += p[i3+1][i2+j2][i1+j1] - p[i3-1][i2+j2][i1+j1];
		break;
	case 'x':
		for (j3=-1; j3<=1; j3++)
		for (j1=-1; j1<=1; j1++)
			a += p[i3+j3][i2+1][i1+j1] - p[i3+j3][i2-1][i1+j1];
	}
	return ((float)(a/18.0));
}

float fd1_2(float **p, int i2, int i3, char id)
/*< first and second order FD >*/
{
        int  j2, j3, j;
	double a=0;
	float c0=0.0;

	switch(id)
	{
	case 'u':
		for (j=-1; j<=1; j++)	
			a += p[i3-j][i2+1] - p[i3-j][i2-1];
		        c0=6.0;
		break;
        case 'v':
		for (j=-1; j<=1; j++)	
			a += p[i3+1][i2-j] - p[i3-1][i2-j];
		        c0=6.0;
		break; 
	case 'U':
		for (j3=-1; j3<=1; j3++)
		for (j2=-1; j2<=1; j2++)
		         a += (-2*pow(-0.5,j2*j2))*(p[i3-j3][i2+j2]);
		         c0=6.0;
		break;
	case 'V':
                for (j3=-1; j3<=1; j3++)
		for (j2=-1; j2<=1; j2++)
		         a += (-2*pow(-0.5,j3*j3))*(p[i3+j3][i2-j2]);
		         c0=6.0;
		break;
        case 'W':
	    a = p[i3+1][i2+1]+p[i3-1][i2-1]-p[i3+1][i2-1]-p[i3-1][i2+1];
	    c0=4.0;

	}
	return ((float)(a/c0));
}    

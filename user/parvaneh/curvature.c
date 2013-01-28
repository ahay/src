
#include <rsf.h>

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



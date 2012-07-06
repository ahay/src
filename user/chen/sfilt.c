/* simple filters  */
#include <rsf.h>

typedef float (*sfilt)(int,float*);
/* generic filter interface */
/*^*/


float sfilt_mean(int n, float *x)
/*< mean filter >*/
{
	float d=0.0;
	do{
		n--;
		d += x[n-1];
	}while(n>0);
	return (d/n);
}

float sfilt_median(int n, float *p)
/*< median filter >*/
{
	int i1, j1, chg;
	float temp;

	for(j1=n-1; j1>=n/2; j1--)
	{
		chg=0;
		for(i1=0; i1<j1; i1++)
		if(p[i1] > p[i1+1])
		{
			temp = p[i1];
			p[i1] = p[i1+1];
			p[i1+1] = temp;
			chg=1;
		}
		if(chg==0) break;
	}
	return (p[n/2]);
}




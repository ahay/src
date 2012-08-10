/* singular value decomposition */
#include <rsf.h>
#include "_lapack.h"


void svd( int n1, int n2, float **a, bool vt,
	float **u, float *s, float **v)
/*< singular value decomposition >*/
{
	int i1, i2;
	int lwork, info;
	float e, *work;

	lwork = n1*n2;
	work = sf_floatalloc(lwork);

	sgesvd_("A", "A", &n1, &n2, a[0], &n1, 
		s, u[0], &n1, v[0], &n2,
		work, &lwork, &info );
	free(work);

	if(vt == false)
	{
		for(i1=0; i1<n2; i1++)
		for(i2=0; i2<i1; i2++)
		{
			e = v[i2][i1];
			v[i2][i1] = v[i1][i2]; 
			v[i1][i2] = e;
		}
	}
}




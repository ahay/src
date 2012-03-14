
#include <rsf.h>

static sf_complex ***mrt;
static int n1, n2, n3; // offset stepout freq

void sf_fart_init(int curv, bool inv,
	float op, float dp, int np,
	float ox, float dx, int nx, 
	float of, float df, int nf)
/*< initialize >*/
{
	int i1,i2,i3;
	double x1, x3, *px, *pp;

	n1 = nx;
	n2 = np;
	n3 = nf;

	px = (double*)sf_alloc(nx, sizeof(double));
	for(i1=0;i1<nx;i1++)
	{
		x1=ox+i1*dx;
		switch(curv)
		{
		case 1:
			px[i1] = x1*x1;
			break;
		default:
			px[i1] = x1;
		}
	}
	pp = (double*)sf_alloc(np, sizeof(double));
	for(i1=0;i1<np;i1++) pp[i1] = op+i1*dp;

	mrt = (sf_complex***) sf_alloc(nf, sizeof(sf_complex**));

	for(i3=0;i3<nf;i3++)
	{
		x3 = 2*M_PI*(of+i3*df);
		mrt[i1] = sf_complexalloc2(nx, np);
		for(i2=0;i2<np;i2++)
		{
			for(i1=0;i1<nx;i1++)
			{
				x1 = x3*px[i1]*pp[i2];
				mrt[i3][i2][i1] = cos(x1) + I*sin(x1);
			}
		}
	}
}

void sf_fart(sf_complex **in, sf_complex **out)
/*< forward Radon transform>*/
{
	int i1, i2, i3;
	for(i3=0;i3<n3;i3++)
	{
		for(i2=0;i2<n2;i2++)
		{
			out[i2][i3] = 0.0;
			for(i1=0;i1<n1;i1++)
				out[i2][i3] += mrt[i3][i2][i1]*in[i1][i3];
		}
	}
}

void sf_ifart(sf_complex **in, sf_complex **out)
/*< inverse Radon transform>*/
{
	int i1, i2, i3;
	for(i3=0;i3<n3;i3++)
	{
		for(i1=0;i1<n1;i1++)
		{
			out[i1][i3] = 0.0;
			for(i2=0;i2<n2;i2++)
				out[i1][i3] += conj(mrt[i3][i2][i1])*in[i2][i3];
		}
	}
}

void sf_fart_release()
/*< release memory >*/
{
	int i3;

	for(i3=0;i3<n3;i3++)
	{
		free(mrt[i3][0]);
		free(mrt[i3]);
	}
	free(mrt);
}


/* step split phase shift */

#include <rsf.h>

static int nw,nx;
static float dz;
static float *vel;
static int *k2;

void sf_psss_init(int nw0,int nx0,int nz,
	float dw,float dx, float dz0, float *v0)
/*< initialize >*/
{
	int i;
	nw=nw0;
	nx=nx0;

	dz=dz0*2.0*M_PI/(dx*nx);

	vel=(float*)sf_floatalloc(nz);
	k2=(int*)sf_intalloc(nx);

	for(i=0; i<nx/2+1; i++)		k2[i] = i*i;
	for(i=nx/2+1; i<nx; i++) 	k2[i] = (i-nx)*(i-nx);

	for(i=0; i<nz; i++)
	{
		vel[i] = dw*dx*nx/v0[i];
		vel[i] *= vel[i];
	}
}

void sf_psss_step(int iz, sf_complex **io) 
/*< step in depth >*/
{
	int iw,ix;
	float t1;
	sf_complex a;

	for(ix=0;ix<nx;ix++)
	{
		for(iw=0;iw<nw;iw++)	
		{
			t1  = (vel[iz]*iw*iw-k2[ix]);
			if(t1>0) 
			{
				t1 = sqrt(t1)*dz;
				a = sf_cmplx(cos(t1),sin(t1));
			}else	a = 0.0;
			io[ix][iw] *= a;
		}
	}
}

void sf_psss_close()
/*< free allocated storage >*/
{
	free(vel);
}


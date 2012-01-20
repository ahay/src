/* step split phase shift */

#include <rsf.h>

static int nw,nx;
static float dz,dx;
static float *vel;

void sf_psss_init(int nw0,int nx0,int nz,
	float dw,float dx, float dz0, float **v0)
/*< initialize >*/
{
	int iz,ix;
	nw=nw0;
	nx=nx0;

	dz=dz0;

	vel=(float*)sf_floatalloc(nz);

	for(iz=0;iz<nz;iz++)
	{
		vel[iz] = 0.0;
		for(ix=0;ix<nx;ix++)	vel[iz] += v0[ix][iz];
		vel[iz] = dw/(vel[iz]*dx*nx);
		vel[iz] *= vel[iz];
	}
}

void sf_psss_step(int iz,kiss_fft_cpx **io) 
/*< step in depth >*/
{
	int iw,ix;
	kiss_fft_cpx a,b;

	for(ix=0;ix<nx;ix++)
	{
		for(iw=0;iw<nw;iw++)	
		{
			a.r = io[ix][iw].r;
			a.i = io[ix][iw].i;
			b.r = sqrt(vel[iz]*iw*iw-ix*ix)*dz/(dx*nx);
			b.i = sin(b.r);
			b.r = cos(b.r);
			
			io[ix][iw].r = a.r*b.r-a.i*b.i;
			io[ix][iw].i = a.r*b.i+a.i*b.r;
		}
	}
}

void sf_psss_exit()
/*< free allocated storage >*/
{
	free(vel);
}


/* step split phase shift */

#include <rsf.h>
#include <malloc.h>

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

	vel=(float*)malloc(nz*sizeof(float));

	for(iz=0;iz<nz;iz++)
	{
		vel[iz] = 0.0;
		for(ix=0;ix<nx;ix++)	vel[iz] += v0[ix][iz];
		vel[iz] = dw/(vel[iz]*dx*nx);
		vel[iz] *= vel[iz];
	}
}

void sf_psss_step(int iz,float **io) 
/*< step in depth >*/
{
	int iw,ix;
	float ar,ai,br,bi;

	for(ix=0;ix<nx;ix++)
	{
		for(iw=0;iw<nw;iw++)	
		{
			ar = io[ix][iw*2];
			ai = io[ix][iw*2+1];
			br = sqrt(vel[iz]*iw*iw-ix*ix)*dz/(dx*nx);
			bi = sin(br);
			br = cos(br);
			
			io[ix][iw*2]   = ar*br-ai*bi;
			io[ix][iw*2+1] = ar*bi+ai*br;
		}
	}
}

void sf_psss_exit()
/*< free allocated storage >*/
{
	free(vel);
}


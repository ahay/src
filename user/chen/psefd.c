/* explicit finite difference phase shift */

#include <rsf.h>

static int nw,nx;
static float f0,dz,dx2;
static float **vel;
static kiss_fft_cpx *buf;

void sf_psefd_init(int nx0,int nw0,
	float dx,float dz0,
	float f00,float **v0)
/*< initialize >*/
{
	nx=nx0;
	nw=nw0;

	dz=dz0;
	f0=f00;
	vel=v0;

	dx2=(dx*dx);
	buf=(kiss_fft_cpx*)sf_complexalloc(nx);
}

void sf_psefd_step3(int iz,kiss_fft_cpx **io) 
/*< step in depth >*/
{
	int iw,ix;
	float t1,a0i,a1i,a2i, br,bi,cr,ci;
	//	0	a0i		1	a1i		0	a2i
//	a0r=0.0;
//	a1r=1.0;
//	a2r=0.0;

	for(iw=1;iw<nw;iw++)  // frequency slice
	{
		t1=dz/(iw*dx2*f0);
		// ix=0 boundary
		ix=0;
		a1i=-t1*vel[ix][iz];
		a2i=0.5*t1*vel[ix+1][iz];
		br=f0*dz*iw/vel[ix][iz];
		bi=sin(br);
		br=cos(br);
		cr=	+ io[iw][ix].r
			- io[iw][ix].i   * a1i
			- io[iw][ix+1].i * a2i;
		ci=	+ io[iw][ix].i
			+ io[iw][ix].r   * a1i
			+ io[iw][ix+1].r * a2i;
		buf[ix].r = cr*br - ci*bi;	
		buf[ix].i = cr*bi + ci*br;		

		for(ix=1;ix<nx-1;ix++)
		{        
			a0i=0.5*t1*vel[ix-1][iz];
			a1i=-t1*vel[ix][iz];
			a2i=0.5*t1*vel[ix+1][iz];
			br=f0*dz*iw/vel[ix][iz];
			bi=sin(br);
			br=cos(br);
			cr=	  io[iw][ix].r 
				- io[iw][ix-1].i * a0i
				- io[iw][ix].i   * a1i
				- io[iw][ix+1].i * a2i;
			ci=	+ io[iw][ix].i
				+ io[iw][ix-1].r * a0i
				+ io[iw][ix].r   * a1i
				+ io[iw][ix+1].r * a2i;
			buf[ix].r = cr*br - ci*bi;	
			buf[ix].i = cr*bi + ci*br;		
		}

		//ix=nx-1 boundary
		a0i=0.5*t1*vel[ix-1][iz];
		a1i=-t1*vel[ix][iz];
		br=f0*dz*iw/vel[ix][iz];
		bi=sin(br);
		br=cos(br);
		cr=	+ io[iw][ix].r
			- io[iw][ix-1].i * a0i
			- io[iw][ix].i   * a1i;
		ci=	+ io[iw][ix].i 
			+ io[iw][ix-1].r * a0i
			+ io[iw][ix].r   * a1i;
		buf[ix].r = cr*br - ci*bi;	
		buf[ix].i = cr*bi + ci*br;	
		
		for(ix=0;ix<nx;ix++)
		{
			io[iw][ix].r = buf[ix].r;
			io[iw][ix].i = buf[ix].i;
		}
	}
	for(ix=0;ix<nx;ix++)
	{
		io[0][ix].r = 0.0;
		io[0][ix].i = 0.0;
	}
}

void sf_psefd_exit()
/*< free allocated storage >*/
{
	free(buf);
}


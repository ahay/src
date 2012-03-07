/* explicit finite difference phase shift */

#include <rsf.h>

static int nw,nx;
static float df,dz,dx2;
static float **vel;
static sf_complex *buf;

void sf_psefd_init(int nx0,int nw0,
	float dx,float dz0,
	float f00,float **v0)
/*< initialize >*/
{
	nx=nx0;
	nw=nw0;

	dz=dz0;
	df=2.0*M_PI*f00;
	vel=v0;

	dx2=(dx*dx);
	buf=(sf_complex*)sf_complexalloc(nx);
}

void sf_psefd_step3(int iz, sf_complex **io) 
/*< step in depth >*/
{
#ifdef SF_HAS_COMPLEX_H			// sf_complex = float complex
	int iw,ix;
	float t1, t2;
	sf_complex a0,a1,a2,b;
	//	0	a0i		1	a1i		0	a2i
//	a0r=0.0;
//	a1r=1.0;
//	a2r=0.0;

	for(iw=1; iw<nw; iw++)  // frequency slice
	{
		t1=dz/(iw*dx2*df);
		// ix=0 boundary
		ix=0;
		a1 = 1.0 - I*t1*vel[ix][iz];
		a2 = 0 + I*0.5*t1*vel[ix+1][iz];
		t2 = df*dz*iw/vel[ix][iz];
		b = cos(t2) + I*sin(t2);
		buf[ix] = b * (io[ix][iw]*a1 + io[ix+1][iw]*a2);

		for(ix=1;ix<nx-1;ix++)
		{
			a0 = 0 + I*0.5*t1*vel[ix-1][iz];
			a1 = 1.0 - I*t1*vel[ix][iz];
			a2 = 0 + I*0.5*t1*vel[ix+1][iz];
			t2 = df*dz*iw/vel[ix][iz];
			b = cos(t2) + I*sin(t2);
			buf[ix] = b * (io[ix-1][iw]*a0 + io[ix][iw]*a1 + io[ix+1][iw]*a2);
		}

		//ix=nx-1 boundary
		a0 = 0 + I*0.5*t1*vel[ix-1][iz];
		a1 = 1.0 - I*t1*vel[ix][iz];
		t2 = df*dz*iw/vel[ix][iz];
		b = cos(t2) + I*sin(t2);
		buf[ix] = b * (io[ix-1][iw]*a0 + io[ix][iw]*a1);
		
		for(ix=0;ix<nx;ix++)
			io[ix][iw] = buf[ix];
	}
//	for(ix=0;ix<nx;ix++)
//		io[ix][0] = 0.0;
#else
	int iw,ix;
	float t1,a0i,a1i,a2i, br,bi,cr,ci;
	//	0	a0i		1	a1i		0	a2i
//	a0r=0.0;
//	a1r=1.0;
//	a2r=0.0;

	for(iw=1;iw<nw;iw++)  // frequency slice
	{
		t1=dz/(iw*dx2*df);
		// ix=0 boundary
		ix=0;
		a1i=-t1*vel[ix][iz];
		a2i=0.5*t1*vel[ix+1][iz];
		br=df*dz*iw/vel[ix][iz];
		bi=sin(br);
		br=cos(br);
		cr=	+ io[ix][iw].r
			- io[ix][iw].i   * a1i
			- io[ix+1][iw].i * a2i;
		ci=	+ io[ix][iw].i
			+ io[ix][iw].r   * a1i
			+ io[ix+1][iw].r * a2i;
		buf[ix].r = cr*br - ci*bi;	
		buf[ix].i = cr*bi + ci*br;		

		for(ix=1;ix<nx-1;ix++)
		{
			a0i=0.5*t1*vel[ix-1][iz];
			a1i=-t1*vel[ix][iz];
			a2i=0.5*t1*vel[ix+1][iz];
			br=df*dz*iw/vel[ix][iz];
			bi=sin(br);
			br=cos(br);
			cr=	  io[ix][iw].r 
				- io[ix-1][iw].i * a0i
				- io[ix][iw].i   * a1i
				- io[ix+1][iw].i * a2i;
			ci=	+ io[ix][iw].i
				+ io[ix-1][iw].r * a0i
				+ io[ix][iw].r   * a1i
				+ io[ix+1][iw].r * a2i;
			buf[ix].r = cr*br - ci*bi;	
			buf[ix].i = cr*bi + ci*br;		
		}

		//ix=nx-1 boundary
		a0i=0.5*t1*vel[ix-1][iz];
		a1i=-t1*vel[ix][iz];
		br=df*dz*iw/vel[ix][iz];
		bi=sin(br);
		br=cos(br);
		cr=	+ io[ix][iw].r
			- io[ix-1][iw].i * a0i
			- io[ix][iw].i   * a1i;
		ci=	+ io[ix][iw].i 
			+ io[ix-1][iw].r * a0i
			+ io[ix][iw].r   * a1i;
		buf[ix].r = cr*br - ci*bi;	
		buf[ix].i = cr*bi + ci*br;	
		
		for(ix=0;ix<nx;ix++)
		{
			io[ix][iw].r = buf[ix].r;
			io[ix][iw].i = buf[ix].i;
		}
	}
	for(ix=0;ix<nx;ix++)
	{
		io[ix][0].r = 0.0;
		io[ix][0].i = 0.0;
	}
#endif
}

void sf_psefd_release()
/*< release allocated memory >*/
{
	free(buf);
}


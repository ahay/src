#include <rsf.h>

void comaz_velocity_slice(int nx, int ny, float **slice);
/*< read a velocity slice [ny][nx] >*/

void comaz_image_slice (int nx, int ny, float **image);
/*< output an image slice [ny][nx] >*/

void comaz_gather_slice (int nh, int nx, int ny, float ***gather);
/*< output a gather slice [ny][nx][nh] >*/

void comaz(float w              /* frequency */,
	   int nz               /* depth samples */,
	   int nh               /* offset samples */, 
	   int nx               /* midpoint inline samples */,
	   int ny               /* midpoint cross-line samples */,
	   float dz             /* depth sampling */,
	   float dh             /* offset sampling */,
	   float h0             /* offset origin */,
	   float dx             /* midpoint inline sampling */,
	   float dy             /* midpoint cross-line sampling */,
	   int jx               /* inline subsampling for gathers */,
	   int jy               /* cross-line subsampling for gathers */,
	   float complex ***data      /* frequency slice [ny][nx][nh] */)
/*< common-azimuth migration >*/
{
    int nkx,nky,iz,ix,iy;
    float kx,dkx,fkx,ky,dky,fky;
    float complex cshift,w2;
    kiss_fft_cfg xforw, xinvs, yforw, yinvs;
    
        /* determine wavenumber sampling, pad by 2 */
    nkx = nx*2;
    nkx *= 2;
    dkx = 2.0*SF_PI/(nkx*dx);
    fkx = -SF_PI/dx;

    nky = ny*2;
    nky *= 2;
    dky = 2.0*SF_PI/(nky*dy);
    fky = -SF_PI/dy;

    xforw = kiss_fft_alloc(nkx,0,NULL,NULL);
    xinvs = kiss_fft_alloc(nkx,1,NULL,NULL);
    yforw = kiss_fft_alloc(nky,0,NULL,NULL);
    yinvs = kiss_fft_alloc(nky,1,NULL,NULL);

    if (NULL == xforw || NULL == xinvs || NULL == yforw || NULL == yinvs) 
	sf_error("%s: KISS FFT allocation error",__FILE__);

    /* loop over migrated depths z */
    for (iz=0; iz<nz-1; iz++) {
	/* accumulate image (summed over frequency and offset) */
    }
}

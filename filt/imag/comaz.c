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
    int iz;

    /* loop over migrated depths z */
    for (iz=0; iz<nz-1; iz++) {
	/* accumulate image (summed over frequency and offset) */
    }
}

/* Two-step Lowrank - 2-D as a Linear operator*/

#include <math.h>
#include <rsf.h>
#include "lr2zortm.h"
#include "fft2.h"


static int nt, nx, nz, nk, nx2, nz2, ntx, nzx, nzx2, pad1, m2;
static float *curr, *prev, *currm;
static float **rht, **lft,**wave2;
static sf_complex cc,*cwave, *cwavem,**wave;

void lr2zortm_init(int nt_in /*trace lenhth */, 
		   int nx_in /* number of traces */,
		   int nz_in /* depth length */,
		   int nk_in /* FFT out size */,
		   int nxpad,
		   int nzpad,
		   int m2_in,
		   float **left,
		   float **right)
/*< initialize >*/
{
    nt = nt_in;
    nx = nx_in;
    nz = nz_in;
    nk = nk_in;
    nz2 = nzpad;
    nx2 = nxpad;

    ntx = nt*nx;
    nzx = nz*nx;
    nzx2 = nz2*nx2;

    m2 = m2_in;

    pad1 = 1;

    lft = left;
    rht = right;

    /* Alloc arrays */
    curr = sf_floatalloc(nzx2);
    currm = sf_floatalloc(nzx2);    
    prev = sf_floatalloc(nzx);

    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);

    wave2  = sf_floatalloc2(nzx2,m2);
    wave = sf_complexalloc2(nk,m2);

}




void lr2zortm_close(void)
/*< clean allocated storage >*/

{
    free(*wave2);
    free(*wave);
    free(wave2);
    free(wave);
    free(*lft);
    free(*rht);
    free(lft);
    free(rht);
    free(curr);
    free(currm);
    free(prev);
    free(cwave);
    free(cwavem);

}




void lr2zortm_lop(bool adj, bool add, int nm, int nd, float *modl, float *data)
/*< linear operator >*/
{


    int it,iz,ix,ik,im,i,j;
    float old;
    float **dat, **img;
    sf_adjnull(adj,add,nm,nd,modl,data);


    /* initialization */
    img = sf_floatalloc2(nz,nx);
    dat = sf_floatalloc2(nt,nx);
 

    for (ik = 0; ik < nk; ik++) {
	for (im = 0; im < m2; im++) {
	    wave[im][ik] = sf_cmplx(0.,0.);
	}
    }

    for (iz=0; iz < nzx; iz++) {
	prev[iz]=0.;
    }

    for (iz=0; iz < nzx2; iz++) {
	curr[iz]=0.;
	currm[iz]=0.;
    }


    /* Main program */

    if(adj){ /* migration */

	/* scheme: p(t-dt) = 2p(t) - p(t+dt) + A'p(t) + src(t) */

	float *c; /* c is used to store 2*p(t) - p(t+dt) */
	sf_complex *ctmp;

	c = sf_floatalloc(nzx2);
	ctmp = sf_complexalloc(nk);

	for(ix=0;ix<nx;ix++){
	    for(it=0;it<nt;it++)
	    {
		dat[ix][it] = data[it+ix*nt];
	    }
	}

	/* for forward ffts of each rank column */
	fft2_allocate(ctmp);
	/* for inverse fft */
	ifft2_allocate(cwave);

	/* time stepping */
	for (it=nt-1;it>-1; it--) {	
	    sf_warning("Remigrate it=%d;",it);

	    /* update P(t) and P(t+dt) */

	    for (ix = 0; ix < nx; ix++) {
		for (iz=0; iz < nz; iz++) {
		    i = ix+iz*nx;  /* original grid */
		    j = ix+iz*nx2; /* padded grid */
			
		    c[j] = curr[j];
		    old = curr[j];
		    c[j] += c[j] - prev[i];
		    prev[i] = old;

		}
	    }

	    /* matrix multiplication */
	    for (im = 0; im < m2; im++) {

		for (ix = 0; ix < nx; ix++) {
		    for (iz=0; iz < nz; iz++) {
			i = ix+iz*nx;  /* original grid */
			j = ix+iz*nx2; /* padded grid */

			currm[j] = lft[im][i]*curr[j];
		    }
		}
		fft2(currm,ctmp);
		for (ik=0;ik<nk;ik++){
		    /* copy data to wave*/
		    wave[im][ik] =  ctmp[ik];

		}
	    }

	    for (ik = 0; ik < nk; ik++) {
		cc = sf_cmplx(0.,0.);
		for (im = 0; im < m2; im++) {
		    cc += wave[im][ik]*rht[ik][im];
		}	 
		cwave[ik] = cc;
	    }

	    ifft2(curr,cwave);
		

	    for (ix = 0; ix < nx; ix++) {
		for (iz=0; iz < nz; iz++) {
		    i = ix+iz*nx;  /* original grid */
		    j = ix+iz*nx2; /* padded grid */

		    /* p(t-dt) = 2p(t) - p(t+dt) + A'p(t)  */
		    curr[j] += c[j];
		}
	    }

	    /* inject data */
	    for (ix=0; ix < nx; ix++) {
		curr[ix] += dat[ix][it];
	    }

	}/* End of time stepping */

	/* collect imag at time 0 */

	for (ix=0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		img[ix][iz] = curr[ix+iz*nx2];
		/* add to existing modl based on adjnull */
		modl[iz+ix*nz] += img[ix][iz];

	    }
	}



    } else{ /* Forward */

	float c;


	for(ix=0;ix<nx;ix++){
	    for(iz=0;iz<nz;iz++)
	    {
		img[ix][iz] = modl[iz+ix*nz];

	    }
	}

	/* transpose & initialize exploding refl*/
	for (ix=0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		curr[ix+iz*nx2]=img[ix][iz];
	    }
	}
	/* Initialize recorded data */

	for (ix=0; ix < nx; ix++) {
	    for (it=0; it < nt; it++) {
		dat[ix][it] = 0.0f;
	    }
	}

	/* time stepping */

	/* Alloc for forward (cwave) only */
	fft2_allocate(cwave);
	/* Alloc for inverse (cwavem) only */
	ifft2_allocate(cwavem);

	for (it=0; it<nt; it++) {
	    sf_warning("Modeling it=%d;",it);

	    /* record data on the surface */
	    for (ix=0; ix < nx; ix++) {
		dat[ix][it] = curr[ix];
	    }

	    /* matrix multiplication */


	    fft2(curr,cwave);

	    /* I use separated variables to do fft / ifft */

	    for (im = 0; im < m2; im++) {
		for (ik = 0; ik < nk; ik++) {

		    cwavem[ik] = cwave[ik]*rht[ik][im];

		}

		ifft2(wave2[im],cwavem);
	    }

	    for (ix = 0; ix < nx; ix++) {
		for (iz=0; iz < nz; iz++) {
		    i = ix+iz*nx;  /* original grid */
		    j = ix+iz*nx2; /* padded grid */
			
		    old = c = curr[j];
		    c += c - prev[i];
		    prev[i] = old;

		    for (im = 0; im < m2; im++) {
			c += lft[im][i]*wave2[im][j];		    
		    }
		    curr[j] = c;
		}
	    }


	}/* End of time stepping */
	for(ix=0;ix<nx;ix++){
	    for(it=0;it<nt;it++)
	    {
		data[it+ix*nt] += dat[ix][it];
	    }
	}
    }

}

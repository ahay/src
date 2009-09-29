#include <rsf.h>
#include "slice.h"
#include "phs.h"
#include "nos.h"

#define LOOP(a) for(ihx=0;ihx<ahx.n;ihx++){ \
	        for(imx=0;imx<amx.n;imx++){ {a} }}

int main (int argc, char *argv[])
{
    bool freq;
    axa az,amx,ahx,aw;
    int iz,imx,ihx,iw;
    float *sz;
    float         **qq;
    float complex **wx;
    float            w;

    sf_file Fs,Fd,Fi;
    fslice slow,data,imag;

    sf_init(argc,argv);
    if (!sf_getbool("freq",&freq)) freq = true;

/*------------------------------------------------------------*/
    Fs = sf_input ("slo");
    iaxa(Fs,&az ,1);  az.l= "z";

    slow = fslice_init(az.n,1,sizeof(float));
    fslice_load(Fs,slow,SF_FLOAT);

    sz=sf_floatalloc(az.n);
    fslice_get(slow,0,sz);
    for (iz=0; iz<az.n; iz++) {
	sf_warning("sz[%d]=%g",iz,sz[iz]);
    }
/*------------------------------------------------------------*/
    Fd = sf_input ( "in");
    Fi = sf_output("out"); sf_settype(Fi,SF_FLOAT);
    if (SF_COMPLEX !=sf_gettype(Fd)) sf_error("Need complex data");
    
    iaxa(Fd,&amx,1); amx.l="mx"; oaxa(Fi,&amx,1);
    iaxa(Fd,&ahx,2); ahx.l="hx"; oaxa(Fi,&ahx,2);
    iaxa(Fd,&aw ,3);  aw.l= "w"; oaxa(Fi,&az ,3);
    
    /* from hertz to radian */
    aw.d *= 2.*SF_PI; 
    aw.o *= 2.*SF_PI;

    data = fslice_init(amx.n*ahx.n, aw.n,sizeof(float complex));
    imag = fslice_init(amx.n*ahx.n, az.n,sizeof(float));

    wx = sf_complexalloc2(amx.n,ahx.n);
    qq = sf_floatalloc2  (amx.n,ahx.n);

    fslice_load(Fd,data,SF_COMPLEX);
/*------------------------------------------------------------*/

    if(freq) {
	phs_init(amx,ahx,aw,az);
    } else {
	nos_init(amx,ahx,aw,az);
    }

/*------------------------------------------------------------*/
    LOOP( qq[ihx][imx] = 0; );
    for (iz=0; iz<az.n; iz++) {
	fslice_put(imag,iz,qq[0]);
    }

    /* w loop */
    for (iw=0; iw<aw.n; iw++) {
	w = aw.o+iw*aw.d;

	fslice_get(data,iw,wx[0]);

	fslice_get(imag,0,qq[0]);
	LOOP(;      qq[ihx][imx] += 
	     crealf(wx[ihx][imx] ); );
	fslice_put(imag,0,qq[0]);

	/* z loop */
	for (iz=0; iz<az.n-1; iz++) {
	    if(freq) {
		sf_warning ("PHS   iw=%3d of %3d:   iz=%3d of %3d",iw+1,aw.n,iz+1,az.n-1);
		/* k-domain phase shift */
		phs(w,sz[iz],wx);
	    } else {
		sf_warning ("NOS   iw=%3d of %3d:   iz=%3d of %3d",iw+1,aw.n,iz+1,az.n-1);
		nos(w,sz[iz],wx);
	    }

	    fslice_get(imag,iz+1,qq[0]);
	    LOOP(;      qq[ihx][imx] += 
		 crealf(wx[ihx][imx] ); );
	    fslice_put(imag,iz+1,qq[0]);
	}	
    }
    
/*------------------------------------------------------------*/
    fslice_dump(Fi,imag,SF_FLOAT);

    fslice_close(data);
    fslice_close(imag);
    fslice_close(slow);

    free( *qq); free( qq);
    free( *wx); free( wx);

    exit (0);
}

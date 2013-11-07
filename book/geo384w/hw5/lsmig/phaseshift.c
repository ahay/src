/* Phase-shift modeling and migration */

#include <rsf.h>

int main(int argc, char* argv[])
{
    bool adj;
    int ik, iw, nk, nw, iz, nz;
    float k, dk, k0, dw, dz, z0, *v, eps;
    sf_complex *dat, *mod, w2, ps, wave;
    sf_file inp, out, vel;

    sf_init(argc,argv);

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag, 0: modeling, 1: migration */

    inp = sf_input("in");
    out = sf_output("out");
    vel = sf_input("vel"); /* velocity file */

    if (!sf_histint(vel,"n1",&nz)) 
	sf_error("No n1= in vel");
    if (!sf_histfloat(vel,"d1",&dz)) 
	sf_error("No d1= in vel");
    if (!sf_histfloat(vel,"o1",&z0)) z0=0.0;
    
    if (!sf_histint(inp,"n2",&nk)) 
	sf_error("No n2= in input");
    if (!sf_histfloat(inp,"d2",&dk)) 
	sf_error("No d2= in input");
    if (!sf_histfloat(inp,"o2",&k0)) 
	sf_error("No o2= in input");

    dk *= 2*SF_PI;
    k0 *= 2*SF_PI;

    if (adj) { /* migration */
	if (!sf_histint(inp,"n1",&nw)) 
	    sf_error("No n1= in input");
	if (!sf_histfloat(inp,"d1",&dw)) 
	    sf_error("No d1= in input");

	sf_putint(out,"n1",nz);
	sf_putfloat(out,"d1",dz);
	sf_putfloat(out,"o1",z0);
    } else {  /* modeling */
	if (!sf_getint("nw",&nw)) sf_error("No nw=");
	if (!sf_getfloat("dw",&dw)) sf_error("No dw=");

	sf_putint(out,"n1",nw);
	sf_putfloat(out,"d1",dw);
	sf_putfloat(out,"o1",0.0);	
    }

    if (!sf_getfloat("eps",&eps)) eps = 1.0f;
    /* stabilization parameter */

    dw *= 2*SF_PI;

    dat = sf_complexalloc(nw);
    mod = sf_complexalloc(nz);
    
    /* read velocity, convert to slowness squared */
    v = sf_floatalloc(nz);
    sf_floatread(v,nz,vel);
    for (iz=0; iz < nz; iz++) { 
	v[iz] = 2.0f/v[iz];
	v[iz] *= v[iz];
    }

    for (ik=0; ik<nk; ik++) { /* wavenumber */
	sf_warning("wavenumber %d of %d;",ik+1,nk);
	k=k0+ik*dk;
	k *= k;
	
	if (adj) {
	    sf_complexread(dat,nw,inp);
	    for (iz=0; iz < nz; iz++) { 
		mod[iz]=sf_cmplx(0.0,0.0);
	    }
	} else {
	    sf_complexread(mod,nz,inp);
	}

	for (iw=0; iw<nw; iw++) { /* frequency */
	    w2 = sf_cmplx(eps*dw,iw*dw);
	    w2 *= w2;

	    if (adj) { /* migration */
		wave = dat[iw];
		for (iz=0; iz < nz; iz++) { 
		    /* !!! FILL MISSING LINES !!! */
		}
	    } else {  /* modeling */
		wave = mod[nz-1];
		for (iz=nz-2; iz >=0; iz--) {
		    ps = cexpf(-csqrt(w2*v[iz]+k)*dz);
		    wave = wave*ps+mod[iz];
		}
		dat[iw] = wave;
	    }
	}
	
	if (adj) {
	    sf_complexwrite(mod,nz,out);
	} else {
	    sf_complexwrite(dat,nw,out);
	}
    } 
    sf_warning(".");

    exit(0);
}

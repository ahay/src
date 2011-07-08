#include <rsf.h>
#include <omp.h>

#define nx 200
#define nz 100
#define k1 2*SF_PI*4

int main(int argc, char *argv[])
{
#ifdef _OPENMP
    int nth;
#pragma omp parallel
    nth=omp_get_num_threads();
    omp_set_num_threads(nth);
    sf_warning("using total %d threads",nth);
#endif

    int N,nk,ik,ix,iz;
    float dx,*x,dk,*k,**u=NULL,**re=NULL,**im=NULL;
    sf_complex **ut=NULL;
    sf_file Fout,Fre,Fim,Fk;
    kiss_fftr_cfg *cfg;

    sf_init(argc,argv);
    Fout=sf_output("out");  /* u(x) */
    Fre=sf_output("real");  /* u(k)=FT(u) */
    Fim=sf_output("imag");
    Fk =sf_output("k");     /* k */

    N=2*kiss_fft_next_fast_size((nx+1)/2);
    nk=N/2+1;
    sf_warning("nfft=%d",N);

    u =sf_floatalloc2(nx,nz);
    ut=sf_complexalloc2(nk,nz);
    x =sf_floatalloc(nx);
    k =sf_floatalloc(nk);
    re=sf_floatalloc2(nk,nz);
    im=sf_floatalloc2(nk,nz);

    cfg = (kiss_fftr_cfg *) sf_alloc(nz,sizeof(*cfg));
    for (iz=0; iz<nz; iz++)
	cfg[iz] =kiss_fftr_alloc(N,0,NULL,NULL);

    for (dx=1.0f/nx,ix=0; ix<nx; ix++) x[ix]=ix*dx;
    for (iz=0; iz<nz; iz++)
	for (ix=0; ix<nx; ix++)
	    u[iz][ix]=sin(k1*x[ix]);

    for (dk=1.0f/N,ik=0; ik<nk; ik++) k[ik]=ik*dk;



#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)\
    private(iz) shared(u,ut,cfg)
#endif
    for (iz=0; iz<nz; iz++) 
	kiss_fftr(cfg[iz],u[iz],(kiss_fft_cpx*)ut[iz]);

    for (iz=0; iz<nz; iz++)
	for(ik=0; ik<nk; ik++) {
	    re[iz][ik]=creal(ut[iz][ik]);
	    im[iz][ik]=cimag(ut[iz][ik]);
	}

    sf_putint(Fout,"n1",nx);
    sf_putfloat(Fout,"o1",0.0f);
    sf_putfloat(Fout,"d1",dx);
    sf_putint(Fout,"n2",nz);
    for (iz=0; iz<nz; iz++)
	sf_floatwrite(u[iz],nx,Fout);

    sf_putint(Fk,"n1",nk);
    sf_floatwrite(k,nk,Fk);

    sf_putint(Fre,"n1",nk);
    sf_putint(Fre,"n2",nz);
    sf_putint(Fim,"n1",nk);
    sf_putint(Fim,"n2",nz);  
    for (iz=0; iz<nz; iz++) {
	sf_floatwrite(re[iz],nk,Fre);
	sf_floatwrite(im[iz],nk,Fim);
    }

    exit(0);
}

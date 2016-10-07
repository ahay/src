/* Multi-threaded multi-dimensional fft interface based on Kiss-FFT (both 2D and 3D, both r2c and c2c) */ 
/* 
   n1 <-> nz, n2 <-> nx, n3 <-> ny.
   3D FFT -> n3>1, n2>1, n1>1;
   2D FFT -> n3=1; n2>1; n1>1;
   */

#include <rsf.h>

#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef SF_HAS_FFTW
#include <fftw3.h>
#endif

static bool r2c;
static bool pio;
static int nz, nx, ny;
static int n1, n2, n3;
static float wt;

static sf_complex *cc=NULL;
#ifdef SF_HAS_FFTW
static sf_complex *dd=NULL;
static fftwf_plan cfg=NULL, icfg=NULL;
#else
static kiss_fft_cfg *cfg1=NULL, *icfg1=NULL, *cfg2=NULL, *icfg2=NULL, *cfg3=NULL, *icfg3=NULL;
static kiss_fft_cpx *tmp=NULL, **ctrace2=NULL, **ctrace3=NULL;
static sf_complex **trace2=NULL, **trace3=NULL;
#endif

int fft_init(int pad1 /* padding on the first axis */,
        int nz_,   int nx_,   int ny_ /* input data size */, 
        int *nz2, int *nx2, int *ny2 /* padded data size */,
        bool rtoc /* real to complex fft flag */,
        bool padio /* inp and out are padded*/)
    /*< initialize >*/
{
#ifdef SF_HAS_FFTW
#ifdef _OPENMP
    fftwf_init_threads();
    sf_warning("Using threaded FFTW3! %d\n",omp_get_max_threads());
    fftwf_plan_with_nthreads(omp_get_max_threads());
#endif
#else
    int i, nth=1;
#endif

    /* real to complex fft flag */
    r2c = rtoc;
    pio = padio;
    nz = nz_;
    nx = nx_;
    ny = ny_;

    if (nx==1 && ny>1) sf_error("For 2D, second axis should be nx!");

    /* axis 1 */
    n1 = kiss_fft_next_fast_size(nz*pad1);
    /* axis 2 */
    n2 = kiss_fft_next_fast_size(nx);
    /* axis 3 */
    if (ny>1) n3 = kiss_fft_next_fast_size(ny);
    else n3 = 1;

    cc = sf_complexalloc(n1*n2*n3);

#ifdef SF_HAS_FFTW

    /* fftw3 initialization */
    dd = sf_complexalloc(n1*n2*n3);

    if (ny>1) {
        cfg = fftwf_plan_dft_3d(n3,n2,n1,
                (fftwf_complex *) cc, 
                (fftwf_complex *) dd,
                FFTW_FORWARD, FFTW_MEASURE);

        icfg = fftwf_plan_dft_3d(n3,n2,n1,
                (fftwf_complex *) dd, 
                (fftwf_complex *) cc,
                FFTW_BACKWARD, FFTW_MEASURE);
    } else {
        cfg = fftwf_plan_dft_2d(n2,n1,
                (fftwf_complex *) cc, 
                (fftwf_complex *) dd,
                FFTW_FORWARD, FFTW_MEASURE);

        icfg = fftwf_plan_dft_2d(n2,n1,
                (fftwf_complex *) dd, 
                (fftwf_complex *) cc,
                FFTW_BACKWARD, FFTW_MEASURE);
    }

    if (NULL == cfg || NULL == icfg) sf_error("FFTW failure.");

#else

    /* kiss-fft initialization */
#ifdef _OPENMP
#pragma omp parallel
    {nth = omp_get_num_threads();}
#endif

    cfg1  = (kiss_fft_cfg *) sf_alloc(nth,sizeof(kiss_fft_cfg));
    icfg1 = (kiss_fft_cfg *) sf_alloc(nth,sizeof(kiss_fft_cfg));
    cfg2  = (kiss_fft_cfg *) sf_alloc(nth,sizeof(kiss_fft_cfg));
    icfg2 = (kiss_fft_cfg *) sf_alloc(nth,sizeof(kiss_fft_cfg));
    if (n3>1) {
        cfg3  = (kiss_fft_cfg *) sf_alloc(nth,sizeof(kiss_fft_cfg));
        icfg3 = (kiss_fft_cfg *) sf_alloc(nth,sizeof(kiss_fft_cfg));
    } else {
        cfg3  = NULL;
        icfg3 = NULL;
    }

    for (i=0; i<nth; i++) {
        cfg1[i] = kiss_fft_alloc(n1,0,NULL,NULL);
        icfg1[i]= kiss_fft_alloc(n1,1,NULL,NULL);
        cfg2[i] = kiss_fft_alloc(n2,0,NULL,NULL);
        icfg2[i]= kiss_fft_alloc(n2,1,NULL,NULL);
        if (n3>1) {
            cfg3[i] = kiss_fft_alloc(n3,0,NULL,NULL);
            icfg3[i]= kiss_fft_alloc(n3,1,NULL,NULL);
        }
    }

    trace2 = sf_complexalloc2(n2,nth);
    ctrace2= (kiss_fft_cpx **) trace2;
    if (n3>1) {
        trace3 = sf_complexalloc2(n3,nth);
        ctrace3= (kiss_fft_cpx **) trace3;
    } else {
        trace3 = NULL;
        ctrace3 = NULL;
    }

    tmp = (kiss_fft_cpx *) sf_alloc(n3*n2*n1,sizeof(kiss_fft_cpx));

#endif

    *nz2 = n1;
    *nx2 = n2;
    *ny2 = n3;

    wt =  1.0/(n3*n2*n1);

    return (n3*n2*n1);
}

void fft(void *inp /* [n1*n2*n3] */, 
        sf_complex *out /* [n1*n2*n3] */)
    /*< 3-D FFT >*/
{
#ifndef SF_HAS_FFTW
    int ith=0;
#endif
    int i1, i2, i3;
    float *inpf;
    sf_complex *inpc;

    if (r2c) {
        inpf = (float *) inp;
        inpc = NULL;
    } else {
        inpf = NULL;
        inpc = (sf_complex *) inp;
    }

    if (n3>1) {

        if (pio) {
            /* FFT centering */
#ifdef _OPENMP
#pragma omp parallel for private(i3,i2,i1) default(shared)
#endif
            for (i3=0; i3<n3; i3++) {
                for (i2=0; i2<n2; i2++) {
                    for (i1=0; i1<n1; i1++) {
                        if (r2c) {
                            cc[(i3*n2+i2)*n1+i1] = sf_cmplx((((i3%2==0)==(i2%2==0))==(i1%2==0))? inpf[(i3*n2+i2)*n1+i1]:-inpf[(i3*n2+i2)*n1+i1],0.);
                        } else {
#ifdef SF_HAS_COMPLEX_H
                            cc[(i3*n2+i2)*n1+i1] = (((i3%2==0)==(i2%2==0))==(i1%2==0))? inpc[(i3*n2+i2)*n1+i1]:-inpc[(i3*n2+i2)*n1+i1];
#else
                            cc[(i3*n2+i2)*n1+i1] = (((i3%2==0)==(i2%2==0))==(i1%2==0))? inpc[(i3*n2+i2)*n1+i1]:sf_cneg(inpc[(i3*n2+i2)*n1+i1]);
#endif
                        }
                    }
                }
            }
        } else {
            /* FFT centering */
#ifdef _OPENMP
#pragma omp parallel for private(i3,i2,i1) default(shared)
#endif
            for (i3=0; i3<n3; i3++) {
                for (i2=0; i2<n2; i2++) {
                    for (i1=0; i1<n1; i1++) {
                        if (r2c) {
                            if (i1<nz && i2<nx && i3<ny) {
                                cc[(i3*n2+i2)*n1+i1] = sf_cmplx((((i3%2==0)==(i2%2==0))==(i1%2==0))? inpf[(i3*nx+i2)*nz+i1]:-inpf[(i3*nx+i2)*nz+i1],0.);
                            } else {
                                cc[(i3*n2+i2)*n1+i1] = sf_cmplx(0.,0.); 
                            }
                        } else {
                            if (i1<nz && i2<nx && i3<ny) {
#ifdef SF_HAS_COMPLEX_H
                                cc[(i3*n2+i2)*n1+i1] = (((i3%2==0)==(i2%2==0))==(i1%2==0))? inpc[(i3*nx+i2)*nz+i1]:-inpc[(i3*nx+i2)*nz+i1];
#else
                                cc[(i3*n2+i2)*n1+i1] = (((i3%2==0)==(i2%2==0))==(i1%2==0))? inpc[(i3*nx+i2)*nz+i1]:sf_cneg(inpc[(i3*nx+i2)*nz+i1]);
#endif
                            } else {
                                cc[(i3*n2+i2)*n1+i1] = sf_cmplx(0.,0.); 
                            }
                        }
                    }
                }
            }
        }

#ifdef SF_HAS_FFTW

    fftwf_execute(cfg);

#ifdef _OPENMP
#pragma omp parallel for private(i1) default(shared)
#endif
    for (i1=0; i1<n3*n2*n1; i1++) {
        out[i1]=dd[i1];
    }

#else

        /* FFT over first axis */
#ifdef _OPENMP
#pragma omp parallel for private(i3,i2,i1,ith) default(shared)
#endif
        for (i3=0; i3<n3; i3++) {
#ifdef _OPENMP
            ith = omp_get_thread_num();
#endif
            for (i2=0; i2<n2; i2++) {
                kiss_fft_stride(cfg1[ith],(kiss_fft_cpx *) (cc+(i3*n2+i2)*n1),tmp+(i3*n2+i2)*n1,1);
            }
        }

        /* FFT over second axis */
#ifdef _OPENMP
#pragma omp parallel for private(i3,i2,i1,ith) default(shared)
#endif
        for (i3=0; i3<n3; i3++) {
#ifdef _OPENMP
            ith = omp_get_thread_num();
#endif
            for (i1=0; i1<n1; i1++) {
                kiss_fft_stride(cfg2[ith],tmp+i3*n2*n1+i1,ctrace2[ith],n1);
                for (i2=0; i2<n2; i2++) {
                    tmp[(i3*n2+i2)*n1+i1]=ctrace2[ith][i2];
                }
            }
        }

        /* FFT over third axis */
#ifdef _OPENMP
#pragma omp parallel for private(i3,i2,i1,ith) default(shared)
#endif
        for (i2=0; i2<n2; i2++) {
#ifdef _OPENMP
            ith = omp_get_thread_num();
#endif
            for (i1=0; i1<n1; i1++) {
                kiss_fft_stride(cfg3[ith],tmp+i2*n1+i1,ctrace3[ith],n2*n1);
                for (i3=0; i3<n3; i3++) {
                    out[(i3*n2+i2)*n1+i1] = trace3[ith][i3];
                }
            }
        }

#endif

    } else {

        if (pio) {
        /* FFT centering */
#ifdef _OPENMP
#pragma omp parallel for private(i2,i1) default(shared)
#endif
            for (i2=0; i2<n2; i2++) {
                for (i1=0; i1<n1; i1++) {
                    if (r2c) {
                        cc[i2*n1+i1] = sf_cmplx(((i2%2==0)==(i1%2==0))? inpf[i2*n1+i1]:-inpf[i2*n1+i1],0.);
                    } else {
#ifdef SF_HAS_COMPLEX_H
                        cc[i2*n1+i1] = ((i2%2==0)==(i1%2==0))? inpc[i2*n1+i1]:-inpc[i2*n1+i1];
#else
                        cc[i2*n1+i1] = ((i2%2==0)==(i1%2==0))? inpc[i2*n1+i1]:sf_cneg(inpc[i2*n1+i1]);
#endif
                    }
                }
            }
        } else {
#ifdef _OPENMP
#pragma omp parallel for private(i2,i1) default(shared)
#endif
            for (i2=0; i2<n2; i2++) {
                for (i1=0; i1<n1; i1++) {
                    if (r2c) {
                        if (i1<nz && i2<nx) {
                            cc[i2*n1+i1] = sf_cmplx(((i2%2==0)==(i1%2==0))? inpf[i2*nz+i1]:-inpf[i2*nz+i1],0.);
                        } else {
                            cc[i2*n1+i1] = sf_cmplx(0.,0.);
                        }
                    } else {
                        if (i1<nz && i2<nx) {
#ifdef SF_HAS_COMPLEX_H
                            cc[i2*n1+i1] = ((i2%2==0)==(i1%2==0))? inpc[i2*nz+i1]:-inpc[i2*nz+i1];
#else
                            cc[i2*n1+i1] = ((i2%2==0)==(i1%2==0))? inpc[i2*nz+i1]:sf_cneg(inpc[i2*nz+i1]);
#endif
                        } else {
                            cc[i2*n1+i1] = sf_cmplx(0.,0.);
                        }
                    }
                }
            }
        }

#ifdef SF_HAS_FFTW

    fftwf_execute(cfg);

#ifdef _OPENMP
#pragma omp parallel for private(i1) default(shared)
#endif
    for (i1=0; i1<n3*n2*n1; i1++) {
        out[i1]=dd[i1];
    }

#else

        /* FFT over first axis */
#ifdef _OPENMP
#pragma omp parallel for private(i2,i1,ith) default(shared)
#endif
        for (i2=0; i2<n2; i2++) {
#ifdef _OPENMP
            ith = omp_get_thread_num();
#endif
            kiss_fft_stride(cfg1[ith],(kiss_fft_cpx *) (cc+i2*n1),tmp+i2*n1,1);
        }

        /* FFT over second axis */
#ifdef _OPENMP
#pragma omp parallel for private(i2,i1,ith) default(shared)
#endif
        for (i1=0; i1<n1; i1++) {
#ifdef _OPENMP
            ith = omp_get_thread_num();
#endif
            kiss_fft_stride(cfg2[ith],tmp+i1,ctrace2[ith],n1);
            for (i2=0; i2<n2; i2++) {
                out[i2*n1+i1] = trace2[ith][i2];
            }
        }

#endif

    } /*if (n3>1)*/

}

void ifft(void *out /* [n1*n2*n3] */, 
        sf_complex *inp /* [n1*n2*n3] */)
    /*< 3-D inverse FFT >*/
{
#ifndef SF_HAS_FFTW
    int ith=0;
#endif
    int i1, i2, i3;
    float *outf;
    sf_complex *outc;

    if (r2c) {
        outf = (float *) out;
        outc = NULL;
    } else {
        outf = NULL;
        outc = (sf_complex *) out;
    }

    if (n3>1) {

#ifdef SF_HAS_FFTW

#ifdef _OPENMP
#pragma omp parallel for private(i1) default(shared)
#endif
    for (i1=0; i1<n3*n2*n1; i1++) {
        dd[i1]=inp[i1];
    }

    fftwf_execute(icfg);

#else

        /* IFFT over third axis */
#ifdef _OPENMP
#pragma omp parallel for private(i3,i2,i1,ith) default(shared)
#endif
        for (i2=0; i2<n2; i2++) {
#ifdef _OPENMP
            ith = omp_get_thread_num();
#endif
            for (i1=0; i1<n1; i1++) {
                kiss_fft_stride(icfg3[ith],(kiss_fft_cpx *) (inp+i2*n1+i1),ctrace3[ith],n2*n1);
                for (i3=0; i3<n3; i3++) {
                    tmp[(i3*n2+i2)*n1+i1] = ctrace3[ith][i3];
                }
            }
        }

        /* IFFT over second axis */
#ifdef _OPENMP
#pragma omp parallel for private(i3,i2,i1,ith) default(shared)
#endif
        for (i3=0; i3<n3; i3++) {
#ifdef _OPENMP
            ith = omp_get_thread_num();
#endif
            for (i1=0; i1<n1; i1++) {
                kiss_fft_stride(icfg2[ith],tmp+i3*n2*n1+i1,ctrace2[ith],n1);
                for (i2=0; i2<n2; i2++) {
                    tmp[(i3*n2+i2)*n1+i1]=ctrace2[ith][i2];
                }
            }
        }

        /* IFFT over first axis */
#ifdef _OPENMP
#pragma omp parallel for private(i3,i2,i1,ith) default(shared)
#endif
        for (i3=0; i3<n3; i3++) {
#ifdef _OPENMP
            ith = omp_get_thread_num();
#endif
            for (i2=0; i2<n2; i2++) {
                kiss_fft_stride(icfg1[ith],tmp+(i3*n2+i2)*n1,(kiss_fft_cpx *) (cc+(i3*n2+i2)*n1),1);
            }
        }

#endif

        if (pio) {
            /* FFT centering and normalization*/
#ifdef _OPENMP
#pragma omp parallel for private(i3,i2,i1) default(shared)
#endif
            for (i3=0; i3<n3; i3++) {
                for (i2=0; i2<n2; i2++) {
                    for (i1=0; i1<n1; i1++) {
                        if (r2c) {
                            outf[(i3*n2+i2)*n1+i1] = ((((i3%2==0)==(i2%2==0))==(i1%2==0))? wt:-wt)*crealf(cc[(i3*n2+i2)*n1+i1]);
                        } else {
#ifdef SF_HAS_COMPLEX_H
                            outc[(i3*n2+i2)*n1+i1] = ((((i3%2==0)==(i2%2==0))==(i1%2==0))? wt:-wt)*cc[(i3*n2+i2)*n1+i1];
#else
                            outc[(i3*n2+i2)*n1+i1] = sf_crmul(cc[(i3*n2+i2)*n1+i1],((((i3%2==0)==(i2%2==0))==(i1%2==0))? wt:-wt));
#endif
                        }
                    }
                }
            }
        } else {
            /* FFT centering and normalization*/
#ifdef _OPENMP
#pragma omp parallel for private(i3,i2,i1) default(shared)
#endif
            for (i3=0; i3<ny; i3++) {
                for (i2=0; i2<nx; i2++) {
                    for (i1=0; i1<nz; i1++) {
                        if (r2c) {
                            outf[(i3*nx+i2)*nz+i1] = ((((i3%2==0)==(i2%2==0))==(i1%2==0))? wt:-wt)*crealf(cc[(i3*n2+i2)*n1+i1]);
                        } else {
#ifdef SF_HAS_COMPLEX_H
                            outc[(i3*nx+i2)*nz+i1] = ((((i3%2==0)==(i2%2==0))==(i1%2==0))? wt:-wt)*cc[(i3*n2+i2)*n1+i1];
#else
                            outc[(i3*nx+i2)*nz+i1] = sf_crmul(cc[(i3*n2+i2)*n1+i1],((((i3%2==0)==(i2%2==0))==(i1%2==0))? wt:-wt));
#endif
                        }
                    }
                }
            }
        }

    } else {

#ifdef SF_HAS_FFTW

#ifdef _OPENMP
#pragma omp parallel for private(i1) default(shared)
#endif
    for (i1=0; i1<n3*n2*n1; i1++) {
        dd[i1]=inp[i1];
    }

    fftwf_execute(icfg);

#else

        /* IFFT over second axis */
#ifdef _OPENMP
#pragma omp parallel for private(i2,i1,ith) default(shared)
#endif
        for (i1=0; i1<n1; i1++) {
#ifdef _OPENMP
            ith = omp_get_thread_num();
#endif
            kiss_fft_stride(icfg2[ith],(kiss_fft_cpx *) (inp+i1),ctrace2[ith],n1);
            for (i2=0; i2<n2; i2++) {
                tmp[i2*n1+i1]=ctrace2[ith][i2];
            }
        }

        /* IFFT over first axis */
#ifdef _OPENMP
#pragma omp parallel for private(i2,i1,ith) default(shared)
#endif
        for (i2=0; i2<n2; i2++) {
#ifdef _OPENMP
            ith = omp_get_thread_num();
#endif
            kiss_fft_stride(icfg1[ith],tmp+i2*n1,(kiss_fft_cpx *) (cc+i2*n1),1);
        }

#endif

        if (pio) {
            /* FFT centering and normalization*/
#ifdef _OPENMP
#pragma omp parallel for private(i2,i1) default(shared)
#endif
            for (i2=0; i2<n2; i2++) {
                for (i1=0; i1<n1; i1++) {
                    if (r2c) {
                        outf[i2*n1+i1] = (((i2%2==0)==(i1%2==0))? wt:-wt)*crealf(cc[i2*n1+i1]);
                    } else {
#ifdef SF_HAS_COMPLEX_H
                        outc[i2*n1+i1] = (((i2%2==0)==(i1%2==0))? wt:-wt)*cc[i2*n1+i1];
#else
                        outc[i2*n1+i1] = sf_crmul(cc[i2*n1+i1],(((i2%2==0)==(i1%2==0))? wt:-wt));
#endif
                    }
                }
            }
        } else {
            /* FFT centering and normalization*/
#ifdef _OPENMP
#pragma omp parallel for private(i2,i1) default(shared)
#endif
            for (i2=0; i2<nx; i2++) {
                for (i1=0; i1<nz; i1++) {
                    if (r2c) {
                        outf[i2*nz+i1] = (((i2%2==0)==(i1%2==0))? wt:-wt)*crealf(cc[i2*n1+i1]);
                    } else {
#ifdef SF_HAS_COMPLEX_H
                        outc[i2*nz+i1] = (((i2%2==0)==(i1%2==0))? wt:-wt)*cc[i2*n1+i1];
#else
                        outc[i2*nz+i1] = sf_crmul(cc[i2*n1+i1],(((i2%2==0)==(i1%2==0))? wt:-wt));
#endif
                    }
                }
            }
        }

    } /*if (n3>1)*/

}

void fft_finalize()
    /*< clean up fft >*/
{
    /* make sure everything is back to its pristine state */
#ifndef SF_HAS_FFTW
    int i, nth=1;
#endif

#ifdef SF_HAS_FFTW
#ifdef _OPENMP
    fftwf_cleanup_threads();
#endif
    fftwf_destroy_plan(cfg);
    if (NULL!=cfg)   cfg=NULL;
    fftwf_destroy_plan(icfg);
    if (NULL!=icfg) icfg=NULL;
    fftwf_cleanup();
    if (NULL!=dd) { free(dd); dd=NULL; }
#else
#ifdef _OPENMP
#pragma omp parallel
    {nth = omp_get_num_threads();}
#endif
    for (i=0; i<nth; i++) {
        if (NULL!=cfg1)  { free(cfg1[i]);   cfg1[i]=NULL; }
        if (NULL!=icfg1) { free(icfg1[i]); icfg1[i]=NULL; }
        if (NULL!=cfg2)  { free(cfg2[i]);   cfg2[i]=NULL; }
        if (NULL!=icfg2) { free(icfg2[i]); icfg2[i]=NULL; }
        if (NULL!=cfg3)  { free(cfg3[i]);   cfg3[i]=NULL; }
        if (NULL!=icfg3) { free(icfg3[i]); icfg3[i]=NULL; }
    }
    if (NULL!=tmp) { free(tmp); tmp=NULL; }
    if (NULL!=trace2) { free(*trace2); free(trace2); trace2=NULL; }
    if (NULL!=trace3) { free(*trace3); free(trace3); trace3=NULL; }
#endif

    if (NULL!=cc) { free(cc); cc=NULL; }
}


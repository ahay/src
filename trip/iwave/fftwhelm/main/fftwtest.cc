/* Helmholtz */

#include <par.h>
#include <gridio.h>
#include <stdio.h>
#include <parser.h>
#include <fftw3.h>
#include <iostream>
using namespace std;
#include "rsf.h"

//#include "fft2.h"

//extern "C" {int fft2_init(bool cmplx1        /* if complex transform */,
//                         int pad1           /* padding on the first axis */,
//                         int nx,   int ny   /* input data size */,
//                          int *nx2, int *ny2 /* padded data size */);}

/* lenwork must be > 6*n1*n2+3*max(n2,2*n1)+21 */

/*extern void pagedoc(void);*/
/*extern void requestdoc(int);*/

const char * sdoc[] = {
    "IWAVE 2D Helmholtz powers via NCAR FFTs",
    "",
    "Usage: helm.x in= out= [power=0.0, datum=0.0, scale1=1.0, scale2=1.0]",
    "",
    "  applies arbitrary real power of elliptic Helmholtz operator",
    "",
    "  u(x_1,x_2) += [w_1 d^2/dx_1^2 + w_2 d^2/dx_2^2]u(x_1,x_2) ",
    "",
    "  to regular grid function in rectangle x_{1,min} <= x_1 <= x_{1,max},",
    "  x_{2,min} <= x_2 <= x_{2,max}. Imposes Dirichlet boundary conditions",
    "  on datum level x_1=max(x_{1,min},x_{1,min}+datum)",
    "  and on sides x_2=x_{2,min}, x_{2,max}, and Neumann condition on",
    "  bottom x_1=x_{1,max}. Input function u is not modified in the slab",
    "  x_{1,min} <= x_1 <= max(x_{1,min},x_{1,min}+datum).",
    "  Input, output files are RSF format header files. If it exists, output",
    "  header file must define same grid as input, and output (binary) data",
    "  file must exist and have correct length (may be copy of input). If",
    "  output header file does not exist, output header and binary files",
    "  are created.",
    "",
    "Command-line arguments:",
    "  in       = (string) filename of input RSF header file",
    "  out      = (string) filename of output RSF header file",
    "  power    = (real) power to which to raise Helmholtz operator,",
    "             default = 0.0 (identity)",
    "  datum    = (real) Dirichlet condition at x_1=x_{1,min}+datum",
    "             output zero for x_1 < x_{1,min}+datum, default = 0.0",
    "  scale1   = (real) scale applied to spectral second deriv in x_1",
    "             default = 1.0,",
    "  scale2   = (real) scale applied to spectral second deriv in x_2",
    "             default = 1.0,",
    NULL};

int main(int argc, char ** argv) {
    
    /* COMMAND-LINE PARAMS */
    char * in   = NULL;            /* input RSF header file */
    char * out  = NULL;            /* output RSF header file */
    //    char * outd = NULL;            /* output RSF data file (for write, if needed) */
    /* subtracted at start, added back at end */
    grid g;
    grid gout;
    
    IPNT n_arr;
    IPNT gs_arr;
    IPNT ord_arr;
    RPNT d_arr;
    RPNT o_arr;
    
    char * data_format;
    float scale1=1.0f, scale2=1.0f;
    
    float power;                   /* power to apply */
    float wt;
    float weightx, weightz, wx, wz;
    int nz, nx, nz2, nx2, ix, iz,nzc;
    
    //    FILE * fp;
    
    PARARRAY *par=ps_new();
    PARARRAY *gpar=ps_new();
    
    float * indata = NULL;
    float * outdata = NULL;
    fftwf_complex * f = NULL;
    fftwf_complex * gg = NULL;
    fftwf_complex *c=NULL;
    fftwf_complex *dd=NULL;
    fftwf_plan cfg=NULL, icfg=NULL;
    //    int i;
    
    xargc = argc;
    xargv = argv;
    requestdoc(1);
    
    if (ps_createargs(par, argc-1, argv+1)) {
        fprintf(stderr,"Error: HELM - failed to parse input \n");
        exit(1);
    }
    
    if (ps_flcstring(*par,"in",&in)) {
        printf("Error: HELM failed to parse in = [input RSF file]. ABORT.\n");
        exit(1);
    }
    
    if (ps_flcstring(*par,"out",&out)) {
        printf("Error: HELM failed to parse out = [output RSF file]. ABORT.\n");
        exit(1);
    }
    
    if (ps_createfile(gpar,in)) {
        fprintf(stderr,"Error: HELM failed to parse file %s\n",in);
        exit(1);
    }
    
    if (par_grid(&g, *gpar, stderr)) {
        fprintf(stderr,"Error: read from par_grid, file = %s\n",in);
        exit(1);
    }
    
    /* check dimension */
    if (g.dim != 2) {
        fprintf(stderr,"Error: input grid is not 2D\n");
        fprint_grid(stderr,g);
        exit(1);
    }

    /* initialize grid quantities */
    get_n(n_arr,g);
    get_gs(gs_arr,g);
    get_d(d_arr,g);
    get_o(o_arr,g);
    get_ord(ord_arr,g);
    data_format=NULL;
    ps_flcstring(*gpar,"data_format",&data_format);
    if (!data_format) {
        data_format=(char *)malloc(100*sizeof(char));
        strcpy(data_format,"native_float");
    }
    
    ps_flfloat(*par,"scale1",&scale1);
    ps_flfloat(*par,"scale2",&scale2);
    ps_flfloat(*par,"power",&power);

    
    /* initialize f2c ints */
    nz = n_arr[0];
    nx = n_arr[1];
    /* allocate data arrays */
    if (!(indata = (float *)malloc(get_datasize_grid(g)*sizeof(float)))) {
        fprintf(stderr,"Error: HELM - failed to allocate %d floats for input data\n",get_datasize_grid(g));
        exit(1);
    }
    if (!(outdata = (float *)malloc(get_datasize_grid(g)*sizeof(float)))) {
        fprintf(stderr,"Error: HELM - failed to allocate %d floats for output data\n",get_datasize_grid(g));
        exit(1);
    }
    
    if (read_grid(&gout, out, stderr)) {
        fprintf(stderr,"HELM: output header file = %s does not exist, or is not RSF header file\n",out);
        fprintf(stderr," - write new RSF header, data files\n");
        exit(1);
    }
    else {
        if (compare_grid(g,gout)) {
            fprintf(stderr,"Error: HELM - input, output grids do not match\n");
            exit(1);
        }
    }
    /* read input */
    if (rsfread(indata,gs_arr,n_arr,in,1.0,stderr,0)) {
        fprintf(stderr,"Error: HELM - failed to read %d floats from %s\n",get_datasize_grid(g),in);
        exit(1);
    }

    //ps_flint(*par,"pad1",&pad1); /* padding factor on the first axis */
    nzc = (nz+4)%2==1? (nz+5):(nz+4); //kiss_fft_next_fast_size((nz+1)*pad1/2)+2;
    fprintf(stderr,"nzc= %d\n", nzc);

    nz2=nzc; //2*(nzc-1);
    nzc=nz2;
    nx2 = (nx+4)%2==1? (nx+5):(nx+4); //kiss_fft_next_fast_size(nx);

    fprintf(stderr,"nz= %d\n", nz);
    fprintf(stderr,"nx= %d\n", nx);

    gg = (fftwf_complex*) malloc(sizeof(fftwf_complex) * nz2 * nx2);
    f = (fftwf_complex*) malloc(sizeof(fftwf_complex) * nz2 * nx2);
    c = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nz2 * nx2);
    dd = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nz2 * nx2);
    
    wt =  1.0/(nz2*nx2);
    for (ix=0; ix < nz2*nx2; ix++){
        f[ix][0]=0.0f;
        f[ix][1]=0.0f;
    }
    for (ix=0; ix < nx; ix++)
        for (iz=0; iz < nz; iz++)
            f[iz+(ix+2)*nz2+2][0]=indata[iz+ix*nz];
    
    for (ix=0; ix < nx2*nz2; ix++){
        gg[ix][0] = 0.0f;
        gg[ix][1] = 0.0f;
    }

    fprintf(stderr,"scale1=%f\n",scale1);
    fprintf(stderr,"scale2=%f\n",scale2);

    fprintf(stderr,"nz2= %d\n", nz2);
    fprintf(stderr,"nx2= %d\n", nx2);

// forward Fourier transform
    if (cfg==NULL) {
        cfg = fftwf_plan_dft_2d(nx2, nz2, f,  c, FFTW_FORWARD, FFTW_MEASURE);

        if (cfg==NULL) fprintf(stderr,"FFTW failure.\n");
    }
    
    fftwf_execute(cfg);

    fprintf(stderr,"c[100][0]= %f\n", c[100][0]);
    fprintf(stderr,"c[100][1]= %f\n", c[100][1]);


// inverse Fourier transform

    if (icfg==NULL){
        icfg = fftwf_plan_dft_2d(nx2, nz2, dd, gg, FFTW_BACKWARD, FFTW_MEASURE);

        if (icfg==NULL) fprintf(stderr,"FFTW failure.\n");
    }
    

    
    weightx = 2*3.14159265/nx2*scale2;
    weightx = weightx * weightx;
    
    weightz = 2*3.14159265/nzc*scale1;
    weightz = weightz * weightz;

    cerr << " scale1 = " << scale1 << endl;
    cerr << " scale2 = " << scale2 << endl;
    cerr << " wtz = " << weightz << endl;
    cerr << " wtx = " << weightx << endl;
    
    for (ix=0; ix<nx2/2; ix++) {
        wx = ix * ix;
        for (iz=0; iz<nzc/2; iz++) {
            wz = iz * iz;
            dd[ix*nzc+iz][0]=c[ix*nzc+iz][0]*(pow(1+weightx * wx + weightz * wz, power));
            dd[ix*nzc+iz][1]=c[ix*nzc+iz][1]*(pow(1+weightx * wx + weightz * wz, power));
        }
        for (iz=nzc/2; iz<nzc; iz++) {
            wz = (iz-nzc) * (iz-nzc);
            dd[ix*nzc+iz][0]=c[ix*nzc+iz][0]*(pow(1+weightx * wx + weightz * wz, power));
            dd[ix*nzc+iz][1]=c[ix*nzc+iz][1]*(pow(1+weightx * wx + weightz * wz, power));
        }

    }
    for (ix=nx2/2; ix<nx2; ix++) {
        wx = (ix-nx2) * (ix-nx2);
        for (iz=0; iz<nzc/2; iz++) {
            wz = iz * iz;
            dd[ix*nzc+iz][0]=c[ix*nzc+iz][0]*(pow(1+weightx * wx + weightz * wz, power));
            dd[ix*nzc+iz][1]=c[ix*nzc+iz][1]*(pow(1+weightx * wx + weightz * wz, power));
        }
        for (iz=nzc/2; iz<nzc; iz++) {
            wz = (iz-nzc) * (iz-nzc);
            dd[ix*nzc+iz][0]=c[ix*nzc+iz][0]*(pow(1+weightx * wx + weightz * wz, power));
            dd[ix*nzc+iz][1]=c[ix*nzc+iz][1]*(pow(1+weightx * wx + weightz * wz, power));
        }
    }
    
    
    fftwf_execute(icfg);

    fprintf(stderr,"gg[100][0]= %f\n", gg[100][0]);
    fprintf(stderr,"gg[100][1]= %f\n", gg[100][1]);
    
    /* copy to output array */
    for (ix=0; ix < nx; ix++) {
        for (iz=0; iz < nz; iz++) {
            outdata[iz+ix*nz]=wt*gg[iz+(ix+2)*nz2+2][0];
        }
    }
    
    /* write out result */
    if (rsfwrite(outdata,gs_arr,n_arr,out,1.0,stderr,0)) {
        fprintf(stderr,"Error: HELM failed to write data on %s\n",out);
        exit(1);
    }
    
    fftwf_destroy_plan(icfg);
    fftwf_destroy_plan(cfg);
    
    free(gg);
    free(f);
    fftwf_free(c);
    fftwf_free(dd);

    fprintf(stderr,"HELM: normal exit\n");
    exit(0);
}


/* Helmholtz */

#include <par.h>
#include <gridio.h>
#include <stdio.h>
#include <parser.h>
#include <f2c.h>
#include <iostream>
using namespace std;

/* helmholtz power function */
extern "C" void helm_(int,         /* bc */
		      integer *,   /* n1 */
		      integer *,   /* n2 */
		      float *,     /* d1 */
		      float *,     /* d2 */
		      float *,     /* w1 */
		      float *,     /* w2 */
		      float *,     /* p */
		      float *,     /* datum */
		      float *,     /* data in */
		      float *,     /* data out */
		      float *,     /* workspace */
		      integer *,   /* length of workspace */
		      integer *    /* error flag */
    );

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
    char * ref  = NULL;            /* reference RSF data file (optional) - if supplied, array is */
    /* subtracted at start, added back at end */
    grid g;
    grid gout;
    grid gref;

    IPNT n_arr;
    IPNT gs_arr;
    IPNT ord_arr;
    RPNT d_arr;
    RPNT o_arr;

    char * data_format;
    int scale;

    RPNT w_arr;                    /* weight array */
    float datum;                   /* datum depth */
    float power;                   /* power to apply */

    integer f2c_n1;
    integer f2c_n2;
    float * work = NULL;
    integer lenwork;
    integer ier;

//    FILE * fp;
  
    PARARRAY *par=ps_new();
    PARARRAY *gpar=ps_new();

    float * indata = NULL;
    float * outdata = NULL;
    float * refdata = NULL;

//    int i;
    int bc;

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

    ps_flcstring(*par,"ref",&ref);

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

    /* read remaining arguments */
    w_arr[0]=REAL_ONE;
    w_arr[1]=REAL_ONE;
    power = REAL_ZERO;
    datum = REAL_ZERO;
    bc = 0; // 0 = Neumann on sides, !0 = Dirichlet
    ps_flfloat(*par,"scale1",&(w_arr[0]));
    ps_flfloat(*par,"scale2",&(w_arr[1]));
    ps_flfloat(*par,"power",&power);
    ps_flfloat(*par,"datum",&datum);
    ps_flint(*par,"DirichletSides",&bc);

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
    ps_flint(*gpar,"scale",&scale);

    /* allocate data arrays */
    if (!(indata = (float *)malloc(get_datasize_grid(g)*sizeof(float)))) {
	fprintf(stderr,"Error: HELM - failed to allocate %d floats for input data\n",get_datasize_grid(g));
	exit(1);
    }
    if (!(outdata = (float *)malloc(get_datasize_grid(g)*sizeof(float)))) {
	fprintf(stderr,"Error: HELM - failed to allocate %d floats for output data\n",get_datasize_grid(g));
	exit(1);
    }
    if (ref) {
	if (!(refdata = (float *)malloc(get_datasize_grid(g)*sizeof(float)))) {
	    fprintf(stderr,"Error: HELM - failed to allocate %d floats for reference data\n",get_datasize_grid(g));
	    exit(1);
	}      
    }

    if (read_grid(&gout, out, stderr)) {
	fprintf(stderr,"HELM: output header file = %s does not exist, or is not RSF header file\n",out);
	fprintf(stderr," - write new RSF header, data files\n");
	exit(1);
	/*
	outd = (char *)malloc(strlen(out)+10);
	strcpy(outd,out);
	strcat(outd,"@");
	fp = iwave_fopen(&out,"w",NULL,stderr);
	fprintf(fp,"n1=%d\nn2=%d\nd1=%g\nd2=%g\no1=%g\no2=%g\ndata_format=%s\nz_axis=%d\nx_axis=%d\nscale=%d\nin=%s\n",n_arr[0],n_arr[1],d_arr[0],d_arr[1],o_arr[0],o_arr[1],data_format,ord_arr[0]+1,ord_arr[1]+1,scale,outd);
	fflush(fp);
	iwave_fclose(fp);
	fp = iwave_fopen(&outd,"w",NULL,stderr);
	for (i=0;i<get_datasize_grid(g);i++) outdata[i]=0.0f;
	fwrite(outdata,sizeof(float),get_datasize_grid(g),fp);
	fflush(fp);
	iwave_fclose(fp);
	*/
    }
    else {
	if (compare_grid(g,gout)) {
	    fprintf(stderr,"Error: HELM - input, output grids do not match\n");
	    exit(1);
	}
    }

    if (ref) {
	if (read_grid(&gref, ref, stderr)) {
	    fprintf(stderr,"Error: HELM - failed to read reference grid from %s\n",ref);
	    exit(1);
	}
	if (compare_grid(g,gref)) {
	    fprintf(stderr,"Error: HELM - input, output grids do not match\n");
	    exit(1);
	}
    }

    /* read input */
    if (rsfread(indata,gs_arr,n_arr,in,1.0,stderr,0)) {
	fprintf(stderr,"Error: HELM - failed to read %d floats from %s\n",get_datasize_grid(g),in);
	exit(1);    
    }

    /* transform */

    /* initialize f2c ints */
    f2c_n1 = n_arr[0];
    f2c_n2 = n_arr[1];

    /* initialize workspace */
    lenwork = 6*n_arr[1]*n_arr[0]+3*iwave_max(n_arr[1],2*n_arr[0])+21;
    if (!(work = (float *)malloc(lenwork*sizeof(float)))) {
	fprintf(stderr,"Error: HELM - failed to allocate %ld floats for work buffer\n",lenwork);
	exit(1);
    }
    ier = 0;

    helm_(bc,&f2c_n1,&f2c_n2,
	  &(d_arr[0]),&(d_arr[1]),
	  &(w_arr[0]),&(w_arr[1]),
	  &power,&datum,
	  indata,
	  outdata,
	  work,
	  &lenwork,
	  &ier);
    fprintf(stderr, "\n indata [100] = %f\n", indata[100]);
    fprintf(stderr, "\n outdata [100] = %f\n", outdata[100]);
    if (ier != 0) {
	fprintf(stderr,"Error: HELM returned error from Helmholtz function = %ld\n",ier);
	exit(1);
    }

    /* write out result */
    if (rsfwrite(outdata,gs_arr,n_arr,out,1.0,stderr,0)) {
	fprintf(stderr,"Error: HELM failed to write data on %s\n",out);
	exit(1);
    }

    fprintf(stderr,"HELM: normal exit\n");
    exit(0);
}


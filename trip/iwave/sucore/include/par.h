/* Copyright (c) Colorado School of Mines, 2002.*/
/* All rights reserved.                       */

/* par.h - include file for getpar, selfdoc, and error handling functions */

#ifndef PAR_H
#define PAR_H


/* INCLUDES */

#include "cwp.h"

/* GLOBAL DECLARATIONS */
extern int xargc; extern char **xargv;

/* ssize_t does not exist on some systems */
#ifndef ssize_t
#define ssize_t int
#endif


/* TYPEDEFS */
typedef enum {cwp_false, cwp_true} cwp_Bool;
typedef char *cwp_String;

typedef enum {BADFILETYPE = -1,
	TTY, DISK, DIRECTORY, TAPE, PIPE, FIFO, SOCKET, SYMLINK} FileType;
	
/* define structures for Hale's modeling */
typedef struct ReflectorSegmentStruct {
	float x;	/* x coordinate of segment midpoint */
	float z;	/* z coordinate of segment midpoint */
	float s;	/* x component of unit-normal-vector */
	float c;	/* z component of unit-normal-vector */
} ReflectorSegment;
typedef struct ReflectorStruct {
	int ns;			/* number of reflector segments */
	float ds;		/* segment length */
	float a;		/* amplitude of reflector */
	ReflectorSegment *rs;	/* array[ns] of reflector segments */
} Reflector;
typedef struct WaveletStruct {
	int lw;			/* length of wavelet */
	int iw;			/* index of first wavelet sample */
	float *wv;		/* wavelet sample values */
} Wavelet;


/* DEFINES */

/* getpar macros */
#define MUSTGETPARINT(x,y) if(!getparint(x,y)) suerr("must specify %s=",x)
#define MUSTGETPARFLOAT(x,y) if(!getparfloat(x,y)) suerr("must specify %s=",x)
#define MUSTGETPARSTRING(x,y) if(!getparstring(x,y)) suerr("must specify %s=",x)

#define STDIN (0)
#define STDOUT (1)
#define STDERR (2)


/* FUNCTION PROTOTYPES */

#ifdef __cplusplus  /* if C++, specify external C linkage */
extern "C" {
#endif

/* getpar parameter parsing */
void initargs (int argc, char **argv);
int getparint (char *name, int *p);
int getparuint (char *name, unsigned int *p);
int getparshort (char *name, short *p);
int getparushort (char *name, unsigned short *p);
int getparlong (char *name, long *p);
int getparulong (char *name, unsigned long *p);
int getparfloat (char *name, float *p);
int getpardouble (char *name, double *p);
int getparstring (char *name, char **p);
int getparstringarray (char *name, char **p);
int getnparint (int n, char *name, int *p);
int getnparuint (int n, char *name, unsigned int *p);
int getnparshort (int n, char *name, short *p);
int getnparushort (int n, char *name, unsigned short *p);
int getnparlong (int n, char *name, long *p);
int getnparulong (int n, char *name, unsigned long *p);
int getnparfloat (int n, char *name, float *p);
int getnpardouble (int n, char *name, double *p);
int getnparstring (int n, char *name, char **p);
int getnparstringarray (int n, char *name, char **p);
int getnpar (int n, char *name, char *type, void *ptr);
int countparname (char *name);
int countparval (char *name);
int countnparval (int n, char *name);

/* For ProMAX */
void getPar(char *name, char *type, void *ptr);

/* errors and warnings */
void suerr(char *fmt, ...);
void syssuerr(char *fmt, ...);
void suwarn(char *fmt, ...);

/* self documentation */
void pagedoc (void);
void requestdoc (int i);

/* system calls with error trapping */
int ecreat(char *path, int perms);
int efork(void);
int eopen(char *path, int flags, int perms);
int eclose(int fd);
int eunlink(char *path);
long elseek(int fd, long offset, int origin);
int epipe(int fd[2]);

ssize_t eread(int fd, char *buf, size_t nbytes);
ssize_t ewrite(int fd, char *buf, size_t nbytes);

/* system subroutine calls with error trapping */
FILE *efopen(const char *file, const char *mode);
FILE *efreopen(const char *file, const char *mode, FILE *stream1);
FILE *efdopen(int fd, const char *mode);
FILE *epopen(char *command, char *type);
int efclose(FILE *stream);
int epclose(FILE *stream);
int efflush(FILE *stream);
int eremove(const char *file);
int erename(const char *oldfile, const char* newfile);
int efseeko(FILE *stream, off_t offset, int origin);
int efseek(FILE *stream, long offset, int origin);
long eftell(FILE *stream);
off_t eftello(FILE *stream);
void erewind(FILE *stream);
FILE *etmpfile(void);
char *etmpnam(char *namebuffer);
void *emalloc(size_t size);
void *erealloc(void *memptr, size_t size);
void *ecalloc(size_t count, size_t size);
size_t efread(void *bufptr, size_t size, size_t count, FILE *stream);
size_t efwrite(void *bufptr, size_t size, size_t count, FILE *stream);

/* some SUN users may need to comment out the next two items, */
/* if your system does not have "fgetpos()" and "fsetpos()" defined. */
/* You will also need to comment out the lines defining the functions */
/* efgetpos() and efsetpos() in CWPROOT/src/par/lib/subcalls.c */

/* Modified: 21 June 1995: */
/* so that setting -DSUN_A in Makefile.config should make this unnecessary */
/* CWP: John Stockwell */

#ifndef SUN_A
int efgetpos(FILE *stream, fpos_t *position);
int efsetpos(FILE *stream, const fpos_t *position);
#endif

/* allocation with error trapping */
void *ealloc1 (size_t n1, size_t size);
void *erealloc1 (void *v, size_t n1, size_t size);
void **ealloc2 (size_t n1, size_t n2, size_t size);
void ***ealloc3 (size_t n1, size_t n2, size_t n3, size_t size);
void ****ealloc4 (size_t n1, size_t n2, size_t n3, size_t n4, size_t size);
void ****ealloc4 (size_t n1, size_t n2, size_t n3, size_t n4, size_t size);
void *****ealloc5 (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t size);
void ******ealloc6 (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5,
                  size_t n6, size_t size);

int *ealloc1int(size_t n1);
int *erealloc1int(int *v, size_t n1);
int **ealloc2int(size_t n1, size_t n2);
int ***ealloc3int(size_t n1, size_t n2, size_t n3);
float *ealloc1float(size_t n1);
float *erealloc1float(float *v, size_t n1);
float **ealloc2float(size_t n1, size_t n2);
float ***ealloc3float(size_t n1, size_t n2, size_t n3);

int ****ealloc4int(size_t n1, size_t n2, size_t n3, size_t n4);
int *****ealloc5int(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5);
float ****ealloc4float(size_t n1, size_t n2, size_t n3, size_t n4);
float *****ealloc5float(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5);
float ******ealloc6float(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5,
       size_t n6);

unsigned short *****ealloc5ushort(size_t n1, size_t n2,
		 size_t n3, size_t n4, size_t n5);
unsigned char *****ealloc5uchar(size_t n1, size_t n2,
		size_t n3, size_t n4, size_t n5);
unsigned short ******ealloc6ushort(size_t n1, size_t n2,
		 size_t n3, size_t n4, size_t n5, size_t n6);

double *ealloc1double(size_t n1);
double *erealloc1double(double *v, size_t n1);
double **ealloc2double(size_t n1, size_t n2);
double ***ealloc3double(size_t n1, size_t n2, size_t n3);

/* string to numeric conversion with error checking */
short eatoh(char *s);
unsigned short eatou(char *s);
int eatoi(char *s);
unsigned int eatop(char *s);
long eatol(char *s);
unsigned long eatov(char *s);
float eatof(char *s);
double eatod(char *s);

/* file type checking */
FileType filestat(int fd);
char *printstat(int fd);

/* Hale's modeling code */
void decodeReflectors (int *nrPtr,
	float **aPtr, int **nxzPtr, float ***xPtr, float ***zPtr);
int decodeReflector (char *string,
	float *aPtr, int *nxzPtr, float **xPtr, float **zPtr);
void breakReflectors (int *nr, float **ar, 
	int **nu, float ***xu, float ***zu);
void makeref (float dsmax, int nr, float *ar, 
	int *nu, float **xu, float **zu, Reflector **r);
void raylv2 (float v00, float dvdx, float dvdz,
	float x0, float z0, float x, float z,
	float *c, float *s, float *t, float *q);
void addsinc (float time, float amp,
	int nt, float dt, float ft, float *trace);
void makericker (float fpeak, float dt, Wavelet **w);

/* upwind eikonal stuff */
void eikpex (int na, float da, float r, float dr, 
	float sc[], float uc[], float wc[], float tc[],
	float sn[], float un[], float wn[], float tn[]);
void ray_theoretic_sigma (int na, float da, float r, float dr, 
	float uc[], float wc[], float sc[],
	float un[], float wn[], float sn[]);
void ray_theoretic_beta (int na, float da, float r, float dr, 
	float uc[], float wc[], float bc[],
	float un[], float wn[], float bn[]);
void eiktam (float xs, float zs, 
	int nz, float dz, float fz, int nx, float dx, float fx, float **vel,
	float **time, float **angle, float **sig, float **beta);

/* smoothing routines */
void dlsq_smoothing (int nt, int nx, int ift, int ilt, int ifx, int ilx,
        float r1, float r2, float rw, float **traces);
void SG_smoothing_filter (int np, int nl, int nr, int ld, int m, float *filter);
void rwa_smoothing_filter (int flag, int nl, int nr, float *filter);
void gaussian2d_smoothing (int nx, int nt, int nsx, int nst, float **data);
void gaussian1d_smoothing (int ns, int nsr, float *data);
void smooth_histogram (int nintlh, float *pdf);

/* function minimization */
void bracket_minimum(float *ax, float *bx, float *cx, float *fa,
                float *fb,float *fc, float (*func)(float));
float golden_bracket(float ax, float bx, float cx,
                        float (*f)(float), float tol,float *xmin);
float brent_bracket(float ax, float bx, float cx,
                float (*f)(float), float tol, float *xmin);

void linmin(float p[],float xi[],int n,float *fret, float (*func)());
void powell_minimization(float p[], float **xi, int n,
                float ftol,int *iter,float *fret,float (*func)());

/***** lincoeff -- linearized reflection coefficients */
/* type definitions */

typedef struct ErrorFlagStructure
    {
      float iso[5];
      float upper[2];
      float lower[2];
      float global[4];
      float angle[4];
  } ErrorFlag;

/* prototypes for functions defined */

float Rp(float ang, float azim, float kappa, float *rpp, ErrorFlag *rp_1st, ErrorFlag *rp_2nd, 
         int count);

float Rs(float ang, float azim, float kappa, float *rps1, float *rps2, 
       float *sv, float *sh, float *cphi, float *sphi, int i_hsp,
       ErrorFlag *rsv_1st, ErrorFlag *rsv_2nd, ErrorFlag *rsh_1st, ErrorFlag *rsh_2nd, int count);

float Iso_exact(int type, float vp1, float vs1, float rho1, 
		  float vp2, float vs2, float rho2, float ang);  

int Phi_rot(float *rs1,float *rs2,int iso_plane,float pb_x,float pb_y,float pb_z,float gs1_x,float gs1_y,
       float gs1_z,float gs2_x,float gs2_y,float gs2_z,float *CPhi1,float *SPhi1,float *CPhi2,float
*SPhi2);

/* there are other functions involved. However, they are declared inside of the corresponding
   routines */
/*** end lincoeff */

#ifdef __cplusplus  /* if C++ (external C linkage is being specified) */
}
#endif

#endif /* PAR_H */

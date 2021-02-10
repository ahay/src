#include <string.h>

/* This file (created by config) defines either RSF64BIT or RSF32BIT */
#include "ptr_sz.h"

/* Handle 32/64 bit architectures; sf_file is really a pointer */
#ifdef  RSF64BIT
#define RSFFILE LONGLONG
#else
#define RSFFILE INT
#endif

/* Support 64 bit file offsets (off_t == long long) */
#define OFFSET64BIT

#ifdef  OFFSET64BIT
#define OFFSETT LONGLONG
#else
#define OFFSETT LONG
#endif

#include "cfortran.h"

#include <rsf.h>

/* Fortran functions for command-line arguments - NON-STANDARD
   use #ifdefs if necessary */

#ifdef GFORTRAN

extern int _gfortran_iargc();
extern void _gfortran_getarg_i4(int *i,char *arg,int len);

#else

PROTOCCALLSFFUN0(INT,IARGC,iargc)
PROTOCCALLSFSUB2(GETARG,getarg,INT,PSTRING)

#define IARGC() CCALLSFFUN0(IARGC,iargc)
#define GETARG(I,N) CCALLSFSUB2(GETARG,getarg,INT,PSTRING,I,N)

#endif

void sf_init_f(void)
{
    int i, len, argc1;
    char** argv1, *argvi, arg[256];

#ifdef GFORTRAN
    argc1 = 1+_gfortran_iargc();
#else
    argc1 = 1+IARGC(); 
#endif   
    argv1 = (char**) sf_alloc(argc1+1,sizeof(char*));

    for (i=0; i < argc1; i++) {

#ifdef GFORTRAN
        _gfortran_getarg_i4(&i,arg,256);
	len = strchr(arg,' ')-arg+1;
#else
	GETARG(i,arg);
	len = strlen(arg)+1;
#endif
	argvi = argv1[i] = sf_charalloc(len);
	memcpy(argvi,arg,len);
        argvi[len-1]='\0';

	/* fprintf(stderr,"got argv[%d]=\"%s\" len=%d\n",i,argvi,len); */
    }    
    argv1[argc1] = NULL;

    sf_init(argc1,argv1);
}

bool sf_getbool_f (const char* key,/*@out@*/ int* ipar)
{
    bool par;
    
    if (!sf_getbool(key,&par)) return false;
    *ipar = par;
    return true;
}

bool sf_getbools_f (const char* key,/*@out@*/ int* ipar, size_t n)
{
    int i;
    bool *par;

    par = sf_boolalloc(n);
    
    if (!sf_getbools(key,par,n)) {
	free (par);
	return false;
    } else {
	for (i=0; i < n; i++) {
	    ipar[i] = par[i];
	}
	free(par);
	return true;
    }
}

#ifdef f2cFortran /* f2c adds an extra underscore */
#ifndef INTEL_COMPILER
#undef fcallsc
#define fcallsc(UN,LN) append_fcallsc(_,_,UN,LN)
#endif
#endif 

FCALLSCSUB0(sf_init_f,SF_INIT,sf_init)
FCALLSCSUB0(sf_parclose,SF_PARCLOSE,sf_parclose)
FCALLSCFUN0(STRING,sf_getprog,SF_GETPROG,sf_getprog)
FCALLSCFUN2(LOGICAL,sf_getint,SF_GETINT,sf_getint,STRING,PINT)
FCALLSCFUN3(LOGICAL,sf_getints,SF_GETINTS,sf_getints,STRING,INTV,INT)
FCALLSCFUN2(LOGICAL,sf_getfloat,SF_GETFLOAT,sf_getfloat,STRING,PFLOAT)
FCALLSCFUN3(LOGICAL,sf_getfloats,SF_GETFLOATS,sf_getfloats,STRING,FLOATV,INT)
FCALLSCFUN1(STRING,sf_getstring,SF_GETSTRING,sf_getstring,STRING)

#define sf_getstrings_STRV_A2 NUM_ELEM_ARG(3)
FCALLSCFUN3(LOGICAL,sf_getstrings,SF_GETSTRINGS,sf_getstrings,STRING,PSTRINGV,INT)

FCALLSCFUN2(LOGICAL,sf_getbool_f,SF_GETBOOL,sf_getbool,STRING,PLOGICAL)
FCALLSCFUN3(LOGICAL,sf_getbools_f,SF_GETBOOLS,sf_getbools,STRING,LOGICALV,INT)

FCALLSCFUN1(RSFFILE,sf_input,SF_INPUT,sf_input,STRING)
FCALLSCFUN1(RSFFILE,sf_output,SF_INPUT,sf_output,STRING)
FCALLSCFUN1(INT,sf_gettype,SF_GETTYPE,sf_gettype,RSFFILE)
FCALLSCFUN1(INT,sf_getform,SF_GETFORM,sf_getform,RSFFILE)
FCALLSCSUB2(sf_settype,SF_SETTYPE,sf_settype,RSFFILE,INT)
FCALLSCSUB2(sf_setformat,SF_SETFORMAT,sf_setformat,RSFFILE,STRING)
FCALLSCSUB1(sf_fileclose,SF_FILECLOSE,sf_fileclose,RSFFILE)
FCALLSCSUB2(sf_fileflush,SF_FILEFLUSH,sf_fileflush,RSFFILE,RSFFILE)

FCALLSCFUN3(LOGICAL,sf_histint,SF_HISTINT,sf_histint,RSFFILE,STRING,PINT)
FCALLSCFUN4(LOGICAL,sf_histints,SF_HISTINTS,sf_histints,RSFFILE,STRING,INTV,INT)
FCALLSCFUN3(LOGICAL,sf_histfloat,SF_HISTFLOAT,sf_histfloat,RSFFILE,STRING,PFLOAT)
FCALLSCFUN4(LOGICAL,sf_histfloats,SF_HISTFLOATS,sf_histfloats,RSFFILE,STRING,FLOATV,INT)
FCALLSCFUN3(LOGICAL,sf_histbool,SF_HISTBOOL,sf_histbool,RSFFILE,STRING,PLOGICAL)
FCALLSCFUN4(LOGICAL,sf_histbools,SF_HISTBOOLS,sf_histbools,RSFFILE,STRING,LOGICALV,INT)
FCALLSCFUN2(STRING,sf_histstring,SF_HISTSTRING,sf_histstring,RSFFILE,STRING)

FCALLSCSUB3(sf_putint,SF_PUTINT,sf_putint,RSFFILE,STRING,INT)
FCALLSCSUB4(sf_putints,SF_PUTINTS,sf_putints,RSFFILE,STRING,INTV,INT)
FCALLSCSUB3(sf_putfloat,SF_PUTFLOAT,sf_putfloat,RSFFILE,STRING,FLOAT)
FCALLSCSUB3(sf_putstring,SF_PUTSTRING,sf_putstring,RSFFILE,STRING,STRING)
FCALLSCSUB2(sf_putline,SF_PUTLINE,sf_putline,RSFFILE,STRING)

FCALLSCFUN1(OFFSETT,sf_bytes,SF_BYTES,sf_bytes,RSFFILE)
FCALLSCFUN1(OFFSETT,sf_tell,SF_TELL,sf_tell,RSFFILE)
FCALLSCSUB3(sf_seek,SF_SEEK,sf_seek,RSFFILE,OFFSETT,INT)
FCALLSCSUB2(sf_unpipe,SF_UNPIPE,sf_unpipe,RSFFILE,OFFSETT)
FCALLSCSUB0(sf_close,SF_CLOSE,sf_close)

FCALLSCSUB3(sf_floatwrite,SF_WRITE,sf_floatwrite,PFLOAT,INT,RSFFILE)
FCALLSCSUB3(sf_floatread,SF_READ,sf_floatread,PFLOAT,INT,RSFFILE)
FCALLSCSUB3(sf_intwrite,SF_WRITE,sf_intwrite,PINT,INT,RSFFILE)
FCALLSCSUB3(sf_intread,SF_READ,sf_intread,PINT,INT,RSFFILE)
FCALLSCSUB3(sf_complexwrite,SF_WRITE,sf_complexwrite,PVOID,INT,RSFFILE)
FCALLSCSUB3(sf_complexread,SF_READ,sf_complexread,PVOID,INT,RSFFILE)

FCALLSCFUN2(INT,sf_filedims,SF_FILEDIMS,sf_filedims,RSFFILE,INTV)
FCALLSCFUN1(OFFSETT,sf_filesize,SF_FILESIZE,sf_filesize,RSFFILE)
FCALLSCFUN2(OFFSETT,sf_leftsize,SF_LEFTSIZE,sf_leftsize,RSFFILE,INT)

FCALLSCSUB1(sf_error,SF_ERROR,sf_error,STRING)
FCALLSCSUB1(sf_warning,SF_WARNING,sf_warning,STRING)

/* Cosine transform */

FCALLSCSUB1(sf_cosft_init,SF_COSFT_INIT,sf_cosft_init,INT)
FCALLSCSUB0(sf_cosft_close,SF_COSFT_CLOSE,sf_cosft_close)
FCALLSCSUB3(sf_cosft_frw,SF_COSFT_FRW,sf_cosft_frw,PFLOAT,INT,INT)
FCALLSCSUB3(sf_cosft_inv,SF_COSFT_INV,sf_cosft_inv,PFLOAT,INT,INT)

/* Half-order differentiation */
FCALLSCSUB3(sf_halfint_init,SF_HALFINT_INIT,sf_halfint_init,LOGICAL,INT,FLOAT)
FCALLSCSUB6(sf_halfint_lop,SF_HALFINT_LOP,sf_halfint_lop,LOGICAL,LOGICAL,INT,INT,PFLOAT,PFLOAT)

/* Stretch */
FCALLSCFUN5(RSFFILE,sf_stretch4_init,SF_STRETCH4_INIT,sf_stretch4_init,INT,FLOAT,FLOAT,INT,FLOAT)
FCALLSCSUB3(sf_stretch4_define,SF_STRETCH4_DEFINE,sf_stretch4_define,RSFFILE,PFLOAT,LOGICAL)
FCALLSCSUB4(sf_stretch4_apply,SF_STRETCH4_APPLY,sf_stretch4_apply,LOGICAL,RSFFILE,PFLOAT,PFLOAT)
FCALLSCSUB4(sf_stretch4_apply_adj,SF_STRETCH4_APPLY_ADJ,sf_stretch4_apply_adj,LOGICAL,RSFFILE,PFLOAT,PFLOAT)

/* Adjnull */
FCALLSCSUB6(sf_adjnull,SF_ADJNULL,sf_adjnull,LOGICAL,LOGICAL,INT,INT,PFLOAT,PFLOAT)

/* 	$Id: fortran.c 2069 2006-07-26 05:16:46Z sfomel $	 */

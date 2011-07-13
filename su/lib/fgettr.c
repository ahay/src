/* Copyright (c) Colorado School of Mines, 2010.*/
/* All rights reserved.                       */

/* FGETTR: $Revision: 1.43 $; $Date: 2006/03/29 22:44:17 $	*/

/*********************** self documentation **********************/
/****************************************************************************
FGETTR - Routines to get an SU trace from a file 

fgettr		get a fixed-length segy trace from a file by file pointer
fvgettr		get a variable-length segy trace from a file by file pointer
fgettra		get a fixed-length trace from disk file by trace number
gettr		macro using fgettr to get a trace from stdin
vgettr		macro using vfgettr to get a trace from stdin
gettra		macro using fgettra to get a trace from stdin by trace number
 
*****************************************************************************
Function Prototype:
int fgettr(FILE *fp, segy *tp);
int fvgettr(FILE *fp, segy *tp);
int fgettra(FILE *fp, segy *tp, int itr);

*****************************************************************************
Returns:
fgettr, fvgettr:
int: number of bytes read on current trace (0 after last trace)

fgettra:
int: number of traces in disk file
 
Macros defined in segy.h
#define gettr(x)	fgettr(stdin, (x))
#define vgettr(x)	fgettr(stdin, (x))

Usage example:
 	segy tr;
 	...
 	while (gettr(&tr)) {
 		tr.offset = abs(tr.offset);
 		puttr(&tr);
 	}
 	...

	*****************************************************************************
Authors: SEP: Einar Kjartansson, Stew Levin CWP: Shuki Ronen, Jack Cohen
****************************************************************************/
/*
 * Revised: 7/2/95 Stewart A. Levin   Mobil
 *     Major rewrite:  Use xdr library for portable su output file
 *     format.   Merge fgettr and fgettra into same source file.
 *     Make input from multiple streams work (at long last!).
 * Revised: 11/22/95 Stewart A. Levin  Mobil
 *     Always set ntr for DISK input.  This fixes susort failure.
 * Revised: 1/9/96  jkc CWP
 *     Set lastfp on nread <=0 return, too.
 * Revised: 28 Mar, 2006 Stewart A. Levin   Landmark Graphics
 *     Reworked XDR to support random seeks on > 2GB files
 *     and to read big endian SHORTPACK data on little endian machines.
 */
/**************** end self doc ********************************/

#include <stdio.h>
/*^*/

#include <rsf.h>

#include "segy.h"
/*^*/

#include "filestat.h"
#include "fgettr.h"

extern int fseeko(FILE *stream, off_t offset, int whence);
extern off_t ftello (FILE *stream);

#ifndef TEST

#ifdef SU_LINE_HEADER

int out_line_hdr=1;

#else

int out_line_hdr=0;

#endif

/* DEFINES */
#define gettr(x)	fgettr(stdin, (x))
#define vgettr(x)	fvgettr(stdin, (x))
#define puttr(x)	fputtr(stdout, (x))
#define vputtr(x)	fvputtr(stdout, (x))
#define gettra(x, y)    fgettra(stdin, (x), (y))

int isebcdic_txt( unsigned char* ,int len);
int isascii_txt( unsigned char* ,int len);

int in_line_hdr=0;

unsigned char su_text_hdr[3200];
bhed su_binary_hdr;

#define HDRBYTES       240     /* Bytes in the trace header */

#ifdef SUXDR
/********************************************
code using XDR
*********************************************/

static struct insegyinfo {
    FILE *infp;		     /* FILE * ptr for search	 */
    struct insegyinfo *nextinfo; /* linked list pointer	 */
    off_t itr;		     /* number of traces read	 */
    int nsfirst;		     /* samples from 1st header	 */
    unsigned short bytesper;     /* bytes per datum		 */
    size_t nsegy; 		     /* segy bytes from nsfirst	 */
    off_t ntr;		     /* traces in input,if known */
    FileType ftype;		     /* file type of input *fp	 */
    XDR *segy_xdr;		     /* allocated XDR structure  */
    char *buf;		     /* buffer for trace I/O     */
    size_t bufstart;	     /* "offset" of start of buf */
} *insegylist = (struct insegyinfo *) NULL;

static FILE *lastfp = (FILE *) NULL;
static struct insegyinfo *infoptr, **oldinfoptr;



static
void searchlist(FILE *fp)
{
    oldinfoptr = &insegylist;
    for(infoptr = insegylist; infoptr != ((struct insegyinfo *) NULL);
	infoptr = infoptr->nextinfo) {
	if(fp == infoptr->infp) break;
	oldinfoptr = &infoptr->nextinfo;
    }
}

static
int dataread(struct insegyinfo *iptr, segy *tp, cwp_Bool fixed_length)
{
    unsigned int nsread = fixed_length?iptr->nsfirst:tp->ns;
    unsigned int databytes = infoptr->bytesper*nsread;
    int nread;
    int itest = 1;
    char *ctest = (char *) (&itest);


    /* read trace data */
    switch(tp->trid) {
	case CHARPACK:
	    nread = efread((char *) (&((tp->data)[0])),1,databytes,
			   iptr->infp);
	case SHORTPACK:
	    nread = efread((char *) (&((tp->data)[0])),1,databytes,
			   iptr->infp);
	    if(ctest[0]) swab((char *) (&((tp->data)[0])),
			      (char *) (&((tp->data)[0])),
			      databytes);
	    break;
	default:
	    nread = efread(((char *) (iptr->buf))+HDRBYTES,1,databytes,
			   iptr->infp);
	    if(nread != databytes || FALSE == xdr_vector(iptr->segy_xdr,
							 (char *) (&((tp->data)[0])),
							 nsread,sizeof(float),(xdrproc_t) xdr_float))
		nread = 0;
	    else
		nread = databytes;

	    break;
    }
	
    if(nread > 0 && nread != databytes) 
	err("%s: on trace #%ld, tried to read %d bytes, "
	    "read %d bytes",
	    __FILE__, (infoptr->itr)+1, databytes, nread);
	
    return(nread);
}


static
int fgettr_internal(FILE *fp, segy *tp, cwp_Bool fixed_length)
{
    int nread;		/* bytes seen by fread calls	*/
    off_t curoff;

    if(fp != lastfp)  /* search linked list for possible alternative */
	searchlist(fp);

    if (infoptr == ((struct insegyinfo *) NULL)) {
	/* initialize new segy input stream */
	unsigned int databytes;	/* bytes from nsfirst		*/

	/* allocate new segy input information table */
	*oldinfoptr = (struct insegyinfo *)
	    malloc(sizeof(struct insegyinfo));
	infoptr = *oldinfoptr;
	infoptr->nextinfo = (struct insegyinfo *) NULL;
	/* save FILE * ptr */
	infoptr->infp = fp;
	infoptr->itr = 0;
	infoptr->ntr = -1;
	/* allocate XDR struct and associate FILE * ptr */
	infoptr->segy_xdr = (XDR *) malloc(sizeof(XDR));

	switch (infoptr->ftype = filestat(fileno(fp))) {
	    case DIRECTORY:
		err("%s: segy input can't be a directory", __FILE__);
	    case TTY:
		err("%s: segy input can't be tty", __FILE__);
		break;
	    default: /* the rest are ok */
		break;
	}
	/* xdrstdio_create(infoptr->segy_xdr,fp,XDR_DECODE); */
	infoptr->buf = ealloc1(sizeof(segy),sizeof(char));
	xdrmem_create(infoptr->segy_xdr, infoptr->buf, sizeof(segy), XDR_DECODE);
	infoptr->bufstart = xdr_getpos(infoptr->segy_xdr);

	/* retrieve segy trace header */
	nread = efread(infoptr->buf ,1 ,HDRBYTES ,infoptr->infp);

	if(nread != HDRBYTES || FALSE == xdrhdrsub(infoptr->segy_xdr,tp))
	    err("%s: bad first header", __FILE__);
		
	/* Have the header, now for the data */
	infoptr->nsfirst = tp->ns;
	if (infoptr->nsfirst > SU_NFLTS)
	    err("%s: unable to handle %d > %d samples per trace",
		__FILE__, infoptr->nsfirst, SU_NFLTS);

	switch(tp->trid) {
	    case CHARPACK:
		infoptr->bytesper = sizeof(char); break;
	    case SHORTPACK:
		infoptr->bytesper = 2*sizeof(char); break;
	    default:
		infoptr->bytesper = BYTES_PER_XDR_UNIT; break;
	}

	databytes = infoptr->bytesper * tp->ns;

	infoptr->nsegy = databytes + HDRBYTES;

	nread = dataread(infoptr, tp, fixed_length);


	switch (nread) {
	    case 0:   err("%s: no data on first trace", __FILE__);
	    default:  if (nread != databytes)
		err("%s: first trace: tried to read %d bytes "
		    "read %d bytes",
		    __FILE__, databytes, nread);
	    else nread += HDRBYTES;
	}

	if(infoptr->ftype == DISK) { /* compute ntr */
	    curoff = ftello(fp); 
	    fseeko(fp,(off_t) 0,SEEK_END);
	    infoptr->ntr = ftello(fp)/infoptr->nsegy;
	    fseeko(fp, curoff, SEEK_SET); /* restore location */
	}


    } else { /* Not first entry */

	xdr_setpos(infoptr->segy_xdr, infoptr->bufstart);
	nread = efread(infoptr->buf ,1 ,HDRBYTES ,infoptr->infp);
	if ( nread != HDRBYTES || FALSE == xdrhdrsub(infoptr->segy_xdr,tp)) nread=0;
	if(nread == HDRBYTES)
	    nread += dataread(infoptr, tp, fixed_length);
	if (nread <= 0) {
	    lastfp = infoptr->infp;
	    return 0;
	}



	if (fixed_length && (tp->ns != infoptr->nsfirst))
	    err("%s: on trace #%ld, "
		"number of samples in header (%d) "
		"differs from number for first trace (%d)", 
		__FILE__, infoptr->itr, tp->ns, infoptr->nsfirst);
    }

    ++(infoptr->itr);
    lastfp = infoptr->infp;
    return (nread);
}

int fgettr(FILE *fp, segy *tp)
/*< get a fixed-length segy trace from a file by file pointer >*/
{
    return(fgettr_internal(fp,tp,cwp_true));
}

int fvgettr(FILE *fp, segy *tp)
{
    return(fgettr_internal(fp,tp,cwp_false));
}
int fgettra(FILE *fp, segy *tp, int itr)
{
    int nread;
    if(lastfp != fp) /* search for match */
	searchlist(fp);

    if(infoptr == (struct insegyinfo *) NULL) {
	/* get first trace */
	if(0 >= fgettr(fp, tp)) return(0); /* error return */

	switch(infoptr->ftype) {
	    case TTY:
		warn("stdin not redirected");
		break;
	    case DISK: /* correct */
		break;
	    default:
		err("%s: input must be disk file",__FILE__);
		break;
	}

	fseeko(fp,(off_t) 0,SEEK_END);
	infoptr->ntr = ftello(fp)/infoptr->nsegy;
    } /* end first entry initialization */

    /* Check on requested trace number */
    if(itr >= infoptr->ntr) err("%s: trying to read off end of file",__FILE__);

    /* Position file pointer at start of requested trace */
    if(0 > fseeko(fp, ((off_t) itr) * ((off_t) infoptr->nsegy), SEEK_SET)) {
	err("%s: unable to seek xdr disk file to trace %d",__FILE__,itr);
    }

    nread=fgettr(fp, tp);
    if(nread != infoptr->nsegy)
	err("%s: read %d bytes with %d bytes in trace",
	    __FILE__,nread,infoptr->nsegy);

    if(tp->ns != infoptr->nsfirst)
	warn("%s: header ns field = %d differs from first trace = %d",
	     __FILE__,tp->ns,infoptr->nsfirst);

    return(infoptr->ntr);
}


#else 
/**********************************************************
code without  XDR 
***********************************************************/

static struct insegyinfo {
    FILE *infp;		  /* FILE * ptr for search	 */
    struct insegyinfo *nextinfo; /* linked list pointer	*/
    off_t itr;	     /* number of traces read	 */
    int nsfirst;		     /* samples from 1st header	 */
    unsigned short bytesper;     /* bytes per datum		 */
    size_t  nsegy; 		     /* segy bytes from nsfirst	 */
    off_t ntr;		     /* traces in input,if known */
    FileType ftype;		     /* file type of input *fp	 */
} *insegylist = (struct insegyinfo *) NULL;

static FILE *lastfp = (FILE *) NULL;
static struct insegyinfo *infoptr, **oldinfoptr;

static
void searchlist(FILE *fp)
{
    oldinfoptr = &insegylist;
    for(infoptr = insegylist; infoptr != ((struct insegyinfo *) NULL);
	infoptr = infoptr->nextinfo) {
	if(fp == infoptr->infp) break;
	oldinfoptr = &infoptr->nextinfo;
    }
}

static
int dataread(segy *tp, struct insegyinfo *iptr, cwp_Bool fixed_length)
{
    unsigned int nsread = fixed_length?iptr->nsfirst:tp->ns;
    unsigned int databytes = infoptr->bytesper*nsread;
    int nread = (int) fread((char *) (&((tp->data)[0])),1, databytes,
			    iptr->infp);

    if(nread > 0 && nread != databytes) 
	sf_error("%s: on trace #%ld, tried to read %d bytes, "
		 "read %d bytes ",
		 __FILE__, (infoptr->itr)+1, databytes, nread);
	
    return(nread);
}


static
int fgettr_internal(FILE *fp, segy *tp, cwp_Bool fixed_length) {
    int nread;  /* bytes seen by fread calls  */
    unsigned char buf[240];  /* buffer for test for line header */

    /* search linked list for possible alternative */
    if(fp != lastfp)  searchlist(fp);

    if (infoptr == ((struct insegyinfo *) NULL)) {
	/* initialize new segy input stream */
	unsigned int databytes; /* bytes from nsfirst   */

	/* allocate new segy input information table */
	*oldinfoptr = (struct insegyinfo *)
	    malloc(sizeof(struct insegyinfo));
	infoptr = *oldinfoptr;
	infoptr->nextinfo = (struct insegyinfo *) NULL;
	infoptr->infp = fp;  /* save FILE * ptr */
	infoptr->itr = 0;
	infoptr->ntr = -1;
	
	switch (infoptr->ftype = filestat(fileno(fp))) {
	    case DIRECTORY:
		sf_error("%s: segy input can't be a directory", __FILE__);

	    case TTY:
		sf_error("%s: segy input can't be tty", __FILE__);

	    default:
		/* all others are ok */
		break;
	}

/*--------------------------------------------------------------------*\
  Check for the presence of a line header and set a flag if one is
  found. The decision of what to do will be delayed until the call
  to fputtr(). This allows us to accept data w/ or w/o a line
  header.

  Reginald H. Beardsley rhb@acm.org
  \*--------------------------------------------------------------------*/

	/* Attempt to get a text header */

	nread=fread(buf ,1 ,HDRBYTES ,infoptr->infp);

	switch( nread ){

	    case 0:   
		return 0; /* no traces; trap in mains */

	    default:  

		if (nread < HDRBYTES ){
		    return 0; 

		}else if( isascii_txt( buf ,HDRBYTES  )
			  || isebcdic_txt( buf ,HDRBYTES  ) ){
		    in_line_hdr = 1;
		    memcpy( su_text_hdr ,buf ,HDRBYTES );
		    nread += fread(&(su_text_hdr[HDRBYTES]) ,1 
				   ,3200-HDRBYTES ,infoptr->infp);

		}else{
		    in_line_hdr=0;
		    memcpy( tp ,buf ,HDRBYTES );

		}
	}		

	
	if( in_line_hdr ){

	    /* Get the binary header */
	    nread = fread(&su_binary_hdr, 1, sizeof(bhed), infoptr->infp);
	    switch( nread ){
		case 0:   
		    return 0; /* no traces; trap in mains */

		default:  
		    if (nread != sizeof(su_binary_hdr)){
			sf_error("%s:%d bad binary header" , __FILE__ ,__LINE__ );
		    }
	    }	 

	    /* Get the first trace header */
	    nread = fread(tp, 1, HDRBYTES, infoptr->infp);
	    switch( nread ){
		case 0:   
		    return 0; /* no traces; trap in mains */

		default:  
		    if (nread != HDRBYTES){ 
			sf_error("%s: bad first header", __FILE__);
		    }
	    }	 


	}

	/* Have the header, now for the data */
	infoptr->nsfirst = tp->ns;
	if (infoptr->nsfirst > SU_NFLTS){
	    sf_error("%s: unable to handle %d > %d samples per trace",
		__FILE__, infoptr->nsfirst, SU_NFLTS);
	}

	switch (tp->trid) {
	    case CHARPACK:
		infoptr->bytesper = sizeof(char); break;
	    case SHORTPACK:
		infoptr->bytesper = 2*sizeof(char); break;
	    default:
		infoptr->bytesper = sizeof(float); break;
	}

	databytes = infoptr->bytesper * tp->ns;

	infoptr->nsegy = HDRBYTES + databytes;


	/* Inconvenient to bump nread here; do it in the switch */
	nread = dataread(tp, infoptr, fixed_length);
	
	switch (nread) {
	    case 0:   
		sf_error("%s: no data on first trace", __FILE__);

	    default:  
		if (nread != databytes){
		    sf_error("%s: first trace: " "read only %d bytes of %u",
			     __FILE__, nread, databytes);

		}else{
		    nread += HDRBYTES;
		}
	}
	
	
	if (infoptr->ftype == DISK) { /* compute ntr */
	    fseeko(fp, (off_t) 0LL,SEEK_END);

	    if( in_line_hdr ){
		infoptr->ntr = (ftello(fp)-3600)/infoptr->nsegy;
		fseeko(fp, (off_t) 3600+infoptr->nsegy,SEEK_SET);

	    }else{
		infoptr->ntr = ftello(fp)/infoptr->nsegy;
		fseeko(fp, (off_t) infoptr->nsegy ,SEEK_SET);

	    }
	 
	}

    }else{

	/* not first trace */

/*--------------------------------------------------------------------*\
  A number of programs seek on the input file using either fseek(3c)
  or rewind(3c) and then expect to read trace data.  As a consequence
  we need to check and offset the filepointer if needed.
  \*--------------------------------------------------------------------*/

	if( in_line_hdr && ftello( infoptr->infp ) == 0 ){

	    fseeko( infoptr->infp ,(off_t)3600L ,SEEK_SET );

	}

	nread = (int) fread(tp, 1, HDRBYTES, infoptr->infp);

	switch( nread ){

	    case 0: 
		lastfp = infoptr->infp;
		return 0; /* finished */

	    default:  
		if (nread != HDRBYTES){
		    sf_error("%s: on trace #%ld read %d bytes expected %d bytes",
			     __FILE__,(infoptr->itr)+1,nread,HDRBYTES);
		}
	}

	nread += dataread(tp, infoptr, fixed_length);

	if (fixed_length && (tp->ns != infoptr->nsfirst)){
	    sf_error("%s: on trace #%ld number of samples in header (%d) differs from number for first trace (%d)"
		     ,__FILE__, (infoptr->itr)+1, tp->ns,

		     infoptr->nsfirst);
	}
    }

    ++(infoptr->itr);
    lastfp = infoptr->infp;
    return (nread);
}


int fgettr(FILE *fp, segy *tp)
{
    return(fgettr_internal(fp,tp,cwp_true));
}

int fvgettr(FILE *fp, segy *tp)
{
    return(fgettr_internal(fp,tp,cwp_false));
}

int fgettra(FILE *fp, segy *tp, int itr)
{
    int nread;
    if(lastfp != fp)  searchlist(fp);  /* search for match */
		
	
    if(infoptr == (struct insegyinfo *) NULL) {
	/* get first trace */
	if(0 >= fgettr(fp, tp)) return(0); /* error return */

	switch(infoptr->ftype) {
	    case TTY:
		sf_warning("stdin not redirected");
		break;
	    case DISK:	/* correct */
		break;
	    default:
		sf_error("%s: input must be disk file",__FILE__);
	}
	


	fseeko(fp,(off_t) 0LL,SEEK_END);
	if( in_line_hdr ){
	    infoptr->ntr = (off_t)((ftello(fp)-3600)/infoptr->nsegy);
	    fseeko(fp, (off_t) 3600+infoptr->nsegy,SEEK_SET);
	}else{
	    infoptr->ntr = (off_t)(ftello(fp)/infoptr->nsegy);
	    fseeko(fp, (off_t) infoptr->nsegy,SEEK_SET);
	}
    } /* end first entry initialization */
	
    /* Check on requested trace number */
    if(itr >= infoptr->ntr)
	sf_error("%s: trying to read off end of file",__FILE__);
	
    /* Position file pointer at start of requested trace */
    if( in_line_hdr ){
	fseeko(fp, (off_t) 3600+itr*infoptr->nsegy,SEEK_SET);
    }else{
	fseeko(fp, (off_t) itr*infoptr->nsegy,SEEK_SET);
    }
	
    nread=fgettr(fp, tp); /* let fgettr do the work */
    if(nread != infoptr->nsegy)
	sf_error("%s: read %d bytes in trace of %d bytes",
		 __FILE__,nread,infoptr->nsegy);
	
    if(tp->ns != infoptr->nsfirst)
	sf_warning("%s: header ns field = %d differs from first trace = %d",
		   __FILE__,tp->ns,infoptr->nsfirst);
	
    return(infoptr->ntr);
}

#endif /* end of XDR choice */
 
/*====================================================================*\
  These functions determine the presence of a SEGY text header based
  on the character set statistics of the first 3200 bytes in the
  file.

  Reginald H. Beardsley			    rhb@acm.org

  \*====================================================================*/

static  unsigned char asciitext[256] = {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 
};

  
static  unsigned char ebcdictext[256] = {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,
    1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,
    1,1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,
    0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,
    0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,
    0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,
    0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,
    0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,
    0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,
    0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,
    1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0 
};

  
int isascii_txt( unsigned char* a ,int n){

    int ascii  = 0;
    int ebcdic = 0;
    int i;

    for (i=0; i<n; i++){

	ascii += asciitext[a[i]];
	ebcdic += ebcdictext[a[i]];
    }

    if( ascii > ebcdic && ascii > (int)(n*0.9) ){
	return(1);
    }else{
	return(0);
    }
}

  
int isebcdic_txt( unsigned char* a ,int n ){

    int ascii  = 0;
    int ebcdic = 0;
    int i;

    for (i=0; i<n; i++){

	ascii += asciitext[a[i]];
	ebcdic += ebcdictext[a[i]];
    }

    if( ascii < ebcdic && ebcdic > (int)(n*0.9) ){
	return(1);
    }else{
	return(0);
    }
}

#endif


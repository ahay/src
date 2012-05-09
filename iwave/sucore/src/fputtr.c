/* Copyright (c) Colorado School of Mines, 2006.*/
/* All rights reserved.                       */

/* FPUTTR: $Revision: 1.31 $; $Date: 2005/02/07 21:16:03 $	*/


/*********************** self documentation **********************/
/****************************************************************************
FPUTTR - Routines to put an SU trace to a file 

fputtr		put a segy trace to a file by file pointer
fvputtr		put a segy trace to a file by file pointer (variable ns)
puttr		macro using fputtr to put a trace to stdin
vputtr		macro using fputtr to put a trace to stdin (variable ns)
 
*****************************************************************************
Function Prototype:
void fputtr(FILE *fp, segy *tp);
void fvputtr(FILE *fp, segy *tp);

*****************************************************************************
Returns:

	void
 
*****************************************************************************
Notes:

The functions puttr(x) vputtr(x) are macros defined in segy.h
#define puttr(x)	fputtr(stdin, (x))
#define vputtr(x)	fputtr(stdin, (x))

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
/**************** end self doc ********************************/


#include "su.h"
#include "segy.h"
#include "header.h"

extern int in_line_hdr;
extern int out_line_hdr;
extern char  su_text_hdr[3200];
extern bhed  su_binary_hdr;


#ifdef SUXDR

/*
 * Revised:  7/3/95  Stewart A. Levin  (Mobil)
 *          Major rewrite:  Use xdr library for portable su file format.
 *          Make multiple output streams work (at long last!).
 */
#include "su_xdr.h"
#include "header.h"

static struct outsegyinfo {
	FILE *outfp;		      /* FILE * ptr for search		*/
	struct outsegyinfo *nextinfo; /* linked list pointer     	*/
	unsigned off_t itr;	      /* number of traces written	*/
	unsigned int nsfirst;         /* nsamp from first trace		*/
	unsigned short bytesper;   /* bytes per datum		*/
	FileType ftype;		      /* file type of output *fp	*/
	XDR *segy_xdr;		      /* allocated XDR structure 	*/
} *outsegylist = (struct outsegyinfo *) NULL;

static FILE *lastfp = (FILE *) NULL;
static struct outsegyinfo *infoptr, **oldinfoptr;

static
void searchlist(FILE *fp)
{
	oldinfoptr = &outsegylist;
	for (infoptr = outsegylist; infoptr != ((struct outsegyinfo *) NULL);
	    infoptr = infoptr->nextinfo) {
		if (fp == infoptr->outfp) break;
		oldinfoptr = &infoptr->nextinfo;
	}
}

static
int datawrite(struct outsegyinfo *iptr, segy *tp, cwp_Bool fixed_length)
{
	int nwritten;
	unsigned int nstobewritten = fixed_length?iptr->nsfirst:tp->ns;
	unsigned int databytes = iptr->bytesper*nstobewritten;
	
	/* write trace data */
	switch(tp->trid) {
	case CHARPACK:
	case SHORTPACK:
		nwritten = efwrite((char *) (&((tp->data)[0])),1,databytes,
				  iptr->outfp);
	break;
	default:
		if(FALSE == xdr_vector(iptr->segy_xdr,
				       (char *) (&((tp->data)[0])),
				       nstobewritten,sizeof(float),(xdrproc_t) xdr_float))
			nwritten = 0;
		else
			nwritten = databytes;
	}

	return(nwritten);
}

void fputtr_internal(FILE *fp, segy *tp, cwp_Bool fixed_length)
{
        unsigned int databytes;         /* bytes from nsfirst           */
        int nwritten;                   /* bytes seen by fwrite calls   */
	
	/* search linked list for possible alternative */
	if(fp != lastfp)  searchlist(fp);
	
	if (infoptr == ((struct outsegyinfo *) NULL)) {
		/* initialize new segy output stream */
		
		/* allocate new segy output information table */
		*oldinfoptr = (struct outsegyinfo *)
			malloc(sizeof(struct outsegyinfo));
		infoptr = *oldinfoptr;
		infoptr->nextinfo = (struct outsegyinfo *) NULL;
		/* save FILE * ptr */
		infoptr->outfp = fp;
		infoptr->itr = 0;
		/* allocate XDR struct and associate FILE * ptr */
		infoptr->segy_xdr = (XDR *) malloc(sizeof(XDR));
		
		switch (infoptr->ftype = filestat(fileno(fp))) {
		case DIRECTORY:
			suerr("%s: segy output can't be a directory", __FILE__);
		case TTY:
			suerr("%s: segy output can't be tty", __FILE__);
		break;
		default:  /* the rest are ok */
		break;
		}
		xdrstdio_create(infoptr->segy_xdr,fp,XDR_ENCODE);
		
		/* Sanity check the segy header */
		infoptr->nsfirst = tp->ns;
		if (infoptr->nsfirst > SU_NFLTS)
			suerr("%s: unable to handle %d > %d samples per trace",
			    __FILE__, infoptr->nsfirst, SU_NFLTS);

		switch(tp->trid) {
		case CHARPACK:
		   infoptr->bytesper = sizeof(char); break;
		case SHORTPACK:
		   infoptr->bytesper = 2*sizeof(char); break;
		default:
		   infoptr->bytesper = BYTES_PER_XDR_UNIT; break;
		}

	}

	databytes = infoptr->bytesper * (fixed_length?infoptr->nsfirst:tp->ns);
	if(FALSE == xdrhdrsub(infoptr->segy_xdr,tp)) 
		suerr("%s: unable to write header on trace #%ld",
		    __FILE__, (infoptr->itr)+1);
	
        nwritten = datawrite(infoptr, tp, fixed_length);

        if (nwritten != databytes)
                suerr("%s: on trace #%ld, tried to write %d bytes, "
                    "wrote %d bytes",
                    __FILE__, (infoptr->itr)+1, databytes, nwritten);
	
	++infoptr->itr;
	lastfp = infoptr->outfp;
}
void fputtr(FILE *fp, segy *tp)
{
 fputtr_internal(fp,tp,cwp_true);
}

void fvputtr(FILE *fp, segy *tp)
{
 fputtr_internal(fp,tp,cwp_false);
}


#else
/**********************************************************************
code without XDR
**********************************************************************/

#include "su.h"
#include "segy.h"
#include "header.h"

/* part of Beardsley's hack - omitted 15.05.09 WWS */
/*
static char hdr_str[88];
static int i=0;
*/

static struct outsegyinfo {
	FILE *outfp;		      /* FILE * ptr for search		*/
	struct outsegyinfo *nextinfo; /* linked list pointer    	*/
	unsigned off_t itr;	      /* number of traces written	*/
	unsigned int nsfirst;         /* nsamp from first trace		*/
	unsigned short bytesper;      /* bytes per datum	 	*/
	FileType ftype;		      /* file type of output *fp	*/
} *outsegylist = (struct outsegyinfo *) NULL;

static FILE *lastfp = (FILE *) NULL;
static struct outsegyinfo *infoptr, **oldinfoptr;

static
void searchlist(FILE *fp)
{
	oldinfoptr = &outsegylist;
	for(infoptr = outsegylist; infoptr != ((struct outsegyinfo *) NULL);
	    infoptr = infoptr->nextinfo) {
		if(fp == infoptr->outfp) break;
		oldinfoptr = &infoptr->nextinfo;
	}
}

static
void datawrite(segy *tp, struct outsegyinfo *iptr, cwp_Bool fixed_length)
{
	unsigned int nstobewritten = fixed_length?iptr->nsfirst:tp->ns;
	unsigned int databytes = iptr->bytesper * nstobewritten;
	int nwritten = (int) efwrite((char *) (&((tp->data)[0])), 1, databytes,
			       iptr->outfp);

	if (nwritten != databytes)
		suerr("%s: on trace #%ld, tried to write %d bytes, "
		    "wrote %d bytes",
		    __FILE__, (infoptr->itr)+1, databytes, nwritten);

	return;
}


void fputtr_internal(FILE *fp, segy *tp, cwp_Bool fixed_length)
{
	/* search linked list for possible alternative */
	if(fp != lastfp)  searchlist(fp);

	if (infoptr == ((struct outsegyinfo *) NULL)) {
		/* initialize new segy output stream */

		/* allocate new segy output information table */
		*oldinfoptr = (struct outsegyinfo *)
			malloc(sizeof(struct outsegyinfo));
		infoptr = *oldinfoptr;
		infoptr->nextinfo = (struct outsegyinfo *) NULL;
		/* save FILE * ptr */
		infoptr->outfp = fp;
		infoptr->itr = 0;

		switch (infoptr->ftype = filestat(fileno(fp))) {
		case DIRECTORY:
			suerr("%s: segy output can't be a directory", __FILE__);
		case TTY:
			suerr("%s: segy output can't be tty", __FILE__);
		default:
			/* the rest are ok */
		break;
		}

		/* Sanity check the segy header */
		infoptr->nsfirst = tp->ns;
		if (infoptr->nsfirst > SU_NFLTS)
			suerr("%s: unable to handle %d > %d samples per trace",
			    __FILE__, infoptr->nsfirst, SU_NFLTS);
		switch(tp->trid) {
		case CHARPACK:
			infoptr->bytesper = sizeof(char); break;
		case SHORTPACK:
			infoptr->bytesper = 2*sizeof(char); break;
		default:
			infoptr->bytesper = sizeof(float); break;
		}

/*--------------------------------------------------------------------*\
   Write out a line header if it has been set as the default or has 
   requested on the caommandline.  Commandline takes precedence over
   the default in all cases.

   Reginald H. Beardsley                            rhb@acm.org
\*--------------------------------------------------------------------*/
              
/* commented out 15.05.09 WWS 

                getparint( "lheader" ,&out_line_hdr );

                if( out_line_hdr ){

                   if( in_line_hdr ){
                     (void) efwrite(&(su_text_hdr[0]), 1 ,3200 
                                    ,infoptr->outfp);

                   }else{
                     memset( su_text_hdr ,0 ,sizeof(su_text_hdr) );
                     sprintf( hdr_str ,"%-80s" 
                            ,"C 1 CLIENT CWP/SU default text header " );
                     strncat( su_text_hdr ,hdr_str ,80 );
                     for( i=1; i<40; i++ ){
                        sprintf( hdr_str ,"%-80s" ,"C" );
                        strncat( su_text_hdr ,hdr_str ,80 );
                     }
                     (void) efwrite(&(su_text_hdr[0]), 1 ,3200 
                                    ,infoptr->outfp);


                   }

                   memset( &su_binary_hdr ,0 ,sizeof(su_binary_hdr) );
                   su_binary_hdr.format = 1;
                   su_binary_hdr.hns = tp->ns;
                   su_binary_hdr.hdt = tp->dt;
 
 
                  (void) efwrite(&su_binary_hdr, 1
                              ,sizeof(su_binary_hdr), infoptr->outfp);
                }
*/
		
	}

	if (tp->ns != infoptr->nsfirst && fixed_length)
		suerr("%s: on trace #%ld, number of samples in header (%d) "
		    "differs from number for first trace (%d)", 
		    __FILE__, (infoptr->itr)+1, tp->ns, infoptr->nsfirst);
	

	(void) efwrite(tp, 1,HDRBYTES, infoptr->outfp);
	datawrite(tp, infoptr, fixed_length);
	
	++infoptr->itr;
	lastfp = infoptr->outfp;
}

void fputtr(FILE *fp, segy *tp)
{
 fputtr_internal(fp,tp,cwp_true);
}

void fvputtr(FILE *fp, segy *tp)
{
 fputtr_internal(fp,tp,cwp_false);
}


#endif



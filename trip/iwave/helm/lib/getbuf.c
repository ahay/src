/* getbuf.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;

/* GETBUF - returns pointers to beginning and end of buffer segmentc */

/* WWS, 23.9.92 */

/* Arguments: */

/* buffer      --->  buffer name */
/* pointer    <---   pointer assigned to beginning of segment */
/* length      --->  length of segment */
/* next       <--->  pointer to beginning of current segment on */
/*                   call, next segment on return */
/* last        --->  length of buffer */
/* ier        <--->  usual error flag */

/* --------------------------------------------------------------- */

/* Subroutine */ int getbuf_(char *buffer, integer *pointer, integer *length, 
	integer *next, integer *last, integer *ipdmp, integer *ier, ftnlen 
	buffer_len)
{
    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static integer idbg;

    /* Fortran I/O blocks */
    static cilist io___2 = { 0, 0, 0, 0, 0 };
    static cilist io___3 = { 0, 0, 0, 0, 0 };
    static cilist io___4 = { 0, 0, 0, 0, 0 };
    static cilist io___5 = { 0, 0, 0, 0, 0 };
    static cilist io___6 = { 0, 0, 0, 0, 0 };
    static cilist io___7 = { 0, 0, 0, 0, 0 };
    static cilist io___8 = { 0, 0, 0, 0, 0 };
    static cilist io___9 = { 0, 0, 0, 0, 0 };
    static cilist io___10 = { 0, 0, 0, 0, 0 };
    static cilist io___11 = { 0, 0, 0, 0, 0 };
    static cilist io___12 = { 0, 0, 0, 0, 0 };
    static cilist io___13 = { 0, 0, 0, 0, 0 };


/* debug flag */
    if (*ier != 0) {
	return 0;
    }
    idbg = 0;
    if (idbg != 0) {
	io___2.ciunit = *ipdmp;
	s_wsle(&io___2);
	do_lio(&c__9, &c__1, " GETBUF:", (ftnlen)8);
	e_wsle();
	io___3.ciunit = *ipdmp;
	s_wsle(&io___3);
	do_lio(&c__9, &c__1, " buffer  = ", (ftnlen)11);
	do_lio(&c__9, &c__1, buffer, buffer_len);
	e_wsle();
	io___4.ciunit = *ipdmp;
	s_wsle(&io___4);
	do_lio(&c__9, &c__1, " pointer = ", (ftnlen)11);
	do_lio(&c__3, &c__1, (char *)&(*pointer), (ftnlen)sizeof(integer));
	e_wsle();
	io___5.ciunit = *ipdmp;
	s_wsle(&io___5);
	do_lio(&c__9, &c__1, " length  = ", (ftnlen)11);
	do_lio(&c__3, &c__1, (char *)&(*length), (ftnlen)sizeof(integer));
	e_wsle();
	io___6.ciunit = *ipdmp;
	s_wsle(&io___6);
	do_lio(&c__9, &c__1, " next    = ", (ftnlen)11);
	do_lio(&c__3, &c__1, (char *)&(*next), (ftnlen)sizeof(integer));
	e_wsle();
	io___7.ciunit = *ipdmp;
	s_wsle(&io___7);
	do_lio(&c__9, &c__1, " last    = ", (ftnlen)11);
	do_lio(&c__3, &c__1, (char *)&(*last), (ftnlen)sizeof(integer));
	e_wsle();
    }
    *pointer = *next;
    *next += *length;
    if (*next > *last) {
	io___8.ciunit = *ipdmp;
	s_wsle(&io___8);
	do_lio(&c__9, &c__1, " Error: GETBUF", (ftnlen)14);
	e_wsle();
	io___9.ciunit = *ipdmp;
	s_wsle(&io___9);
	do_lio(&c__9, &c__1, " ran off end of buffer ", (ftnlen)23);
	do_lio(&c__9, &c__1, buffer, buffer_len);
	e_wsle();
	io___10.ciunit = *ipdmp;
	s_wsle(&io___10);
	do_lio(&c__9, &c__1, " pointer = ", (ftnlen)11);
	do_lio(&c__3, &c__1, (char *)&(*pointer), (ftnlen)sizeof(integer));
	e_wsle();
	io___11.ciunit = *ipdmp;
	s_wsle(&io___11);
	do_lio(&c__9, &c__1, " length  = ", (ftnlen)11);
	do_lio(&c__3, &c__1, (char *)&(*length), (ftnlen)sizeof(integer));
	e_wsle();
	io___12.ciunit = *ipdmp;
	s_wsle(&io___12);
	do_lio(&c__9, &c__1, " next    = ", (ftnlen)11);
	do_lio(&c__3, &c__1, (char *)&(*next), (ftnlen)sizeof(integer));
	e_wsle();
	io___13.ciunit = *ipdmp;
	s_wsle(&io___13);
	do_lio(&c__9, &c__1, " last    = ", (ftnlen)11);
	do_lio(&c__3, &c__1, (char *)&(*last), (ftnlen)sizeof(integer));
	e_wsle();
	*ier = 99;
	return 0;
    }
    return 0;
} /* getbuf_ */


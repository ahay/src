/* Copyright (c) Colorado School of Mines, 2010.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/***************************************************************************** 
SWAPBYTE - Functions to SWAP the BYTE order of binary data 

swap_short_2		swap a short integer
swap_u_short_2		swap an unsigned short integer
swap_int_4		swap a 4 byte integer
swap_u_int_4		swap an unsigned integer
swap_long_4		swap a long integer
swap_u_long_4		swap an unsigned long integer
swap_float_4		swap a float
swap_double_8		swap a double

******************************************************************************
Function Prototypes:
void swap_short_2(short *tni2);
void swap_u_short_2(unsigned short *tni2);
void swap_int_4(int *tni4);
void swap_u_int_4(unsigned int *tni4);
void swap_long_4(long *tni4);
void swap_u_long_4(unsigned long *tni4);
void swap_float_4(float *tnf4);
void swap_double_8(double *tndd8);

******************************************************************************
Notes:
These routines are necessary for reversing the byte order of binary data
for transportation between big-endian and little-endian machines. Examples
of big-endian machines are IBM RS6000, SUN, NeXT. Examples of little
endian machines are PC's and DEC.

These routines have been tested with PC data and run on PC's running
several PC versions of UNIX, but have not been tested on DEC.

Also, the number appended to the name of the routine refers to the
number of bytes that the item is assumed to be.

******************************************************************************
Authors: Jens Hartmann,   Institut fur Geophysik, Hamburg, Jun 1993
	 John Stockwell, CWP, Colorado School of Mines, Jan 1994
***************************************************************************/
/**************** end self doc ********************************/

void swap_short_2(short *tni2)
/*<*************************************************************************
swap_short_2		swap a short integer
**************************************************************************>*/
{
 *tni2=(((*tni2>>8)&0xff) | ((*tni2&0xff)<<8));  
}

void swap_u_short_2(unsigned short *tni2)
/*<*************************************************************************
swap_u_short_2		swap an unsigned short integer
**************************************************************************>*/
{
 *tni2=(((*tni2>>8)&0xff) | ((*tni2&0xff)<<8));  
}

void swap_int_4(int *tni4)
/*<*************************************************************************
swap_int_4		swap a 4 byte integer
**************************************************************************>*/
{
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |
	    ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));  
}

void swap_u_int_4(unsigned int *tni4)
/*<*************************************************************************
swap_u_int_4		swap an unsigned integer
**************************************************************************>*/
{
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |
	    ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));  
}

void swap_long_4(long *tni4)
/*<*************************************************************************
swap_long_4		swap a long integer
**************************************************************************>*/
{
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |
	    ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));  
}

void swap_u_long_4(unsigned long *tni4)
/*<*************************************************************************
swap_u_long_4		swap an unsigned long integer
**************************************************************************>*/
{
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |
	    ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));  
}

void swap_float_4(float *tnf4)
/*<*************************************************************************
swap_float_4		swap a float
**************************************************************************>*/
{
 int *tni4=(int *)tnf4;
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |
	    ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));  
}

void swap_double_8(double *tndd8)
/*<*************************************************************************
swap_double_8		swap a double
**************************************************************************>*/
{
  char *tnd8=(char *)tndd8;
  char tnc;

  tnc= *tnd8;
  *tnd8= *(tnd8+7);
  *(tnd8+7)=tnc;

  tnc= *(tnd8+1);
  *(tnd8+1)= *(tnd8+6);
  *(tnd8+6)=tnc;

  tnc= *(tnd8+2);
  *(tnd8+2)= *(tnd8+5);
  *(tnd8+5)=tnc;

  tnc= *(tnd8+3);
  *(tnd8+3)= *(tnd8+4);
  *(tnd8+4)=tnc;
}

#ifdef TEST
#define N 0251    /* this is a made-up octal number */
main()
{
	short xs=N;
	unsigned short xus=N;
	long xl=N;
	unsigned long xul=N;
	float xf=N;
	double xd=N;

	/* swap short */
	fprintf (stdout, "Output is in octal\n\n");
	fprintf (stdout, "short integer:\n");
	fprintf (stdout, "unswapped xs = %o\n", xs);  
	swap_short_2(&xs);
	fprintf (stdout, "swapped xs = %o\n", xs);  
	swap_short_2(&xs);
	fprintf (stdout, "swapped back xs = %o\n\n", xs);  

	/* swap u_short */
	fprintf (stdout, "unsigned short integer:\n");
	fprintf (stdout, "unswapped xus = %o\n", xus);
	swap_u_short_2(&xus);
	fprintf (stdout, "swapped xus = %o\n", xus);  
	swap_u_short_2(&xus);
	fprintf (stdout, "swapped back xus = %o\n\n", xus);

	/* swap long */
	fprintf (stdout, "long integer:\n");
	fprintf (stdout, "unswapped xl = %o\n", xl);
	swap_long_4(&xl);
	fprintf (stdout, "swapped xl = %o\n", xl);
	swap_long_4(&xl);
	fprintf (stdout, "swapped back xl = %o\n\n", xl);

	/* swap u_long */
	fprintf (stdout, "unsigned long integer:\n");
	fprintf (stdout, "unswapped xul = %o\n", xul);
	swap_u_long_4(&xul);
	fprintf (stdout, "swapped xul = %o\n", xul);
	swap_u_long_4(&xul);
	fprintf (stdout, "swapped back xul = %o\n\n", xul);

	/* swap float */
	fprintf (stdout, "float:\n");
	fprintf (stdout, "unswapped xf = %o\n", xf);
	swap_float_4(&xf);
	fprintf (stdout, "swapped xf = %o\n", xf);
	swap_float_4(&xf);
	fprintf (stdout, "swapped back xf = %o\n\n", xf);

	/* swap double */
	fprintf (stdout, "double:\n");
	fprintf (stdout, "unswapped xd = %o\n", xd);
	swap_double_8(&xd);
	fprintf (stdout, "swapped xd = %o\n", (float) xd);
	swap_double_8(&xd);
	fprintf (stdout, "swapped back xd = %o\n", xd);
}
#endif

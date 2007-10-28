/* Copyright (c) Colorado School of Mines, 2007.*/
/* All rights reserved.                       */

#ifdef __GNUC__
#	include <time.h>
#else
#	include <sys/time.h>
#endif

float wallsec(void)
/*< return elapsed time (wall clock time) in seconds >*/
/******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 04/29/89
*****************************************************************************/
{
#ifdef __GNUC__
	static int    firsttime = 1;
	static time_t firstsec     ;

	if ( firsttime ) {
		firsttime = 0;
		firstsec  = time ((time_t *)NULL);
		return (0.0);
	} else
		return ((float)(time((time_t *)NULL) - firstsec));
#else
	struct timeval tp;
	struct timezone tzp;
	static int firsttime=1;
	static long firstsec,firstusec;
	long sec,usec;

	gettimeofday(&tp,&tzp);
	if (firsttime) {
		firsttime=0;
		firstsec = tp.tv_sec;
		firstusec = tp.tv_usec;
		return(0.0);
	} else {
		sec = tp.tv_sec-firstsec;
		usec = tp.tv_usec-firstusec;
		return((float)sec+1.0e-6*(float)usec);
	}
#endif
}

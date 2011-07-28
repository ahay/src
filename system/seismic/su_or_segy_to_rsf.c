/* Procedures used by segy2rsf and su2rsf. Excerpted from the code of segyread.
*/
/*
  Copyright (C) 2004 University of Texas at Austin

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE
#endif
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <rsf.h>

#include "segy.h"
/*^*/

off_t
readfilesize (FILE * file)
/*< Returns file size in bytes. Cannot be used on stdin. >*/
{
  off_t pos;
  extern int fseeko (FILE * stream, off_t offset, int whence);
  extern off_t ftello (FILE * stream);

  fseeko (file, 0, SEEK_END);
  pos = ftello (file);		/* pos is the filesize in bytes */
  fseeko (file, 0, SEEK_SET);

  return pos;
}

void
su_or_segy_to_rsf (bool verb, bool su, int ntr, int format, int ns,
		   int itrace[SF_NKEYS], off_t nsegy, FILE * file, float dt)
/*< Common block of code for both segy2rsf and su2rsf >*/
{
  /* 
     ntr = number of traces in input file
     format = input parameter for SEG-Y data (1=IBM float, 5=IEEE float, etc)
     ns = number of samples in one trace
     nsegy = number of bytes in an input trace plus its header
     file = handle to input
     dt = time sampling of traces
   */

  const char *headname;		/* holds name of file with trace headers */
  char *trace;
  const char *read;		/* for similarly-named input parameter */
  sf_file out;			/* tag for stdout */
  sf_file hdr;			/* tag for file with trace headers */
  sf_file msk = NULL;		/* tag of mask file */
  int itr;			/* index to count of traces */
  int n2;			/* number of traces that pass through mask */
  int *mask;			/* array holding mask values */
  float *ftrace;		/* for writing out data coming from SEG-Y file */

  /* Why declare it by hand? It is included in stdio.h */
  extern int fseeko (FILE * stream, off_t offset, int whence);

  if (verb)
    sf_warning ("Expect %d traces", ntr);

  if (NULL == (read = sf_getstring ("read")))
    read = "b";
  /* what to read: h - header, d - data, b - both (default) */

  if (NULL != sf_getstring ("mask"))
    {
      /* optional header mask for reading only selected traces */
      msk = sf_input ("mask");
      if (SF_INT != sf_gettype (msk))
	sf_error ("Need integer mask");

      mask = sf_intalloc (ntr);
      sf_intread (mask, ntr, msk);
      sf_fileclose (msk);

      for (n2 = itr = 0; itr < ntr; itr++)
	{
	  if (mask[itr])
	    n2++;
	}
    }
  else
    {
      mask = NULL;
      n2 = ntr;
    }

  if (read[0] != 'h')
    {				/* not only header */
      out = sf_output ("out");
      sf_putint (out, "n1", ns);
      sf_putint (out, "n2", n2);
      sf_putfloat (out, "d1", dt);
      sf_putfloat (out, "o1", 0.);
      sf_setformat (out, "native_float");
      ftrace = su ? NULL : sf_floatalloc (ns);
    }
  else
    {
      out = NULL;
      ftrace = NULL;
    }

  if (read[0] != 'd')
    {				/* not only data */
      hdr = sf_output ("tfile");
      sf_putint (hdr, "n1", SF_NKEYS);
      sf_putint (hdr, "n2", n2);
      sf_setformat (hdr, "native_int");

      if (NULL == (headname = sf_getstring ("tfile")))
	headname = "tfile";
      /* output trace header file */
      if (NULL != out)
	sf_putstring (out, "head", headname);
    }
  else
    {
      hdr = NULL;
    }

  if (NULL != out)
    sf_fileflush (out, NULL);

  switch (read[0])
    {
    case 'h':			/* header only */
      trace = sf_charalloc (SF_HDRBYTES);
      nsegy -= SF_HDRBYTES;

      for (itr = 0; itr < ntr; itr++)
	{
	  if (NULL != mask && !mask[itr])
	    {
	      fseeko (file, SF_HDRBYTES, SEEK_CUR);
	      continue;
	    }
	  else if (SF_HDRBYTES != fread (trace, 1, SF_HDRBYTES, file))
	    {
	      sf_error ("Error reading trace header %d", itr + 1);
	    }
	  fseeko (file, nsegy, SEEK_CUR);

	  segy2head (trace, itrace, SF_NKEYS);
	  sf_intwrite (itrace, SF_NKEYS, hdr);
	}

      break;
    case 'd':			/* data only */
      nsegy -= SF_HDRBYTES;
      trace = sf_charalloc (nsegy);

      for (itr = 0; itr < ntr; itr++)
	{
	  fseek (file, SF_HDRBYTES, SEEK_CUR);
	  if (NULL != mask && !mask[itr])
	    {
	      fseeko (file, nsegy, SEEK_CUR);
	      continue;
	    }
	  else if (nsegy != fread (trace, 1, nsegy, file))
	    {
	      sf_error ("Error reading trace data %d", itr + 1);
	    }

	  if (su)
	    {
	      sf_charwrite (trace, ns * sizeof (float), out);
	    }
	  else
	    {
	      segy2trace (trace, ftrace, ns, format);
	      sf_floatwrite (ftrace, ns, out);
	    }
	}

      break;
    default:			/* both header and data */
      trace = sf_charalloc (nsegy);

      for (itr = 0; itr < ntr; itr++)
	{
	  if (NULL != mask && !mask[itr])
	    {
	      fseeko (file, nsegy, SEEK_CUR);
	      continue;
	    }
	  else if (nsegy != fread (trace, 1, nsegy, file))
	    {
	      sf_error ("Error reading trace header %d", itr + 1);
	    }

	  segy2head (trace, itrace, SF_NKEYS);
	  sf_intwrite (itrace, SF_NKEYS, hdr);

	  if (su)
	    {
	      sf_charwrite (trace + SF_HDRBYTES, ns * sizeof (float), out);
	    }
	  else
	    {
	      segy2trace (trace + SF_HDRBYTES, ftrace, ns, format);
	      sf_floatwrite (ftrace, ns, out);
	    }
	}
      break;
    }


  exit (0);
}

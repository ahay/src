/* Convert a SEG-Y dataset to RSF.
Same as sfsegyread su=n. Split from sfsuread to decrease the number of
options of sfsuread and to improve readability 
of both codes. See sfsegyread help for more info.*/
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
#include "su_or_segy_to_rsf.h"

int
main (int argc, char *argv[])
{
  const bool su = false;
  bool verb, xdr;
  char ahead[SF_EBCBYTES], bhead[SF_BNYBYTES];
  char *headname, *filename;
  int format, ns, ntr, itrace[SF_NKEYS];
  off_t pos, nsegy;
  FILE *head, *file;
  float dt;

  /**** Start common code with su2rsf (16 lines) ****/
  sf_init (argc, argv);

  if (!sf_getbool ("verb", &verb))
    verb = false;
  /* Verbosity flag */

  if (!sf_getbool ("endian", &xdr))
    xdr = true;
  /* Whether to automatically estimate endianness or not */
  if (xdr)
    endian ();

  if (NULL == (filename = sf_getstring ("tape")))	/* input data */
    sf_error ("Need to specify tape=");

  if (NULL == (file = fopen (filename, "rb")))
    sf_error ("Cannot open \"%s\" for reading:", filename);

  pos = readfilesize (file);
/**** End common code with su2rsf ****/
  if (SF_EBCBYTES != fread (ahead, 1, SF_EBCBYTES, file))
    sf_error ("Error reading ebcdic header");

  ebc2asc (SF_EBCBYTES, ahead);

  if (NULL == (headname = sf_getstring ("hfile")))
    headname = "header";
  /* output text data header file */

  if (NULL == (head = fopen (headname, "w")))
    sf_error ("Cannot open file \"%s\" for writing ascii header:", headname);

  if (SF_EBCBYTES != fwrite (ahead, 1, SF_EBCBYTES, head))
    sf_error ("Error writing ascii header");
  fclose (head);

  if (verb)
    sf_warning ("ASCII header written to \"%s\"", headname);

  if (SF_BNYBYTES != fread (bhead, 1, SF_BNYBYTES, file))
    sf_error ("Error reading binary header");

  if (NULL == (headname = sf_getstring ("bfile")))
    headname = "binary";
  /* output binary data header file */

  if (NULL == (head = fopen (headname, "wb")))
    sf_error ("Cannot open file \"%s\" for writing binary header:", headname);

  if (SF_BNYBYTES != fwrite (bhead, 1, SF_BNYBYTES, head))
    sf_error ("Error writing binary header");
  fclose (head);

  if (verb)
    sf_warning ("Binary header written to \"%s\"", headname);

  if (!sf_getint ("format", &format))
    format = segyformat (bhead);
  /* [1,2,3,5] Data format. The default is taken from binary header.
     1 is IBM floating point
     2 is 4-byte integer
     3 is 2-byte integer
     5 is IEEE floating point
   */

  switch (format)
    {
    case 1:
      if (verb)
	sf_warning ("Assuming IBM floating point format");
      break;
    case 2:
      if (verb)
	sf_warning ("Assuming 4 byte integer format");
      break;
    case 3:
      if (verb)
	sf_warning ("Assuming 2 byte integer format");
      break;
    case 5:
      if (verb)
	sf_warning ("Assuming IEEE floating point format");
      break;
    default:
      sf_error ("Nonstandard format: %d", format);
      break;
    }

  if (!sf_getint ("ns", &ns))
    ns = segyns (bhead);
  /* Number of samples. The default is taken from binary header */
  if (0 >= ns)
    sf_error ("Number of samples is not set in binary header");

  if (verb)
    sf_warning ("Detected trace length of %d", ns);

  dt = segydt (bhead);
  nsegy = SF_HDRBYTES + ((3 == format) ? ns * 2 : ns * 4);
  ntr = (pos - SF_EBCBYTES - SF_BNYBYTES) / nsegy;

  su_or_segy_to_rsf (verb, su, ntr, format, ns, itrace, nsegy, file, dt);

  exit (0);
}

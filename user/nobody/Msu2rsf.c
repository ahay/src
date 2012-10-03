/* Convert a SU dataset to RSF.
Same as sfsegyread su=n. Split from sfsegyread to decrease the number of
options of sfsuread, and to improve readability 
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
  const bool su = true;
  bool verb, xdr;
  char *filename, *trace;
  int format=0, ns, ntr, itrace[SF_NKEYS];
  off_t pos, nsegy;
  float dt;
  FILE *file;
  extern int fseeko (FILE * stream, off_t offset, int whence);

/**** Start common code with segy2rsf (16 lines) ****/
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
/**** End common code with segy2rsf ****/

  /* Figure out ns and ntr */

  trace = sf_charalloc (SF_HDRBYTES);
  if (SF_HDRBYTES != fread (trace, 1, SF_HDRBYTES, file))
    sf_error ("Error reading first trace header");
  fseeko (file, 0, SEEK_SET);

  segy2head (trace, itrace, SF_NKEYS);
  ns = itrace[segykey ("ns")];
  dt = itrace[segykey ("dt")] / 1000000.;
  free (trace);

  nsegy = SF_HDRBYTES + ns * 4;
  ntr = pos / nsegy;

  su_or_segy_to_rsf (verb, su, ntr, format, ns, itrace, nsegy, file, dt);

  exit (0);
}

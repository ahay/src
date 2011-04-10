/* Remove RSF files together with their data.

Takes: file1.rsf [file2.rsf ...] [-i] [-v] [-f] 

Mimics the standard Unix rm command.

See also: sfmv, sfcp.
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

#include <string.h>

#include <rsf.h>

int main (int argc, char *argv[])
{
    int i;
    char *arg;
    bool force = false, verb = false, inquire = false;
   
    sf_init(argc,argv);

    for (i=1; i < argc; i++) {
	arg = argv[i];
	if ('-' == arg[0]) { /* it is an option */
	    if (NULL != strchr(arg,'f')) force = true;
	    if (NULL != strchr(arg,'v')) verb = true;
	    if (NULL != strchr(arg,'i')) inquire = true;
	} else { /* it is a file */
	    sf_rm(arg, force, verb, inquire);
	}
    }

    exit (0);
}

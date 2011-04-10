/* Copy or move a dataset.

Takes: in.rsf out.rsf

sfcp - copy, sfmv - move.
Mimics standard Unix commands.
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

int main(int argc, char* argv[])
{
    int i;
    char *infile=NULL, *prog;
    sf_file in, out;

    sf_init (argc,argv);

    prog = sf_getprog();

    if (sf_stdin() && NULL != strstr (prog,"cp")) { 
	/* cp with input file in stdin */
	in = sf_input("in");
	out = sf_output("out");
    } else {
	in = NULL;
	out = NULL;

	/* the first two non-parameters are in and out files */
	for (i=1; i< argc; i++) {
	    if (NULL == strchr(argv[i],'=')) {
		if (NULL == in) {
		    infile = argv[i];
		    in = sf_input (infile);
		} else {
		    out = sf_output (argv[i]);
		    break;
		}
	    }
	}

	if (NULL == in || NULL == out)
	    sf_error ("not enough input");
    }

    sf_setformat(out,sf_histstring(in,"data_format"));

    sf_cp(in,out);
    if (NULL != strstr (prog,"mv"))
	sf_rm(infile,false,false,false);

    exit (0);
}

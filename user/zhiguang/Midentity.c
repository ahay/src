/* Copy dataset.
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
    bool adj;
    sf_file in, out;

    sf_init (argc,argv);
    if(!sf_getbool("adj", &adj)) adj=true;
    in = sf_input("in");
    out = sf_output("out");

    sf_setformat(out,sf_histstring(in,"data_format"));

    sf_cp(in,out);

    exit (0);
}

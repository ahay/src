/* Creates just the ascii header from parameters
Wrapper for sf_fileflush (copy RSF header of a file to another) */
/*
  Copyright (C) 2007 Ioan Vlad

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

#include <rsf.h>

int main (int argc, char*argv[]) {

    sf_file in=NULL;
    sf_file out=NULL;

    sf_init(argc, argv);

    in  = sf_input("in");
    out = sf_output("out");

    sf_fileflush(out,in);


    exit(0);
}

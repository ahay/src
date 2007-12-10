/* Creates just the ascii header from parameters
Wrapper for sf_fileflush (creating RSF header from params) */
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

    int i, n;
    char key[4];
    sf_file out;

    sf_init(argc, argv);
    out = sf_output("out");

    /* After the model in sfspike: */
    for (i=0; i < SF_MAX_DIM; i++) {
        snprintf(key,3,"n%d",i+1);
    if (!sf_getint(key,&n)) break;
    /*( n# size of #-th axis )*/  
        sf_putint(out,key,n);
    }

    if (0==i) sf_error("Need n1=");

    sf_setformat(out, "native_float");

    sf_fileflush(out, NULL);

    exit(0);
}

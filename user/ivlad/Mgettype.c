/* Displays numerical type of a dataset
Output can be be: SF_CHAR, SF_COMPLEX, SF_FLOAT, SF_INT, SF_SHORT, SF_UCHAR
Shell/tester for sf_gettype*/
/*
  Copyright (C) 2004 University of Texas at Austin
  Copyright (C) 2009 Ioan Vlad

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

int main(int argc, char* argv[])
{
    sf_file in=NULL;
    sf_datatype type;

    sf_init (argc,argv);

    in = sf_input("in");

    type = sf_gettype (in);

    switch (type) {
        case SF_FLOAT:
            printf("SF_FLOAT\n");
            break;
        case SF_INT:
            printf("SF_INT\n");
            break;
        case SF_SHORT:
            printf("SF_SHORT\n");
            break;
        case SF_COMPLEX:
            printf("SF_COMPLEX\n");
            break;
        case SF_UCHAR:
            printf("SF_UCHAR\n");
            break;
        case SF_CHAR:
            printf("SF_CHAR\n");
            break;
        case SF_DOUBLE:
            printf("SF_DOUBLE\n");
            break;
        default:
            printf("Unknown!\n");
            break;
    }

    exit (0);
}

/* Extend a dataset by duplicating in the specified axis dimension.
   This operation is adjoint to sfstack.
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
#include <stdio.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int j, n1, n2, i2, ni, axis, esize;
    float f;
    off_t i3, n3;
    size_t n;
    sf_file in, out;
    char key1[12], key2[12], *val, *trace;

    sf_init (argc, argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"esize",&esize)) sf_error("No esize= in input");

    if (!sf_getint("axis",&axis)) axis=2;
    /* which axis to spray */

    n = (size_t) esize;
    for (j=0; j < axis-1; j++) {
	snprintf(key2,12,"n%d",j+1);
	if (!sf_histint(in,key2,&ni)) break;
	n *= ni;
    }

    if (!sf_getint("n",&n2)) sf_error("Need n=");
    /* Size of the newly created dimension */    
    snprintf(key1,12,"n%d",axis);
    sf_putint(out,key1,n2);

    if (sf_getfloat("d",&f)) {
	/* Sampling of the newly created dimension */ 
	snprintf(key1,4,"d%d",axis);
	sf_putfloat(out,key1,f);
    }
    
    if (sf_getfloat("o",&f)) {
	/* Origin of the newly created dimension */
	snprintf(key1,4,"o%d",axis);
	sf_putfloat(out,key1,f);
    }

    if (NULL != (val = sf_getstring("label"))) {
	/* Label of the newly created dimension */
	snprintf(key1,8,"label%d",axis);
	sf_putstring(out,key1,val);
    }

    if (NULL != (val = sf_getstring("unit"))) {
	/* Units of the newly created dimension */
	snprintf(key1,7,"unit%d",axis);
	sf_putstring(out,key1,val);
    }

    n3 = sf_shiftdim(in, out, axis);
    
    sf_fileflush(out,in);
    sf_setform(in,SF_NATIVE);
    sf_setform(out,SF_NATIVE);

    trace = sf_charalloc (n);
    
    for (i3=0; i3 < n3; i3++) {
	sf_charread (trace, n, in);
	for (i2=0; i2 < n2; i2++) {
	    sf_charwrite(trace, n, out);
	} 
    }


    exit (0);
}

/* 	$Id$	 */

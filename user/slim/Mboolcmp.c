/* Element-wise boolean comparison of values. For int/float/complex data-sets.

Written by: C. Brown, UBC
Created: Nov 2007
*/
/*
  Copyright (C) 2006 The University of British Columbia at Vancouver
  
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

int main(int argc, char *argv[])
{
  int n[SF_MAX_DIM],n_r[SF_MAX_DIM],*qq;
    size_t nsiz,nsiz_r,i,dim,dim_r;
    float eps,fl,fr;
    char *right,*sign,*pl,*pr;
    sf_complex c;
    sf_file in,in_r,out;
    sf_datatype type;

    sf_init(argc,argv);

    if (NULL == (right=sf_getstring("right"))) sf_error("No right parameter set.");
    /* the rsf file you will be comparing to */

    if (NULL == (sign=sf_getstring("sign"))) sign="eq";
    /* 'eq'(default),'gt','ge','lq','lt','ne'
        sign=   'eq' equal-to ( == )
        sign=   'gt' greater-than ( > )
        sign=   'ge' greater-than or equal-to ( >= )
        sign=   'lq' less-than or equal-to ( <= )
        sign=   'lt' less-than ( < )
        sign=   'ne' not-equal ( != )
    */ 

    if (!sf_getfloat("eps",&eps)) eps=0;
    /* comparing within this range epsilon */

    in = sf_input("in");
    in_r = sf_input(right);
    out = sf_output("out");

    dim = (size_t) sf_filedims(in,n);
    for (nsiz=1, i=0; i < dim; i++) nsiz *= n[i];

    dim_r = (size_t) sf_filedims(in_r,n_r);
    for (nsiz_r=1, i=0; i < dim_r; i++) nsiz_r *= n_r[i];

    type = sf_gettype(in);

    if (type != sf_gettype(in_r)) sf_error("Type of input and right files do not match.");
    if (nsiz != nsiz_r) sf_error("Size of input and right files do not match.");


    switch (type) {
      case SF_FLOAT:
	pl = (char*)sf_floatalloc(nsiz);
	pr = (char*)sf_floatalloc(nsiz);
	sf_floatread((float*) pl,nsiz,in);
	sf_floatread((float*) pr,nsiz,in_r);
	break;
      case SF_INT:
	pl = (char*)sf_intalloc(nsiz);
	pr = (char*)sf_intalloc(nsiz);
	sf_intread((int*) pl,nsiz,in);
	sf_intread((int*) pr,nsiz,in_r);
	break;
      case SF_COMPLEX:
	pl = (char*)sf_complexalloc(nsiz);
	pr = (char*)sf_complexalloc(nsiz);
	sf_complexread((sf_complex*) pl,nsiz,in);
	sf_complexread((sf_complex*) pr,nsiz,in_r);
	break;
      default:
	sf_error("Type not understood.");
	break;
    }
    qq = sf_intalloc(nsiz);

    for (i=0; i<nsiz; i++) {
      switch (type) {
        case SF_FLOAT:
	  fl = ((float*)pl)[i];
	  fr = ((float*)pr)[i];
	  break;
        case SF_INT:
	  fl = (float) ((int*)pl)[i];
	  fr = (float) ((int*)pr)[i];
	  break;
        case SF_COMPLEX:
	  c=((sf_complex*)pl)[i];  
	  fl=cabsf(c);
	  c=((sf_complex*)pr)[i];  
	  fr=cabsf(c);
	  break;
        default:
	  sf_error("Type not understood.");
	  break;
      }

      if      (0==strcmp(sign,"ge")) qq[i] = ((fl-fr) >= -eps);
      else if (0==strcmp(sign,"gt")) qq[i] = ((fl-fr) > -eps);
      else if (0==strcmp(sign,"eq")) qq[i] = (abs(fl-fr) <= eps);
      else if (0==strcmp(sign,"lt")) qq[i] = ((fl-fr) < eps);
      else if (0==strcmp(sign,"lq")) qq[i] = ((fl-fr) <= eps);
      else if (0==strcmp(sign,"ne")) qq[i] = (abs(fl-fr) > eps);
      else sf_error("Sign not recognized. Please specify: gt,ge,eq,lq,lt,ne");
    }

    sf_settype(out,SF_INT);
    sf_intwrite(qq,nsiz,out);
    exit(0);
}

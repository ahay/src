/* Element-wise boolean comparison of values. For int/float/complex data-sets.
This program will solve the solution to this problem:
    - [input] [sign] [right]
    - sfboolcmp <left.rsf sign=ge right=right.rsf 
    - left.rsf >= right.rsf
This will return a vector of same length as left and return 0's or 1's depending on the result of the inequality.  Optionally you can supply a right_f parameter to compare the input data to a single value.

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
    off_t n[SF_MAX_DIM],n_r[SF_MAX_DIM], nsiz,nsiz_r=0,nleft;
    int qq[BUFSIZ];
    char buf[BUFSIZ],buf_r[BUFSIZ],*right=0,*sign;
    float eps,fl=0,fr;
    size_t bufsiz=BUFSIZ,dim,dim_r,i,nbuf;
    sf_complex c;
    sf_file in,in_r=0,out;
    sf_datatype type;
    bool cmp_num=false;

    sf_init(argc,argv);

    cmp_num = sf_getfloat("right_f",&fr);
    /* compare input (left) to a single float value (right) */

    if (!cmp_num && NULL == (right=sf_getstring("right"))) sf_error("No right or right_f parameter set.");
    /* the rsf file you will be comparing to */

    if (NULL == (sign=sf_getstring("sign"))) sign="eq";
    /* 'eq'(default),'gt','ge','lq','lt','ne'
        sign=   'eq' equal-to ( == )
        sign=   'gt' greater-than ( > )
        sign=   'ge' greater-than or equal-to ( >= )
        sign=   'lq' less-than or equal-to ( <= )
        sign=   'lt' less-than ( < )
        sign=   'ne' not-equal ( != )
	sign=   'and' the values are both non-zero ( && )
	sign=   'or' one value is non-zero ( !! )
    */ 

    if (!sf_getfloat("eps",&eps)) eps=0;
    /* comparing within this range epsilon */

    in = sf_input("in");
    out = sf_output("out");
    sf_settype(out,SF_INT);

    dim = sf_largefiledims(in,n);
    for (nsiz=1, i=0; i < dim; i++) nsiz *= n[i];

    if (!cmp_num) {
      in_r = sf_input(right);
      dim_r = (size_t) sf_largefiledims(in_r,n_r);
      for (nsiz_r=1, i=0; i < dim_r; i++) nsiz_r *= n_r[i];
    }

    bufsiz /= sf_esize(in);
    type = sf_gettype(in);

    if (!cmp_num && type != sf_gettype(in_r)) sf_error("Type of input and right files do not match.");
    if (!cmp_num && nsiz != nsiz_r) sf_error("Size of input and right files do not match.");


    for (nleft=nsiz;nleft>0;nleft -= nbuf) {
      nbuf = (bufsiz < nleft)? bufsiz: nleft;
      switch (type) {
        case SF_FLOAT:
	  sf_floatread((float*) buf,nbuf,in);
	  if (!cmp_num) sf_floatread((float*) buf_r,nbuf,in_r);
	  break;
        case SF_INT:
	  sf_intread((int*) buf,nbuf,in);
	  if (!cmp_num) sf_intread((int*) buf_r,nbuf,in_r);
	  break;
        case SF_COMPLEX:
	  sf_complexread((sf_complex*) buf,nbuf,in);
	  if (!cmp_num) sf_complexread((sf_complex*) buf_r,nbuf,in_r);
	  break;
        default:
	  sf_error("Type not understood.");
	  break;
      }
      for (i=0; i<nbuf; i++) {
	switch (type) {
          case SF_FLOAT:
	    fl = ((float*)buf)[i];
	    if (!cmp_num) fr = ((float*)buf_r)[i];
	    break;
          case SF_INT:
	    fl = (float) ((int*)buf)[i];
	    if (!cmp_num) fr = (float) ((int*)buf_r)[i];
	    break;
          case SF_COMPLEX:
	    c=((sf_complex*)buf)[i];  
	    fl=cabsf(c);
	    if (!cmp_num) {
	      c=((sf_complex*)buf_r)[i];  
	      fr=cabsf(c);
	    }
	    break;
          default:
	    sf_error("Type not understood.");
	    break;
	}

	if      (0==strcmp(sign,"ge")) qq[i] = ((fl-fr) >= -eps);
	else if (0==strcmp(sign,"gt")) qq[i] = ((fl-fr) > -eps);
	else if (0==strcmp(sign,"eq")) qq[i] = (fabs(fl-fr) <= eps);
	else if (0==strcmp(sign,"lt")) qq[i] = ((fl-fr) < eps);
	else if (0==strcmp(sign,"lq")) qq[i] = ((fl-fr) <= eps);
	else if (0==strcmp(sign,"ne")) qq[i] = (fabs(fl-fr) > eps);
	else if (0==strcmp(sign,"and")) qq[i] = ((fabs(fl) > eps) && (fabs(fr) > eps));
	else if (0==strcmp(sign,"or")) qq[i] = ((fabs(fl) > eps) || (fabs(fr) > eps));
	else sf_error("Sign not recognized. Please specify: gt,ge,eq,lq,lt,ne,and,or");
      }
      sf_intwrite(qq,nbuf,out);
    }

    exit(0);
}

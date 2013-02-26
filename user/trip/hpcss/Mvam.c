/* Create a layered model. */

/*************************************************************************

Copyright Rice University, 2009.
All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, provided that the above copyright notice(s) and this
permission notice appear in all copies of the Software and that both the
above copyright notice(s) and this permission notice appear in supporting
documentation.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT OF THIRD PARTY
RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR HOLDERS INCLUDED IN THIS
NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL INDIRECT OR CONSEQUENTIAL
DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.

Except as contained in this notice, the name of a copyright holder shall
not be used in advertising or otherwise to promote the sale, use or other
dealings in this Software without prior written authorization of the
copyright holder.

**************************************************************************/

/* Modified for distribution with Madagascar */

#include <rsf.h>

int main(int argc, char ** argv){
  int nz,nx;
  float dz,dx;
  int ix;
  int iz;
  float *v;
  sf_file vfile;

  sf_init(argc,argv);
  vfile = sf_output("out");
  sf_setformat(vfile,"native_float");

  if (!sf_getint("nz",&nz)) sf_error("Need nz=");
  /* depth grid */
  if (!sf_getint("nx",&nx)) sf_error("Need nx=");
  /* distance grid */
  if (!sf_getfloat("dz",&dz)) sf_error("Need dz=");
  /* depth sampling */
  if (!sf_getfloat("dx",&dx)) sf_error("Need dx=");
  /* distance sampling */
  
  sf_putint(vfile,"n1",nz);
  sf_putint(vfile,"n2",nx);
  sf_putfloat(vfile,"d1",dz);
  sf_putfloat(vfile,"d2",dx);
  sf_putfloat(vfile,"o1",0.);
  sf_putfloat(vfile,"o2",0.);
  sf_putstring(vfile,"label1","Depth");
  sf_putstring(vfile,"label2","Offset");
  sf_putstring(vfile,"unit1","km");
  sf_putstring(vfile,"unit2","km");

  sf_putstring(vfile,"label","Velocity");
  sf_putstring(vfile,"unit","km/s");

  v = sf_floatalloc(nz);
  
  for (ix=0;ix<nx;ix++){
      iz=0;
      while(iz<nz/3){
	  v[iz]=1.5;
	  iz++;
      }
      while(iz<nz*2/3){
	  v[iz]=1.6;
	  iz++;
      }
      while(iz<nz){
	  v[iz]=1.55;
	  iz++;
      }
      sf_floatwrite(v,nz,vfile);
  }
  
  exit(0);
}

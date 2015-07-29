/* Rice HPCSS seismic modeling and migration. */

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

/*
#define VERB 
*/

void step24_forward(float * tgt_p, 
		    float * src_p, 
		    float * v,
		    int nz, int nx,
		    float rz, float rx, float s) 
/*< step forward using 2,4 scheme >*/
{

  int iz;
  int ix;
  int ioff; 
  float two=2.0;
  float mot=-1.0/12.0;
  float fth=4.0/3.0;
  float mfh=-5.0/4.0;

  float motz=rz*mot;
  float motx=rx*mot;
  float fthz=rz*fth;
  float fthx=rx*fth;
  mfh=mfh*s;

  for (ix=2;ix<nx-2;ix++) {
    for (iz=2;iz<nz-2;iz++) {
      ioff=iz+ix*nz;
      tgt_p[ioff]=two*src_p[ioff]-tgt_p[ioff] +
	v[ioff]*(motz*(src_p[ioff+2]+   src_p[ioff-2]) +
		 motx*(src_p[ioff+2*nz]+src_p[ioff-2*nz]) +
		 fthz*(src_p[ioff+1]+   src_p[ioff-1]) +
		 fthx*(src_p[ioff+nz]+  src_p[ioff-nz]) +
		 mfh*src_p[ioff]);    
    }
  }
  
  for (ix=0;ix<nx;ix++) {
    tgt_p[ix*nz]=-tgt_p[ix*nz+2];
    tgt_p[ix*nz+nz-1]=-tgt_p[ix*nz+nz-3];
  }
  for (iz=0;iz<nz;iz++) {
    tgt_p[iz]=-tgt_p[iz+2*nz];
    tgt_p[iz+(nx-1)*nz]=-tgt_p[iz+(nx-3)*nz];
  }
  
}


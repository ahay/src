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

#ifdef IWAVE_USE_MPI
#include <mpi.h>
#endif

void step_forward(float * tgt_p, 
		  float * src_p, 
		  float * v,
		  int nz, int nx,
		  float rz, float rx, float s) 
/*< step forward >*/
{

    int iz;
    int ix;
    int ioff; 
    float two=2.0;
    
    for (ix=1;ix<nx-1;ix++) {
	for (iz=1;iz<nz-1;iz++) {
	    ioff=iz+ix*nz;
	    tgt_p[ioff]=two*src_p[ioff]-tgt_p[ioff]
		+v[ioff]*(rz*(src_p[ioff+1]+src_p[ioff-1]) +
			  rx*(src_p[ioff+nz]+src_p[ioff-nz]) - 
			  s*src_p[ioff]);    
	}
    }
}

#ifdef IWAVE_USE_MPI
void mpi_step_forward(float * p0,
		      float * p1,
		      float * v,
		      int nz,
		      int nx,
		      float rz,
		      float rx,
		      float s,
		      int rank,
		      int size,
		      MPI_Status * stat) {

  /*
   * BLOCKING DATA EXCHANGE
   * see notes on implementation dependent deadlock potential -
   * safer to use nonblocking calls here!
   */
  
#ifdef VERB
    fprintf(stderr,"rank=%d send left\n",rank);
#endif
    
    /* send updated p1(second column) to the left */
    if (rank>0) 
	MPI_Send(&(p1[nz]),nz,MPI_FLOAT,rank-1,0,MPI_COMM_WORLD);
    
#ifdef VERB
    fprintf(stderr,"rank=%d recv right\n",rank);
#endif

    /* receive updated p1(last column) from the right */
    if (rank<size-1) 
	MPI_Recv(&(p1[(nx-1)*nz]),nz,MPI_FLOAT,rank+1,0,MPI_COMM_WORLD,stat);
    
#ifdef VERB
    fprintf(stderr,"rank=%d send right\n",rank);
#endif
    
    /* send updated p1(next-to-last column) to the right */
    if (rank<size-1) 
	MPI_Send(&(p1[(nx-2)*nz]),nz,MPI_FLOAT,rank+1,0,MPI_COMM_WORLD);
    
#ifdef VERB
    fprintf(stderr,"rank=%d recv left\n",rank);
#endif
    
    /* receive updated p1(first column) from the left */
    if (rank>0) 
	MPI_Recv(&(p1[0]),nz,MPI_FLOAT,rank-1,0,MPI_COMM_WORLD,stat);
    
    /* construct next time step, overwrite on p0 */
    
    step_forward(p0,p1,v,nz,nx,rz,rx,s);
    
}

#endif



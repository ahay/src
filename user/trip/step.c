#ifdef IWAVE_USE_MPI
#include <mpi.h>
#endif

void step_forward(float * restrict tgt_p, 
		  float * restrict src_p, 
		  float * restrict v,
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
void mpi_step_forward(float * restrict p0,
		      float * restrict p1,
		      float * restrict v,
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



/*This is to test vectrization of  for loop*/

#include <assert.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <sys/times.h>
#include <sys/time.h>

#define         NZ_COEFFS                       14
#define         NY_COEFFS                       16
#define         NX_COEFFS                       16

int main() 
{
  int nz = 352;
  int ny = 704;
  int nx = 704;
  
  float** inArray = new float*[nz];
  float** outArray = new float*[nz];
  
  float outArrayScale;
  
  for (int iz=0; iz<nz; iz++) {
    inArray[iz] = new float[nx*ny+8];
    outArray[iz] = new float[nx*ny+8];
    inArray[iz] = (float*) ( ((unsigned long)inArray[iz])/32*32);
    outArray[iz] = (float*) ( ((unsigned long)outArray[iz])/32*32);
    assert( ((unsigned long)inArray[iz]) %32 == 0 );
    assert( ((unsigned long)outArray[iz]) %32 == 0 );
    
    for (int ixy=0; ixy<nx*ny; ixy++)
      inArray[iz][ixy] = ixy-iz;
    for (int ixy=0; ixy<nx*ny; ixy++)
      outArray[iz][ixy] = ixy-iz;
  }
  
  printf("Allocated %.1f Gbytes \n\n", 1.0e-9*nx*ny*nz*(double)4*2);
  
  //warm up threads.
  for (int iterNum=0; iterNum<3; iterNum++) {
    for (int iz=0; iz<nz; iz++) {
      for (int ixy=0; ixy<nx*ny; ixy++) {
	outArray[iz][ixy] = sqrt(1.11f/inArray[iz][ixy]) ;
      }
    }
  }
  
  // Make an array long enough to test COEFF arrays of length 65
  float coeffs[]=
    { -3.240983757744433369, 1.951616255829886848, -0.4530812397614908393, 0.1776361197426963301, -0.08345825062463696487,
      0.04202539276056183049, -0.0214584977599150406, 0.01070136434836005401, -0.005031015276634687172, 0.002130331875264679375,
      -0.0007451160804013701178, 0.0001565338330123990096, 0.001, -0.001, 0.001, -0.001, 0.001, -0.001, 0.001, -0.001, 0.001,
      -0.001, 0.001, -0.001, 0.001, -0.001, 0.001, -0.001, 0.001, -0.001, 0.001, -0.001, 0.001, -0.001, 0.001, -0.001, 0.001, -0.001,
      0.001, -0.001, 0.001, -0.001, 0.001, -0.001, 0.001, -0.001, 0.001, -0.001, 0.001, -0.001, 0.001, -0.001, 0.001, -0.001, 0.001,
      -0.001, 0.001, -0.001, 0.001, -0.001, 0.001, -0.001, 0.001, -0.001, 0.001, -0.001, 0.001, -0.001 } ;
  
  // do Z derivative without re-ordering memory
  {
    long begin = clock();
    
    struct timeval time1;
    gettimeofday( &time1, NULL);
    
    int nIterations = 4;
    for ( int iterNum=0; iterNum < nIterations; iterNum++ ) {
      
      for (int iy=0; iy<ny; iy++) {
	int ixy= (iy*nx);
	for ( int iz= NZ_COEFFS-1; iz< nz- NZ_COEFFS; iz++) {
	  for ( int ix=0; ix< nx; ix++ ) {
	    //outArray[iz][ixy + ix]= inArray[iz][ixy+ ix]* coeffs[0];
	    outArrayScale = inArray[iz][ixy+ ix]* coeffs[0];
	    for ( int ic=1; ic< NZ_COEFFS; ic++ ) {
	      //outArray[iz][ixy + ix]+= coeffs[ic]* ( inArray[iz+ ic][ixy+ ix] + inArray[iz- ic][ixy+ ix]);
	      // Test for vectorization
	      outArrayScale+= coeffs[ic]* ( inArray[iz+ ic][ixy+ ix] + inArray[iz- ic][ixy+ ix]);
	    }
	    outArray[iz][ixy + ix] = outArrayScale;
	    
	  } // ic
	} //iz
      } // iy
    }// iterNum
    
    long end= clock();
    
    struct timeval time2;
    gettimeofday( &time2, 0 );
    double wallTimeSeconds= time2.tv_sec - time1.tv_sec + time2.tv_usec*1.e-6 -time1.tv_usec*1.e-6 ;
    
    double nFLOPS= nx*ny*(double)nz*NZ_COEFFS*3*nIterations;
    double nReadRegBytes= nx*ny*(double)nz*NZ_COEFFS*2*nIterations*4;
    double nReadL2Bytes= nx*ny*(double)nz*nIterations*4;
    double nWriteBytes= nx*ny*(double)nz*nIterations*4;
    
    printf("\nZ derivative, no mem reorder: \n");
    printf("Performance:  %.1f%% threadding efficiency.  \n",  100.0*((end-begin)/(float)CLOCKS_PER_SEC)/wallTimeSeconds );
    printf("Performance:  %.1f Gflops per core.  \n", nFLOPS*1.e-9/ ((end-begin)/(float)CLOCKS_PER_SEC) );
    printf("Performance:  %.1f Gflops total.  \n", nFLOPS*1.e-9/ wallTimeSeconds );
    printf("              IO for all cores:  %.1f Gbytes/sec read from (& written to) main mem  \n",
	   nReadL2Bytes*1.e-9/ wallTimeSeconds );
    printf("\n");
  } // {
  

  // do x derivative
  {
    long begin= clock();

    struct timeval   time1;
    gettimeofday( &time1, NULL );

    int nIterations= 12;
    for ( int iterNum=0; iterNum < nIterations; iterNum++ ) {

      for ( int iz= 0; iz< nz; iz++) {
	for ( int iy= 0; iy< ny; iy++ )  {
	  float* in2= &inArray[iz][iy*nx+  0 ];
	  float* out2= &outArray[iz][iy*nx+  0 ];
	  // Since we are only operating on a single X vector at a time, memory access should be very fast.      
	  for (int ix= NX_COEFFS; ix< nx- NX_COEFFS; ix++) {
            out2[ ix]= in2[ ix]* coeffs[0];
	    for ( int ic=1; ic< NX_COEFFS; ic++ ) {
	      out2[ ix]+= coeffs[ic]* ( in2[ ix+ ic] + in2[ ix- ic] );
	    }
	  } // ixBlock
	} //iy
      }// iz
    } // iterNum
    
    long end= clock();
    
    struct timeval time2;
    gettimeofday( &time2, 0 );
    double wallTimeSeconds= time2.tv_sec - time1.tv_sec + time2.tv_usec*1.e-6 -time1.tv_usec*1.e-6 ;
    
    double nFLOPS= nx*ny*(double)nz*NX_COEFFS*3*nIterations;
    double nReadRegBytes= nx*ny*(double)nz*NX_COEFFS*2*nIterations*4;
    double nReadL2Bytes= nx*ny*(double)nz*nIterations*4;
    double nWriteBytes= nx*ny*(double)nz*nIterations*4;
    
    printf("X derivative, \n" );
    printf("Performance:  %.1f%% threadding efficiency.  \n",  100.0*((end-begin)/(float)CLOCKS_PER_SEC)/wallTimeSeconds );
    printf("Performance:  %.1f Gflops per core.  \n", nFLOPS*1.e-9/ ((end-begin)/(float)CLOCKS_PER_SEC) );
    printf("Performance:  %.1f Gflops total.  \n", nFLOPS*1.e-9/ wallTimeSeconds );
    printf("              IO for all cores:  %.1f Gbytes/sec read from (& written to) main mem  \n",
	   nReadL2Bytes*1.e-9/ wallTimeSeconds );
    printf("\n");
    
    
  } // {

  // Y derivative 
  {
    long begin= clock();

    struct timeval time1;
    gettimeofday( &time1, NULL );


    int nIterations= 8;
    for ( int iterNum=0; iterNum < nIterations; iterNum++ ) {
      for ( int iz= 0; iz< nz; iz++ ) {
	for ( int iy= NY_COEFFS-1; iy< ny- NY_COEFFS; iy++){
	  int ixy= (iy*nx);
#pragma ivdep
	  for ( int ix=0; ix<nx; ix++ ) {
	    outArray[iz][ixy + ix]= inArray[iz][ixy+ ix]* coeffs[0];
          }

	  for ( int ic=1; ic< NY_COEFFS; ic++ ) {
#pragma ivdep
	   for ( int ix=0; ix<nx; ix++ ) {
	     outArray[iz][ixy + ix]+= coeffs[ic]* ( inArray[iz][ixy+ ic*nx+ ix] + inArray[iz][ixy- ic*nx+ ix]);
            }
	  } // ix
	} // iy
      }  // iz
    } // iterNum
    
    long end= clock();

    struct timeval time2;
    gettimeofday( &time2, 0 );
    double wallTimeSeconds= time2.tv_sec - time1.tv_sec + time2.tv_usec*1.e-6 -time1.tv_usec*1.e-6 ;


    //printf("Final CPU time:  %.2f seconds .   Elapsed time: %.2f seconds\n", (end-begin)/(float)CLOCKS_PER_SEC, wallTimeSeconds );                                                               


    double nFLOPS= nx*nz*(double)ny*NY_COEFFS*3*nIterations;
    double nReadRegBytes= nx*nz*(double)ny*NY_COEFFS*2*nIterations*4;
    double nReadL2Bytes= nx*nz*(double)ny*nIterations*4;
    double nWriteBytes= nx*nz*(double)ny*nIterations*4;

    printf("Y derivative, no mem reorder: \n");
    printf("Performance:  %.1f%% threadding efficiency.  \n",  100.0*((end-begin)/(float)CLOCKS_PER_SEC)/wallTimeSeconds );
    printf("Performance:  %.1f Gflops per core.  \n", nFLOPS*1.e-9/ ((end-begin)/(float)CLOCKS_PER_SEC) );
    printf("Performance:  %.1f Gflops total.  \n", nFLOPS*1.e-9/ wallTimeSeconds );
    printf("              IO for all cores:  %.1f Gbytes/sec read from (& written to) main mem  \n",
	   nReadL2Bytes*1.e-9/ wallTimeSeconds );
    printf("\n");


  } //{
    
}
  

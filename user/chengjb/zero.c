#include <rsf.h>
#include "_cjb.h"

void zero1int(int *p, int n1)
/*< free a 1-d array of int >*/
{
     int i;
	 #ifdef OPENMP
	 #pragma omp parallel for private(i) 
	 #endif
     for(i=0;i<n1;i++) p[i]=0;
}

void zero2int(int **p, int n1, int n2)
/*< free a 2-d array of int >*/
{
     int i, j;
	 #ifdef OPENMP
	 #pragma omp parallel for private(i,j) 
	 #endif
     for(i=0;i<n2;i++) 
       for(j=0;j<n1;j++)
         p[i][j]=0;
}

void zero3int(int ***p, int n1, int n2, int n3)
/*< free a 3-d array of int >*/
{
     int i, j, k;
	 #ifdef OPENMP
	 #pragma omp parallel for private(i,j,k) 
	 #endif
     for(i=0;i<n3;i++) 
       for(j=0;j<n2;j++)
          for(k=0;k<n1;k++)
            p[i][j][k]=0;
}

void zero3ushort(unsigned short ***p, int n1, int n2, int n3)
/*< free a 3-d array of unsigned short >*/
{
     int i, j, k;
	 #ifdef OPENMP
	 #pragma omp parallel for private(i,j,k) 
	 #endif
     for(i=0;i<n3;i++) 
       for(j=0;j<n2;j++)
          for(k=0;k<n1;k++)
            p[i][j][k]=0;
}

void zero2ushort(unsigned short **p, int n1, int n2)
/*< free a 2-d array of unsigned short >*/
{
     int i, j;
	 #ifdef OPENMP
	 #pragma omp parallel for private(i,j) 
	 #endif
     for(i=0;i<n2;i++) 
       for(j=0;j<n1;j++)
            p[i][j]=0;
}

void zero1float(float *p, int n1)
/*< free a 1-d array of float >*/
{
     int i;
	 #ifdef OPENMP
	 #pragma omp parallel for private(i) 
	 #endif
     for(i=0;i<n1;i++) p[i]=0.0;
}

void zero2float(float **p, int n1, int n2)
/*< free a 2-d array of float >*/
{
     int i, j;
	 #ifdef OPENMP
	 #pragma omp parallel for private(i,j) 
	 #endif
     for(i=0;i<n2;i++) 
       for(j=0;j<n1;j++)
         p[i][j]=0.0;
}

void zero3float(float ***p, int n1, int n2, int n3)
/*< free a 3-d array of float >*/
{
     int i, j, k;
	 #ifdef OPENMP
	 #pragma omp parallel for private(i,j,k) 
	 #endif
     for(i=0;i<n3;i++) 
       for(j=0;j<n2;j++)
          for(k=0;k<n1;k++)
            p[i][j][k]=0.0;
}

void zero4float(float ****p, int n1, int n2, int n3, int n4)
/*< free a 4-d array of float >*/
{
     int m, i, j, k;
	 #ifdef OPENMP
	 #pragma omp parallel for private(m,i,j,k) 
	 #endif
     for(m=0;m<n4;m++) 
       for(i=0;i<n3;i++) 
          for(j=0;j<n2;j++)
            for(k=0;k<n1;k++)
               p[m][i][j][k]=0.0;
}

void zero1double(double *p, int n1)
/*< free a 1-d array of double >*/
{
     int i;
	 #ifdef OPENMP
	 #pragma omp parallel for private(i) 
	 #endif
     for(i=0;i<n1;i++) p[i]=0.0;
}

void zero2double(double **p, int n1, int n2)
/*< free a 2-d array of double >*/
{
     int i, j;
	 #ifdef OPENMP
	 #pragma omp parallel for private(i,j) 
	 #endif
     for(i=0;i<n2;i++) 
       for(j=0;j<n1;j++)
         p[i][j]=0.0;
}

void zero3double(double ***p, int n1, int n2, int n3)
/*< free a 3-d array of double >*/
{
     int i, j, k;
	 #ifdef OPENMP
	 #pragma omp parallel for private(i,j,k) 
	 #endif
     for(i=0;i<n3;i++) 
       for(j=0;j<n2;j++)
          for(k=0;k<n1;k++)
            p[i][j][k]=0.0;
}

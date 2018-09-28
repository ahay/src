#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif

int main() 
{

    int niter=100000, sum=0, i;

#pragma omp parallel for shared(sum)
    for (i=0; i<=niter; i++) {
       sum += i; 
    }

    sf_warning("sum=%d",sum);

    exit(0);
}

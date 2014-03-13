#include "rsf.hh"
#include "vai.hh"
using namespace std;

// constructor
VAI::VAI(int m0, int m1)
{
    n[0]=m0;
    n[1]=m1;
}

// () operator 
int VAI::operator() (int i0,int i1)
{
    int ii;
    ii = i1*n[0] + i0;
    return ii;
}

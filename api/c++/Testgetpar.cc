#include <assert.h>

#include "rsf.hh"

int main (void)
{
    int i, argc=5;
    float f, d[4];
    bool yes, no[2];
    const char *argv[] = {"prog","a=5","b=as","a=100","par=Testgetpar.cc"};
//    char *str;
    iRSF par(0);

//      a=100 Xa=5
//      float=5.625 cc=fgsg
//      dd=1,2x4.0,2.25 true=yes false=2*no label="Time (sec)"

    sf_init(argc,(char**) argv);
    par.get("a",i);
    assert (i==100);
    par.get("c",i,0);
    assert (i==0);
    par.get("float",f);
    assert (f==5.625);
    par.get("dd",4,d);
    assert (d[0]==1 && d[1]==4 && d[2]==4 && d[3]==2.25);
    par.get("true",yes);
    assert (yes);
    par.get("false",2,no);
    assert (!no[0] && !no[1]);
//   assert (NULL != (str = sf_getstring("label")) && 
//	    0==strcmp(str,"Time (sec)"));

    exit (0);
}

// 	$Id: Testgetpar.cc 982 2005-01-30 23:38:22Z shan $	

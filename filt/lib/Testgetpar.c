#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "getpar.h"
#include "_bool.h"

int main (void)
{
    int i, argc=5;
    float f, d[4];
    bool yes, no[2];
    char *argv[] = {"prog","a=5","b=as","a=100","par=Testgetpar.c"};
    char *str;

    /*
      a=100 Xa=5
      float=5.625 cc=fgsg
      dd=1,2x4.0,2.25 true=yes false=2*no label="Time"
    */

    sf_init(argc,argv);
    assert (sf_getint("a",&i) && i==100);
    assert (!sf_getint("c",&i));
    assert (sf_getfloat("float",&f) && f==5.625);
    assert (sf_getfloats("dd",d,4) && 
	    d[0]==1 && d[1]==4 && d[2]==4 && d[3]==2.25);
    assert (sf_getbool("true",&yes) && yes);
    assert (sf_getbools("false",no,2) && !no[0] && !no[1]);
    assert (NULL != (str = sf_getstring("label")) && 
	    0==strcmp(str,"Time"));

    exit (0);
}

/* 	$Id$	 */


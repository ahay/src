/*  Make a synthtic two-layer CMP gather with known t0
*/
/*
  Copyright (C) 2014 Jilin University
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include <rsf.h>

int main(int argc, char* argv[])
{ 
    int i,j,n1,n2,evt1,evt2;
    float o1,o2,dt,dx,t01,t02,v01,v02,f1,f2;
    float *trace;
    bool verb;
    sf_file in, spcmp;
    
    sf_init (argc,argv);
    
    if (!sf_stdin()) { /* no input file in stdin */
	in = NULL;
    } else {
	in = sf_input("in");
    }

    spcmp = sf_output("out");
    
    if (NULL == in) {
	sf_setformat(spcmp,"native_float");
    } else if (SF_FLOAT != sf_gettype(in)) {
	sf_error("Need float input");
    }
    /*get dimension parameter*/
    if(!sf_getint("n1",&n1)) n1=1000;
    /* number of n1*/
    if(!sf_getint("n2",&n2)) n2=20;
    /* number of n2*/
    if(!sf_getfloat("dt",&dt)) dt=0.001;
    /* sampling on 1-th axis(time)*/
    if(!sf_getfloat("dx",&dx)) dx=50;
    /* sampling on 2-th axis(offset)*/
    if(!sf_getfloat("o1",&o1)) o1=0;
    if(!sf_getfloat("o2",&o2)) o2=0;
    if(!sf_getfloat("v01",&v01)) v01=1000;
    /* first event rms vel */
    if(!sf_getfloat("v02",&v02)) v02=1000;
    /* second event rms vel */
    if(!sf_getfloat("t01",&t01)) t01=0.4;
    /* t01 start point */
    if(!sf_getfloat("t02",&t02)) t02=0.8;
    /* t02 start point */
    if (!sf_getbool("verb",&verb)) verb = false;
    
    /* dimensions */
    sf_putint(spcmp,"n1",n1);
    sf_putint(spcmp,"n2",n2);
    sf_putfloat(spcmp, "d1", dt);
    sf_putfloat(spcmp, "d2", dx);
    sf_putfloat(spcmp, "o1", o1);
    sf_putfloat(spcmp, "o2", o2);
    
    trace = sf_floatalloc(n1*n2);
    
    for(i = 0; i < n1*n2 ;i++){
	trace[i]=0;
    }
    
    for(j = 0; j < n2; j++){
	for(i = 0; i < n1; i++){
	    f1 = t01*t01+powf(j*dx,2)/powf(v01,2);
	    f2 = t02*t02+powf(j*dx,2)/powf(v02,2);
	    evt1 = sqrtf(f1)/dt;
	    evt2 = sqrtf(f2)/dt;
	    
	    if(i==evt1){
		trace[i+j*n1]=1;
		if(verb) sf_warning("t01= %d;\n",evt1);
		    };
	    if(i==evt2){
		trace[i+j*n1]=1;
		if(verb) sf_warning("t02= %d;\n",evt2);
	    }
	}
    }
    
    sf_floatwrite(trace,n1*n2,spcmp);
    
    exit(0);
}

/* 	$Id$	 */

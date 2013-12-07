/* Mask design for dip angle gathers. */
/*
  Copyright (C) 2013 University of Texas at Austin
  
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
    int time_sample_number, *mleft, *mright, maxangle,*mlength;
    float mask_width, **mask=NULL;
    sf_file out=NULL; /* output file */


	max_dip_angle=par.float("maxangle",20) # maximum dip angle

	time_sample_number=100
	mask_width=3./4
   	if (!sf_getfloat("mask_width",&mask_width)) sf_error("Need mask width=");
    	if (!sf_getfloat("time_sample_number",&time_sample_number)) sf_error("Need time sample number=");
    	if (!sf_getfloat("max_dip_angle",&max_dip_angle)) sf_error("Need max dip angle=");
    	if (!sf_getfloat("mask_width",&mask_width)) sf_error("Need mask width=");

	for (i1=0,i1<time_sample_number,i1++){
	mlength[i1]=i1;
	}

	a=time_sample_number;
	b=floor(mask_width*max_dip_angle);

	width=2*max_dip_angle+1;
	
	for (i1=0,i1<time_sample_number,i1++){
	mleft[i1]=floor(b*(1-(mlength[i1])**2/a**2)**.5);
	}


	for (i1=0,i1<time_sample_number,i1++){
	mlength[i1]=i1;
	}


	mleft=floor(b*(1-(mlength)**2/a**2)**.5)
	mleft=mleft.astype(int)
	mright=-mleft+width

	for (i1=0,i1<time_sample_number+1,i1++){
		for (j1=0,j1<width+1, j1++){
		if ()

		}

	}
	


    /* Initialize RSF */
    sf_init(argc,argv);
    /* standard input */
    in = sf_input("in");
    /* standard output */
    out = sf_output("out");

    /* check that the input is float */
    if (SF_FLOAT != sf_gettype(in))
	sf_error("Need float input");

    /* parameter from the command line (i.e. clip=1.5 ) */
    if (!sf_getfloat("clip",&clip)) sf_error("Need clip=");

    /* allocate floating point array */
    trace = sf_floatalloc (n1);

    /* loop over traces */
    for (i2=0; i2 < n2; i2++) {

	/*read a trace */
	sf_floatread(trace,n1,in);

	/* loop over samples */
	for (i1=0; i1 < n1; i1++) {
	    if      (trace[i1] >  clip) trace[i1]= clip;
	    else if (trace[i1] < -clip) trace[i1]=-clip;
	}

	/* write a trace */
	sf_floatwrite(trace,n1,out);
    }

    sf_close();
    exit(0);
}

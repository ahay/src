/* Mask for Dip Angle Gathers */
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

#include <stdio.h>
#include <math.h>
#include <rsf.h>

int main(int argc, char* argv[])
{

    int n1,n2, i1, j1, time_sample_number, maxangle,a,b;
    float mask_width;




    sf_file  out=NULL; /* Input and output files */

    /* Initialize RSF */
    sf_init(argc,argv);

    /* standard output */
    out = sf_output("out");

    /* check that the input is float */
//    if (SF_FLOAT != sf_gettype(in))
//	sf_error("Need float input");

    /* parameter from the command line (i.e. clip=1.5 ) */
    if (!sf_getfloat("mask_width",&mask_width)) sf_error("Need mask width=");

    if (!sf_getint("time_sample_number",&time_sample_number)) sf_error("Need time sample number=");
    if (!sf_getint("maxangle",&maxangle)) sf_error("Need maxangle=");




    int mlength[time_sample_number];
    int mleft[time_sample_number];
    int mright[time_sample_number];
    a=time_sample_number;
    b=round(mask_width*maxangle);
    int width=2*maxangle+1;
    sf_putint (out, "n1", time_sample_number);
    sf_putint (out, "n2", width);    
    float** mask = sf_floatalloc2 (n1, n2);


    for (i1=0;i1<time_sample_number;i1++){
    mlength[i1]=i1;
    }

    for (i1=time_sample_number-1;i1>=0;i1--){
    mleft[i1]= round(b*sqrt((1-pow(mlength[i1],2)/pow(a,2))));
    mright[i1]= width-mleft[i1];
    }

    for (i1=0;i1<time_sample_number;i1++){
        for (j1=0;j1<width; j1++){
        if (j1>=mleft[i1] && j1<= mright[i1]){
            mask[i1][j1]=1.;
        }
        else{
            mask[i1][j1]=0.;
        }

        }

    }
    sf_floatwrite(mask,n1,out);
    sf_close();
    return 0;
}

/* Hungs a 1d velocity function from the Water bottom.
   Should work for 2D models



   stdin    1D velocity function to be used 
   file mask [required]   The water bottom is read from the mask file.
                          1 above the WB
                          0 bellow the WB

   stdout The output velocity model has dimensions of the mask file.
   vel [1.5]    velocity to use above the horizon (usually water velocity) 



   

*/
/*
  Copyright (C) 2011 Colorado School of Mines

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

int main (int argc, char* argv[])
{

    int i1,i2,i3,j1,j2;

    int n1,n2,n3;
    float o1,o2,o3;
    float d1,d2,d3,x;
    
    sf_axis ax1,m1,m2,m3;

    int n1d;     // samples of 1d function
    float o1d;   // origin of 1d function (has to be zero)
    float d1d;   // 1D function smapling (has to be the same as d1)

    float vel;
    float *v1d,**mask,**vmod,*wb;

    sf_file in=NULL, out=NULL, mask1=NULL, wbot=NULL;



    sf_init (argc,argv);
    in= sf_input("in");
    out = sf_output("out");

    //=====================================================
    //Get parameters from stdin file

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    
    
   /* parameters from input file*/


    m1=sf_iaxa(in,1);
    m2=sf_iaxa(in,2);
    m3=sf_iaxa(in,3);

    n1=sf_n(m1);  o1=sf_o(m1); d1=sf_d(m1); 
    n2=sf_n(m2);  o2=sf_o(m2); d2=sf_d(m2); 

	sf_oaxa(out,m1,1);
	sf_oaxa(out,m2,2);

    




//  =======================================
    //allocate input mask file
    mask=sf_floatalloc2(n1,n2);

    sf_floatread(mask[0],n1*n2,in);

//  =======================================
    printf("hola %d%n",1);
    for (i2=0; i2<n2; i2++) {
        for (i1=1; i1<n1; i1++){
           if(mask[i2][i1]-mask[i2][i1-1] != 0.0) wb[i2]=i1*d1+o1;
        }
    }

    x=0.0;
    for (i2=0; i2<n2 ; i2++){
        x=x+d2;
        printf("%15.6f%15.2f%n",x,wb[i2]);
    }




    sf_floatwrite(mask[0],n1*n2,out);

    exit(0);
}

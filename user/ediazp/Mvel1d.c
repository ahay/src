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

    int i1,i2,j1,j2;

    int n1,n2;
    float o1;
    float d1;
    
    sf_axis ax1,m1,m2;

    int n1d;     // samples of 1d function
    float d1d;   // 1D function smapling (has to be the same as d1)

    float vel;
    float *v1d,**mask,**vmod,*wb;

    sf_file in=NULL, out=NULL, mask1=NULL, wbot=NULL;



    sf_init (argc,argv);
    in = sf_input("in");
    mask1= sf_input("mask");
    out = sf_output("out");
    wbot=sf_output("wb");

    //=====================================================    
    //Get parameters from command line:
    
    if (! sf_getfloat("vel", &vel)) vel=1.5;
    
    //=====================================================
    //Get parameters from stdin file

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    
    
   /* parameters from input file*/
    ax1=sf_iaxa(in,1);
    n1d=sf_n(ax1); d1d=sf_d(ax1);


    m1=sf_iaxa(mask1,1);
    m2=sf_iaxa(mask1,2);

    n1=sf_n(m1);  o1=sf_o(m1); d1=sf_d(m1); 
    n2=sf_n(m2); 

	sf_oaxa(out,m1,1);
	sf_oaxa(out,m2,2);
	sf_oaxa(wbot,m2,1);

    
    if( d1!=d1d)  sf_error("1d vel has different sampling par than mask");




//  =======================================
    //allocate input mask file
    mask=sf_floatalloc2(n1,n2);

    //allocate output model file
    vmod=sf_floatalloc2(n1,n2);

    // water bottom vector
    wb= sf_floatalloc(n2);

    // stdin 1d function
    v1d= sf_floatalloc(n1);

    //read 1d function
    sf_floatread(v1d,n1d,in);

    //read mask
    sf_floatread(mask[0],n1*n2,mask1);

//  =======================================

    for (i2=0; i2<n2; i2++) {
        for (i1=1; i1<n1; i1++){
           if(mask[i2][i1]-mask[i2][i1-1] != 0.0) wb[i2]=i1*d1+o1;
        }
    }




    for (i2=0; i2<n2; i2++) {
        j1=(int) ((wb[i2]-o1)/d1);
        for (i1=0; i1<j1; i1++){
             vmod[i2][i1]=vel;
        }
        for (i1=j1; i1<n1 ;i1++){
            j2=i1-j1;
            if(j2>n1d-1){ 
                vmod[i2][i1]=v1d[n1d-1];
            }else {
                vmod[i2][i1]=v1d[j2];
            }
            
        }
    }

    sf_floatwrite(vmod[0],n1*n2,out);
    sf_floatwrite(wb,n2,wbot);
    free (*mask); free (mask); 
    free (*vmod); free (vmod); 
    free (wb); 
    free (v1d);

    exit(0);
}

/* 
Applies a 1,2, or 3D convolution kernel or its adjoint
The filter is composed by n coefficients.
 
example: 2d laplacian

     1
  1 -4  1
     1 

filter: 1 1 -4 1 1
lag1 (vertical lag)  :  0  1  0  -1  0 
lag2 (horizontal lag): -1  0  0   0  1  

*/
/*
  Copyright (C) 2007 Colorado School of Mines
  
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
#include "convkernel.h"
void set_axis(sf_file file, sf_axis ax, int axis);
int main(int argc, char* argv[])
{
  bool adj; 
  float *input1, *output1, *filter;
  float **input2, **output2;
  float ***input3, ***output3;
  int *lag1,*lag2,*lag3,nd;
  int n1,n2,n3,nl,ndim,  n[SF_MAX_DIM];
  int ntraces, itr;
  sf_axis axlag; 

  sf_file Fin=NULL; /* velocity  */
  sf_file Fout=NULL; /* density */
  sf_file Flag1=NULL; /* filter lags */
  sf_file Flag2=NULL; /* filter lags */
  sf_file Flag3=NULL; /* filter lags */
  sf_file Ffilter=NULL; /* filter lags */

  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  /* init RSF */
  sf_init(argc,argv);

  if(! sf_getbool("adj",&adj)) adj=false;
  if(! sf_getint("n",&nd)) nd=1;
  
  /*------------------------------------------------------------*/

  /*------------------------------------------------------------*/
  /* I/O files */
  
  Fin = sf_input ("in" ); /* wavelet   */
  Ffilter = sf_input ("filter"); /* filter lags */
  Fout = sf_output ("out"); /* velocity  */
  
  ndim =sf_filedims (Fin,n);

  if(ndim==1 || nd==1){
    sf_warning("1d conv,%d",ndim);
    n1 = n[0];
    input1 = sf_floatalloc(n1);
    output1 = sf_floatalloc(n1);
    Flag1 = sf_input ("lag"); /* filter lags */

    axlag = sf_iaxa(Flag1,1); 
    nl = sf_n(axlag);
    filter = sf_floatalloc(nl);
    lag1 = sf_intalloc(nl);
    sf_floatread(filter,nl,Ffilter);
    sf_intread(lag1,nl,Flag1);

    ntraces = sf_leftsize(Fin,1);

    sf_warning("read n1=%d tr=%d",n1,ntraces);
    convkernel1_init(n1,nl,lag1,filter);
    for (itr=0; itr<ntraces; ++itr){
      sf_floatread(input1,n1,Fin);
      if (adj){
        convkernel1_apply(output1,input1,adj);
      }else{
        convkernel1_apply(input1,output1,adj);
      }
      sf_floatwrite(output1,n1,Fout);
    }


    free(input1);
    free(output1);
    free(lag1);
    free(filter);

  }else if(ndim==2){ 
    sf_warning("2d conv");
    n1 = n[0];
    n2 = n[1];
    input2 = sf_floatalloc2(n1,n2);
    output2 = sf_floatalloc2(n1,n2);
    sf_floatread(input2[0],n1*n2,Fin);
    Flag1 = sf_input ("lag1"); /* filter lags */
    Flag2 = sf_input ("lag2"); /* filter lags */

    axlag = sf_iaxa(Flag1,1); 
    nl = sf_n(axlag);
    filter = sf_floatalloc(nl);
    lag1 = sf_intalloc(nl);
    lag2 = sf_intalloc(nl);
    sf_floatread(filter,nl,Ffilter);
    sf_intread(lag1,nl,Flag1);
    sf_intread(lag2,nl,Flag2);
    convkernel2_init(n1,n2,nl,lag1,lag2,filter);
    if (adj){
      convkernel2_apply(output2,input2,adj);
    }else{
      convkernel2_apply(input2,output2,adj);
    }
    sf_floatwrite(output2[0],n1*n2,Fout);
    free(*input2); free(input2);
    free(*output2); free(output2);
    free(lag1);
    free(lag2);
    free(filter);

  }else if(ndim==3){
    sf_warning("3d conv");
    n1 = n[0];
    n2 = n[1];
    n3 = n[2];

    input3 = sf_floatalloc3(n1,n2,n3);
    output3 = sf_floatalloc3(n1,n2,n3);
    sf_floatread(input3[0][0],n1*n2*n3,Fin);
    Flag1 = sf_input ("lag1"); /* filter lags */
    Flag2 = sf_input ("lag2"); /* filter lags */
    Flag3 = sf_input ("lag3"); /* filter lags */


    axlag = sf_iaxa(Flag1,1); 
    nl = sf_n(axlag);
    filter = sf_floatalloc(nl);
    lag1 = sf_intalloc(nl);
    lag2 = sf_intalloc(nl);
    lag3 = sf_intalloc(nl);
    sf_floatread(filter,nl,Ffilter);
    sf_intread(lag1,nl,Flag1);
    sf_intread(lag2,nl,Flag2);
    sf_intread(lag3,nl,Flag3);
    convkernel3_init(n1,n2,n3,nl,lag1,lag2,lag3,filter);
    if (adj){
      convkernel3_apply(output3,input3,adj);
    }else{
      convkernel3_apply(input3,output3,adj);
    }
    sf_floatwrite(output3[0][0],n1*n2*n3,Fout);
    free(**input3); free(*input3); free(input3);
    free(**input3); free(*output3); free(output3);
    free(lag1);
    free(lag2);
    free(lag3);
    free(filter);  
  }else{
      sf_warning("X conv");
      exit(1);
  }


  exit (0);
}


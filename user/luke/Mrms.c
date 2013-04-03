/* Local RMS Determination for an array of arbitrary dimension.  The absolute value of a unit is used for indicies falling  between within the edge and the sampling length for each dimension.
*/
/*
  Copyright (C) 2004 University of Texas at Austin

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
    int dim, dim1, i, k, l, m, n[SF_MAX_DIM], rect[SF_MAX_DIM], s[SF_MAX_DIM], t[SF_MAX_DIM], test;
    int n1, n2, n3, indexalloc, call, indexalloc1, dataindexpoint;
    char key[6];
    float* data;
    float* attr;
    int* dataindex;
    int* dataindex1;
    float anal, holder, mean, rootmean, coordnum;
    sf_file in, out;

    sf_init (argc, argv);
    in  = sf_input ("in");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    dim = sf_filedims (in,n);
    dim1 = -1;
    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
	/*( rect#=(1,1,...) number of samples to use on #-th axis )*/ 
	if (rect[i] > 0) dim1 = i;

    }//end i loop
    n1 = n2 = 1;
    for (i=0; i < dim; i++) {
        t[i] = 0;
	if (i <= dim1) {
	    s[i] = n1;
	    n1 *= n[i];
            t[i] = n1;
	} else {
	    n2 *= n[i];
	}

    }
//                      sf_warning("dim1=%i",dim1);
    n3=n1*n2;
    

   //define index for shift for each dimension
/*   index[1]=1
   for (i=1; i< dim; i++){
   index[i]=index[i-1]*n[i-1];
   }*/




//Determine Array of values for calling desired cells to perform attribue analysis on.
                  indexalloc = 1;
                  dataindex = sf_intalloc (indexalloc);
                  dataindex[0]=0;
                  indexalloc1 = 1;
                  
          for (k=0; k <= dim1; k++){//loop through rectified dimensions
   
                  indexalloc1 *= 2*rect[k];
 //                     sf_warning("rect[k]=%i",rect[k]);
                   dataindex1 = sf_intalloc (indexalloc1);
              for (m = 0; m<indexalloc; m++){
                  for (l=0; l < 2*rect[k]; l++){
                         dataindexpoint = 2*m*rect[k]+l;
                         if (k > 0){
                         dataindex1[dataindexpoint] = dataindex[m] + (l-rect[k])*t[k-1];
                         }else{
                         dataindex1[dataindexpoint] = dataindex[m]+ l - rect[k];
                         }//end dataindex conditional
                       } //end l loop
                  }//end m loop
                  indexalloc=indexalloc1;


                  free (dataindex);
                  dataindex = sf_intalloc (indexalloc);   

                  //write values to other index
                       for (l=0; l<indexalloc; l++){
                           dataindex[l] = dataindex1[l];

                           } // end l loop
 
//                      sf_warning("indexalloc=%i",indexalloc);

                  free (dataindex1);

              }//end k loop
 
//Read Data
data = sf_floatalloc (n3);
sf_floatread(data,n3,in);
//attr = sf_floatalloc (n3);    
attr = sf_floatalloc (1);  

//figure out which values you want to call to use attribute analysis on.
//let's do this the easy way:  boundary conflict regions are unaltered


int coord[dim]; 
//loop through all data points
for (i=0; i<n3; i++){
//first, determine the coordinate of each called point for boundary conflict test

          //determine coordinate
          coordnum = i;
//                      sf_warning("initial index!%g",coordnum);
          test = 0;
          for (m=0; m<dim; m++){
               k = dim - m - 1; // make it a falling value loop
               if (k > 0){
                     coord[k] = floor(coordnum/t[k-1]);
                     coordnum = coordnum - floor(coordnum/t[k-1])*t[k-1];
               }else{
                     coord[k] = coordnum;
               }//end declare coordinate conditional
//                      sf_warning("i=%i",i);
//                      sf_warning("k!%i",k);
//                       sf_warning("t[k-1]!%i",t[k-1]);
//                       sf_warning("Coord!%i",coord[k]);
//                       sf_warning("Coordnum!%g",coordnum);


               //test boundary conditions
               if ((coord[k] - rect[k] >=0) && (coord[k] + rect[k] <= n[k])){
                  test += 1;
                  
               }else{
                    break;
               }//end BC test
             }//end coord k loop
          attr[0] = 0;
          holder = 0;
//          sf_warning("tst!%i",test);
  //                     sf_warning("dim!%i",dim);
 //         sf_warning("dim %i",dim);
 //         sf_warning("test %i",test);
          if (test == dim){//go ahead, BC ARE A-O-K!
              //loop through dataindex values to call data for the attribute analysis.
//             sf_warning("in loop"); 
              attr[0] = 0; //zero out attribute
              m = 0; //zero out mean number

//sf_warning("indexalloc1 %i",indexalloc1);
              for (k=0; k < indexalloc1; k++){
                   call = i + dataindex[k];
//                   sf_warning("call %i",call);
                   anal = data[call];//value used for attribute analysis.

                   holder += anal*anal; //add the square of the analysis number
                   m += 1.; //add index for mean, is a float for division purposes
//                           sf_warning("m %i",m);             
//                           sf_warning("holder %g",holder); 
                   }//end k loop
              if (m > 0){
                 mean = holder/m;

//                sf_warning("m %i",m);
                 }else{
                 mean = data[i]*data[i];
                 }

              rootmean = sqrtf(mean);
              //and finally the result
              attr[0] = rootmean;
//             sf_warning("rootmean %g",attr[0]);
             }//end determine attribute if
          else{
             attr[0]=sqrtf(data[i]*data[i]);//call edge effect regions their initial value.
              }//end attribute else for boundary conditions.
   sf_floatwrite(attr,1,out);//write data 
}//end i loop

    // write data
//   sf_floatwrite(attr,n3,out); //no longer writing completed data
                                 //writing partial data saves memory.

    sf_fileclose (in);
    sf_fileclose (out);

    free (attr);
    free (data);
}//end


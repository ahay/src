/* Local absolute value determination for an array of arbitrary dimension.  The absolute value of a unit is used for indicies falling  between within the edge and the sampling length for each dimension.
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
    int dim, dim1, i, i2, k, l, m, n[SF_MAX_DIM], rect[SF_MAX_DIM], s[SF_MAX_DIM], test;
    int n1, n2, indexalloc, call, indexalloc1, dataindexpoint;
    char key[6];
    float* data;
    float* attr;
    int* dataindex;
    int* dataindex1;
    float anal, holder, mean, coordnum;
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
	if (rect[i] > 1) dim1 = i;

    }//end i loop


//this loop determines the number of points to read and write.  n1 is the number of points in rectified dimensions.  n2 is the number of unrectified dimensional points, the number of loops around which you must run n1 iterations to complete the analysis. n1*n2=total number of poitns in the read data.
    n1 = n2 = 1;
    for (i=0; i < dim; i++) {
	if (i <= dim1) {
	    s[i] = n1;
	    n1 *= n[i];
	} else {
	    n2 *= n[i];
	}

    }
    
//Determine Array of values for calling desired cells to perform attribue analysis on.
                  indexalloc = 1;
                  dataindex = sf_intalloc (indexalloc);
                  dataindex[0]=0;
                  indexalloc1 = 1;
                  
          for (k=0; k<dim1; k++){//loop through rectified dimensions
   
                  indexalloc1 *= 2*rect[k];

                   dataindex1 = sf_intalloc (indexalloc1);
              for (m = 0; m<indexalloc; m++){
                  for (l=0; l < 2*rect[k]; l++){
                         dataindexpoint = 2*m*rect[k]+l;
                         if (k > 0){
                         dataindex1[dataindexpoint] = dataindex[m] + (l-rect[k])*s[k-1];
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
                  free (dataindex1);
              }//end k loop
//Read Data

//allocate data matrix.  This has value of the size of all datapoints in rectified dimensions.  Unrectified dimensions are read in loop.
data = sf_floatalloc (n1);
attr = sf_floatalloc (1);  //only need to read & write item at a time

//figure out which values you want to call to use attribute analysis on.
//let's do this the easy way:  boundary conflict regions are unaltered

int coord[dim1]; //only need to know coordinates in rectified dimensions

for (i2=0; i2 < n2; i2++) { //loop through unrectified dimensions for greater memory efficiency
	sf_floatread(data,n1,in);  //read data each time looping through an unrectified dimension


//loop through rectified data points.  Remember n1*n2=total points
for (i=0; i<n1; i++){
//first, determine the coordinate of each called point for boundary conflict test

          //determine coordinate
          coordnum = i;
//                      sf_warning("initial index!%g",coordnum);
          test = 0;
          for (m=0; m<dim1; m++){
               k = dim1 - m - 1; // make it a falling value loop
               if (k > 0){
                     coord[k] = floor(coordnum/n[k-1]);
                     coordnum = coordnum - floor(coordnum/n[k-1])*n[k-1];
               }else{
                     coord[k] = coordnum;
               }//end declare coordinate conditional

               //test boundary conditions
               if ((coord[k] - rect[k] >=0) && (coord[k] + rect[k] <= n[k])){
                  test += 1;
                  
               }else{
                    break;
               }//end BC test
             }//end coord k loop
          attr[0] = 0;
          holder = 0;

          if (test == dim1){//go ahead, BC ARE A-O-K!
              //loop through dataindex values to call data for the attribute analysis.

              attr[0] = 0; //zero out attribute
              m = 0; //zero out mean number


//TIME TO COOK UP THE ACTUAL ATTRIBTUE ANALYSIS
//                           ))
//                          ((
//                    ___o___)
//      ___           |     |====O
//     (0 0)          |_____|
//--ooO-(_)-Ooo--------------------

              for (k=0; k < indexalloc1; k++){
                   call = i + dataindex[k];

                   anal = data[call];//value used for attribute analysis.

                   holder += fabsf(anal); //add the amplitude of the analysis number
                   m += 1.; //add index for mean, is a float for division purposes

                   }//end k loop
              if (m > 0){
                 mean = holder/m;


                 }else{
                 mean = fabsf(data[i]);
                 }
              //and finally the result
              attr[0] = mean;
             }//end determine attribute if
          else{
             attr[0]=fabsf(data[i]);//call edge effect regions their initial value.
              }//end attribute else for boundary conditions.

//                   - -
//                 { 0 0 }
// +===========oOO===(_)===OOo=========+
// |_____|_____|_____|_____|_____|_____|
// |__|_____|_____|_____|_____|_____|__|
// |_____| And this is the end   |_____|
// |__|__|      of actual        |__|__|
// |_____| Attribute Analysis    |_____|
// |__|__|_______________________|__|__|
// |_____|_____|_____|_____|_____|_____|
// +===================================+


   sf_floatwrite(attr,1,out);//write data 
}//end i loop

}//end i2 loop
    //free up used files and arrays
    sf_fileclose (in);
    sf_fileclose (out);

    free (attr);
    free (data);
}//end


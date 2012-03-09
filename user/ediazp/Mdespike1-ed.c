/* Despike filter:
   move outliers values to the tolerance
   parameter. The mean is calculated with moving
   windows 

   Example:
   
  
   if (a>3sigma) a=3sigma
 
                        outlier
          * *              ^
         *    *            |
        *      *           |
       *        *          *
   ****          * * * * *  *****
 
 
   
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
    int n1;
    int i3,i2,i1,k;
    int arr_size, t2loop;
    int m;
    int klo, khi;
    float sigma;
    float *uo,*der,*der_tmp;

    float sum; 
    
    sf_file in=NULL, out=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    
    
    //=====================================================    
    //Get parameters from command line:
    
    if (! sf_getint("window",&m)) m=20;
    if (! sf_getfloat("sigma", &sigma)) sigma=3.0;
    
    
    
    //=====================================================
    //Get parameters from input file

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    arr_size=n1;
    t2loop=sf_leftsize(in,1) ;

    uo  = sf_floatalloc(arr_size);
    der = sf_floatalloc(arr_size);
    der_tmp = sf_floatalloc(arr_size);
    fprintf(stderr,"%4d \n",t2loop);

    for (i3=0; i3<t2loop; i3++) {
        sf_floatread(uo,arr_size,in);
 
        for (i1=0; i1<arr_size; i1++){
            
            sum=0.0;
            klo=i1-m ; if(klo<0) klo=0;
            khi=i1+m ; if(khi>=arr_size) khi=arr_size;
            for (i2=klo ; i2<=khi ; i2++){
                k=i2;
                sum= sum+ uo[k];
            }
            sum *= 1.0/(2.0*m+1.0);
            der_tmp[i1]=sum;
        }
        
        for (i1=0; i1<arr_size; i1++){
            sum=0.0;
            klo=i1-m ; if(klo<0) klo=0;
            khi=i1+m ; if(khi>=arr_size) khi=arr_size;
            for (i2=klo ; i2<=khi ; i2++){
                k=i2;
                sum= sum+ (uo[k]-der_tmp[k])*(uo[k]-der_tmp[k]);                
            }
            sum *= 1.0/(2.0*m+1.0);
            if(fabsf(uo[i1]-der_tmp[i1])/sum >sigma){
                der[i1]=der_tmp[i1] +(uo[i1]-der_tmp[i1])*sigma*sum/(fabsf(uo[i1]-der_tmp[i1]));
            }else {
                der[i1]=uo[i1];
            }                                                             
            
        }
        sf_floatwrite(der,arr_size,out);
    }

    free (der); 
    free(der_tmp); 
    free(uo);
    exit(0);
}

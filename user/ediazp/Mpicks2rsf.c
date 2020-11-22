/* Creates a mask from horizons:
   
   horizon format:

   x1 h1
   x2 h2
   x3 h3
   x4 h4
    .
    .
    .
   xn hn



   xn> ... >x4 >x3 >x2>x1

   picks (file)  ascii file with two columns (x and h(x))
                 the x values must be increasing order,
                 you can easily achieve that by doing:

                 sort -k 1  unsorted_picks.txt > sorted_picks.txt
            
   stdin             2D file from which the axes will be read
   extend [false]    Extends picks to the boundaries of the axis
                          n Do not extend
                          y Extend to boundary

   tmask [true]     write a mask (1 if z>h(x))
          false     put a tic on the horizon

   above [false] put 1 above the horizon
          true   put 1 below the horizon

   ntic [1]     works with tmask=false; put 1 around ntic grid points
                above and below the horizon.   

   stdout       It writes a file with the same dimensions as stdin 
                with a mask function, 1 below the horizon 0 above 
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
#include "linted.h"

int main (int argc, char* argv[])
{

    int i1,i2,i3,i4,lb,hb=0;

    int n1,n2,ntic;
    float o1,o2;
    float d1,d2;
    float x,z;
    float *xcoord,*zcoord,*interp, **mask;
    int nl;
    char  *PICKS;
    char line[80];
    FILE *fp; 
    bool tmask, extend,above;
    sf_axis m1,m2;
    sf_file in=NULL, out=NULL;



    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    //=====================================================    
    //Get parameters from command line:
    
    if (!sf_getint("ntic", &ntic)) ntic=5;
    if (!sf_getbool("tmask", &tmask)) tmask=true;
    if (!sf_getbool("extend", &extend)) extend=false;
    if (!sf_getbool("above", &above)) above=false;
 
    //=====================================================
    //Get parameters from stdin file

    //if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");


    if (NULL == (PICKS=sf_getstring("picks"))) sf_error("need pick file");
    


   /* parameters from input file*/
    m1=sf_iaxa(in,1);
    m2=sf_iaxa(in,2);

    n1=sf_n(m1);  o1=sf_o(m1); d1=sf_d(m1); 
    n2=sf_n(m2);  o2=sf_o(m2); d2=sf_d(m2); 

	sf_oaxa(out,m1,1);
	sf_oaxa(out,m2,2);





 
   fp = fopen(PICKS, "rt");
   /* open the file for reading */
   /* elapsed.dta is the name of the file */
   /* "rt" means open the file for reading text */

   nl=0;
   while(fgets(line, 80, fp) != NULL)
   {
	 /* get a line, up to 80 chars from fr.  done if NULL */
	 sscanf (line, "%f %f", &x,&z);

     if ( x>o2 && x < (o2+(n2-1)*d2)) {
        if (z>o1 && z < (o1+(n1-1)*d1)) {
            nl=nl+1;
        }
     }
   }
   fclose(fp);  /* close the file prior to exiting the routine */ 
   
    
   fp = fopen(PICKS, "rt");
   /* open the file for reading */
   /* elapsed.dta is the name of the file */
   /* "rt" means open the file for reading text */

   xcoord=sf_floatalloc(nl);
   zcoord=sf_floatalloc(nl);
   interp=sf_floatalloc(n2);
   nl=0;
   while(fgets(line, 80, fp) != NULL)
   {
 	 /* get a line, up to 80 chars from fr.  done if NULL */
    sscanf (line, "%f %f", &x,&z);

    if ( x>o2 && x < (o2+(n2-1)*d2)) {
        if (z>o1 && z < (o1+(n1-1)*d1)) {
            xcoord[nl]=x;
            zcoord[nl]=z;
            nl=nl+1;
         }
      }
    }
    fclose(fp);  /* close the file prior to exiting the routine */ 



   // Use GEE's linear interpolation:


    lint1_init(o2,d2,xcoord);

    lint1_interp(nl, n2, zcoord, interp,extend);


    //allocate output model file
    mask=sf_floatalloc2(n1,n2);

    if(tmask){
        if (above){

            for (i2=0; i2<n2 ; i2++ ){
                x=o1;
                for ( i1=0; i1<n1; i1++){
                    if(x<interp[i2]) {
                        mask[i2][i1]=1.0;
                    } else {
                        mask[i2][i1]=0.0;
                    }
                    if(interp[i2]==-999.0){
                        mask[i2][i1]=0.0;
                    }
                    x+=d1;
                }
            }
            
        }else{
            for (i2=0; i2<n2 ; i2++ ){
                x=o1;
                for ( i1=0; i1<n1; i1++){
                    if(x>interp[i2]) {
                        mask[i2][i1]=1.0;
                    } else {
                        mask[i2][i1]=0.0;
                    }
                    if(interp[i2]==-999.0){
                        mask[i2][i1]=0.0;
                    }
                    x+=d1;
                }
            }
        }
    }else{
    
        for (i2=0; i2<n2 ; i2++ ){
            x=o1;
            for ( i1=0; i1<n1; i1++){
                mask[i2][i1]=0.0;
            }
            if(interp[i2]!=-999.0){
                i3 = (int) ((interp[i2]-o1)/d1);
                if(i3-ntic<0) lb=0;
                if(i3+ntic >n1) hb=n1;
                for (i4=lb; i4<hb; i4++){
                    mask[i2][i4]=1.0;
                }
                x+=d1;
            }
        }
    }
    sf_floatwrite(mask[0],n1*n2,out);



   exit(0);
}

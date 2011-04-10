/* Add swell noise to the data.*/
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

#include <math.h>
#include <time.h>
#include <rsf.h>

int main (int argc, char* argv[])
{
    float den, inten, max, *dat, slope;
    int n1, n2, n3, point, i, j, k;
    int dx, tt, width, length, num, temp1, temp2, temp3, temp4, ii;
    bool rep, noise;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input ("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);
	/* get the trace length (n1) and the number of traces (n2) and n3*/

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");


    if (!sf_getfloat("den",&den)) den=10.;
    /* noise density (percent, default=10, Min=0, Max=100) */
    if ( den < 0. ){
       den =0.;
    } else {
           if ( den > 100.) {
              den = 100.;
           }
    }
    den *= 0.01;

    if (!sf_getfloat("inten",&inten)) inten=0.1;
    /* noise intensity (multiple peak value of data, default=0.1) */

    if (!sf_getfloat("slope",&slope)) slope=0.1;
    /* noise slope (default=0.1) */
 
   if (!sf_getint("width",&width)) width=4;
    /* max noise width (default=4) */

   if (!sf_getint("length",&length)) length=30;
    /* max noise length (default=30) */

   if (!sf_getint("num",&num)) num=5;
    /* noise number (default=5) */

    if (!sf_getbool("rep",&rep)) rep=false;
    /* if y, replace data with noise */

    if (!sf_getbool("noise",&noise)) noise=false;
    /* if y, output noise only */

    dat = sf_floatalloc (n1*n2);

    srand((unsigned)time(0));
    for(ii=0;ii<n3;ii++)
    {

    if (noise) {
        max=0.;
	sf_floatread(dat,n1*n2,in);
        for (i=0; i<n1*n2; i++) {
            if (max < fabs(dat[i])) {
               max = fabs(dat[i]);
            }
        }
        for (i=0; i<n1*n2; i++) {
             dat[i] = 0.;
        }
        for (i=0; i< num; i++) {
            tt = rand()%n2;
            temp1 = rand()%100;
            temp2 = temp1+rand()%(n1/4);
            for (j=temp1; j< temp2; j++) {
                point= (int)(j*1.0/den)+rand()%10;
                dx=(int)(j*slope);
                temp4 = -1*width;
                for(k=temp4; k<(temp4+rand()%(width)); k++){
                    temp3 = point + rand()%(length);
                    if((tt+k+dx)<n2 && temp3<n1 && (tt+k+dx)>0 && temp3>0){
                        dat[(tt+k+dx)*n1+temp3] += inten*max*(0.01*(rand()%100)*2.-1.);
                    }
                    temp3=0;
                }
                temp4=0;
                point=0;
            }
            temp1=0;
            temp2=0;
        }

    } else {
        max=0.;
	sf_floatread(dat,n1*n2,in);
        for (i=0; i<n1*n2; i++) {
            if (max < fabs(dat[i])) {
               max = fabs(dat[i]);
            }
        }

        if (rep) {
            for (i=0; i< num; i++) {
                tt = rand()%n2;
                temp1 = rand()%100;
                temp2 = temp1+rand()%(n1/4);
                for (j=temp1; j< temp2; j++) {
                    point= (int)(j*1.0/den)+rand()%10;
                    dx=(int)(j*slope);
                    temp4 = -1*width;
                    for(k=temp4; k<(temp4+rand()%(width)); k++){
                        temp3 = point + rand()%(length);
                        if((tt+k+dx)<n2 && temp3<n1 && (tt+k+dx)>0 && temp3>0){
                            dat[(tt+k+dx)*n1+temp3] = inten*max*(0.01*(rand()%100)*2.-1.);
                        }
                        temp3=0;
                    }
                    temp4=0;
                    point=0;
                }
                temp1=0;
                temp2=0;
            }
        } else {
            for (i=0; i< num; i++) {
                tt = rand()%n2;
                temp1 = rand()%100;
                temp2 = temp1+rand()%(n1/4);
                for (j=temp1; j< temp2; j++) {
                    point= (int)(j*1.0/den)+rand()%10;
                    dx=(int)(j*slope);
                    temp4 = -1*width;
                    for(k=temp4; k<(temp4+rand()%(width)); k++){
                        temp3 = point + rand()%(length);
                        if((tt+k+dx)<n2 && temp3<n1 && (tt+k+dx)>0 && temp3>0){
                            dat[(tt+k+dx)*n1+temp3] += inten*max*(0.01*(rand()%100)*2.-1.);
                        }
                        temp3=0;
                    }
                    temp4=0;
                    point=0;
                }
                temp1=0;
                temp2=0;
            }
        }            
    }

    sf_floatwrite(dat,n1*n2,out);  
    }

    exit (0);
}

/* 	$Id$	 */

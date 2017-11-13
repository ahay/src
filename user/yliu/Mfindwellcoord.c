/* Find well location by using well coordinates. */
/*
  Copyright (C) 2017 Jilin University
  
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
#include <stdio.h>
#include <math.h>

int main (int argc, char* argv[]) 
{
    int nx1, nx2, nx3, ny1, ny2, ny3, n123; 
    int i1, i2, i3, coord1, coord2, coord3;

    float wx, wy, distance, temp;
    float *seisx, *seisy;
    sf_file xcoord, ycoord;

    sf_init (argc, argv); 
    xcoord = sf_input("xcoord");
    ycoord = sf_input("ycoord");
    
    if (!sf_histint(xcoord,"n1",&nx1)) sf_error("No n1= in input");
    if (!sf_histint(xcoord,"n2",&nx2)) sf_error("No n2= in input");
    nx3 = sf_leftsize(xcoord,2);

    if (!sf_histint(ycoord,"n1",&ny1)) sf_error("No n1= in ycoord");
    if (!sf_histint(ycoord,"n2",&ny2)) sf_error("No n2= in ycoord");
    ny3 = sf_leftsize(ycoord,2);

    if(nx1!=ny1 || nx2!=ny2 || nx3!=ny3) 
	sf_error("Need same size between input and ycoord file");
    
    if (!sf_getfloat("wellx",&wx)) sf_error("Need float well X coordinate");
    /* X coordinate for well */

    if (!sf_getfloat("welly",&wy)) sf_error("Need float well Y coordinate");
    /* Y coordinate for well */

    n123=nx1*nx2*nx3;

    seisx = sf_floatalloc(n123);
    seisy = sf_floatalloc(n123);

    sf_floatread(seisx,n123,xcoord);
    sf_floatread(seisy,n123,ycoord);

    coord1 = 0;
    coord2 = 0;
    coord3 = 0;
    distance = sqrtf((seisx[0]-wx)*(seisx[0]-wx)+(seisy[0]-wy)*(seisy[0]-wy));

    for(i3=0; i3 < nx3; i3++) {
	for(i2=0; i2 < nx2; i2++) {
	    for(i1=0; i1 < nx1; i1++) {
		temp = sqrtf((seisx[i3*nx2*nx1+i2*nx1+i1]-wx)*
			     (seisx[i3*nx2*nx1+i2*nx1+i1]-wx)+
			     (seisy[i3*nx2*nx1+i2*nx1+i1]-wy)*
			     (seisy[i3*nx2*nx1+i2*nx1+i1]-wy));
		if(temp<=distance) {
		    coord1 = i1;
		    coord2 = i2;
		    coord3 = i3;
		    distance = temp;
		}
	    }
	}
    }
                            
    printf("***************************************\n");
    printf("coordinate 1 in data grid = %d \n", coord1);
    printf("coordinate 2 in data grid = %d \n", coord2);
    printf("coordinate 3 in data grid = %d \n", coord3);
    printf("well x coordinate = %f \n", wx);
    printf("well y coordinate = %f \n", wy);
    printf("data x coordinate = %f \n", seisx[coord3*nx2*nx1+coord2*nx1+coord1]);
    printf("data y coordinate = %f \n", seisy[coord3*ny2*ny1+coord2*ny1+coord1]);
    printf("***************************************\n");
 
    exit (0);
}

/* 	$Id$	 */

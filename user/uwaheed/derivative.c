/* Final working version

 This program computes first derivative using 2nd and fourth order accuracy 
 To call this program use the following syntax:
 sf_derivative_2D(v,dtdz,dtdy,d2,d1,n2,n1,4);
                   or
 sf_derivative_3D(v,dtdz,dtdy,dx,d3,d2,d1,n3,n2,n1,4);

 v : field you want to take derivative of
 dtdz : derivative of the field v along axis 1
 dtdy : derivative of the field v along axis 2
 dtdx : derivative of the field v along axis 3
 d1 : Sampling along axis 1
 d2 : Sampling along axis 2
 d3 : Sampling along axis 3
 n1 : Number of grid points along axis 1
 n2 : Number of grid points along axis 2
 n3 : Number of grid pointsalong axis 3
 4 : its the accuracy order of the derivative you want (choose from 2,4)

*/



/*
  Copyright (C) 2013 KAUST
  
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
#include <rsf.h>

void sf_derivative_2D (float *t,float *dtdz, float *dtdx,
                       double d2, double d1,
                       int n2, int n1, int accuracy)

{
    int i1, i2;

    // For second order accuracy
    if (accuracy==2){

        if (n1<3 || n2<3){
            sf_error("Need more points for second order accuracy");
        }

        for(i2 = 0; i2 < n2; i2++){
            for(i1 = 1; i1 < n1-1; i1++){

            dtdz[n1*i2 + i1] = (t[n1*i2+(i1+1)]-t[n1*i2+(i1-1)])/(2*d1);

            }

        dtdz[n1*i2+0] = (-t[n1*i2 + 2] +4*t[n1*i2 + 1] -3*t[n1*i2 + 0])/(2*d1);

        dtdz[n1*i2+(n1-1)] = (3*t[n1*i2 + n1-1]- 4*t[n1*i2 + n1-2] + t[n1*i2 + n1-3])/(2*d1);

        }


        for(i1 = 0; i1 < n1; i1++){
            for(i2 = 1; i2 < n2-1; i2++){

            dtdx[n1*i2 + i1] = (t[n1*(i2+1)+i1] - t[n1*(i2-1)+i1])/(2*d2);

            }

        dtdx[n1*0+i1] = (-t[n1*2 + i1] +4*t[n1*1 + i1] -3*t[n1*0 + i1])/(2*d2);

        dtdx[n1*(n2-1)+i1] = (3*t[n1*(n2-1) + i1]- 4*t[n1*(n2-2) + i1] + t[n1*(n2-3) + i1])/(2*d2);

        }

    }


    //For fourth order accuracy
    if (accuracy==4){


        if (n1<5 || n2<5){
        sf_error("Need more points for fourth order accuracy");
        return;
        }
        for(i2 = 0; i2 < n2; i2++){
            for(i1 = 2; i1 < n1-2; i1++){

            dtdz[n1*i2 + i1] = (-t[n1*i2+(i1+2)] + 8*t[n1*i2+(i1+1)]-8*t[n1*i2+(i1-1)]+t[n1*i2+(i1-2)])/(12*d1);

            }

        dtdz[n1*i2+0] = (-3*t[n1*i2 + 4] + 16* t[n1*i2 + 3] -36* t[n1*i2 + 2] +48*t[n1*i2 + 1] -25*t[n1*i2 + 0])/(12*d1);

        dtdz[n1*i2+1] = (t[n1*i2 + 4] - 6* t[n1*i2 + 3] +18* t[n1*i2 + 2]-10*t[n1*i2 + 1] -3*t[n1*i2 + 0])/(12*d1);

        dtdz[n1*i2+(n1-1)] = (25*t[n1*i2 + n1-1]- 48*t[n1*i2 + n1-2] + 36*t[n1*i2 + n1-3] - 16*t[n1*i2 + n1-4] + 3*t[n1*i2 + n1-5])/(12*d1);

        dtdz[n1*i2+(n1-2)] = (3*t[n1*i2 + n1-1] + 10*t[n1*i2 + n1-2]- 18*t[n1*i2 + n1-3] + 6*t[n1*i2 + n1-4] - t[n1*i2 + n1-5])/(12*d1);

        }


        for(i1 = 0; i1 < n1; i1++){
            for(i2 = 2; i2 < n2-2; i2++){

            dtdx[n1*i2 + i1] = (-t[n1*(i2+2)+i1] + 8* t[n1*(i2+1)+i1] -8*t[n1*(i2-1)+i1] + t[n1*(i2-2)+i1])/(12*d2);

            }

        dtdx[n1*0+i1] = (-3*t[n1*4 + i1] +16*t[n1*3 + i1] -36*t[n1*2 + i1] +48*t[n1*1 + i1] -25*t[n1*0 + i1])/(12*d2);

        dtdx[n1*1+i1] = (t[n1*4 + i1] - 6*t[n1*3 + i1] + 18*t[n1*2 + i1] -10*t[n1*1 + i1] -3*t[n1*0 + i1])/(12*d2);


        dtdx[n1*(n2-1)+i1] = (25*t[n1*(n2-1) + i1]- 48*t[n1*(n2-2) + i1] + 36*t[n1*(n2-3) + i1]-16*t[n1*(n2-4) + i1]+3*t[n1*(n2-5) + i1])/(12*d2);

        dtdx[n1*(n2-2)+i1] = (3*t[n1*(n2-1) + i1] + 10*t[n1*(n2-2) + i1]- 18*t[n1*(n2-3) + i1] + 6*t[n1*(n2-4) + i1]- t[n1*(n2-5) + i1])/(12*d2);


        }

    }

   
}



void sf_derivative_3D (float *t,float *dtdz, 
                       float *dtdy, float *dtdx,
                       double d3, double d2, double d1,
                       int n3, int n2, int n1, int accuracy)

{
    int i1, i2, i3;
    int n12 = n1*n2;

    // For second order accuracy
    if (accuracy==2){

        if (n1<3 || n2<3 || n3<3){
            sf_error("Need more points for second order accuracy");
        }

        for(i3 = 0; i3 < n3; i3++){
            for(i2 = 0; i2 < n2; i2++){
                for(i1 = 1; i1 < n1-1; i1++){

                dtdz[n12*i3 + n1*i2 + i1] = (t[n12*i3+n1*i2+(i1+1)]-t[n12*i3+n1*i2+(i1-1)])/(2*d1);

                }

            dtdz[n12*i3+n1*i2+0] = (-t[n12*i3+n1*i2 + 2] +4*t[n12*i3+n1*i2 + 1] -3*t[n12*i3+n1*i2 + 0])/(2*d1);

            dtdz[n12*i3+n1*i2+(n1-1)] = (3*t[n12*i3+n1*i2 + n1-1]- 4*t[n12*i3+n1*i2 + n1-2] + t[n12*i3+n1*i2 + n1-3])/(2*d1);

            }
         }


        for(i3 = 0; i3 < n3; i3++){
            for(i1 = 0; i1 < n1; i1++){
                for(i2 = 1; i2 < n2-1; i2++){

                dtdy[n12*i3 + n1*i2 + i1] = (t[n12*i3 + n1*(i2+1) + i1] - t[n12*i3 + n1*(i2-1) + i1])/(2*d2);
                
                }

            dtdy[n12*i3 + n1*0 + i1] = (-t[n12*i3 + n1*2 + i1] +4*t[n12*i3 + n1*1 + i1] -3*t[n12*i3 + n1*0 + i1])/(2*d2);

            dtdy[n12*i3 + n1*(n2-1) + i1] = (3*t[n12*i3+n1*(n2-1) + i1]- 4*t[n12*i3+n1*(n2-2) + i1] + t[n12*i3+n1*(n2-3) + i1])/(2*d2);

            }

        }

        for(i2 = 0; i2 < n2; i2++){
            for(i1 = 0; i1 < n1; i1++){
                for(i3 = 1; i3 < n3-1; i3++){

                dtdx[n12*i3 + n1*i2 + i1] = (t[n12*(i3+1) + n1*i2 + i1] - t[n12*(i3-1) + n1*i2 + i1])/(2*d3);

                }

            dtdx[n12*0 + n1*i2 + i1] = (-t[n12*2 + n1*i2 + i1] +4*t[n12*1 + n1*i2 + i1] -3*t[n12*0 + n1*i2 + i1])/(2*d3);

            dtdx[n12*(n3-1) + n1*i2 + i1] = (3*t[n12*(n3-1)+n1*i2 + i1]- 4*t[n12*(n3-2)+n1*i2 + i1] + t[n12*(n3-3)+n1*i2 + i1])/(2*d3);

            }

        }
  
    }


    //For fourth order accuracy
    if (accuracy==4){


        if (n1<5 || n2<5 || n3<5){
        sf_error("Need more points for fourth order accuracy");
        }

        for(i3 = 0; i3 < n3; i3++){
            for(i2 = 0; i2 < n2; i2++){
                for(i1 = 2; i1 < n1-2; i1++){

                dtdz[n12*i3 + n1*i2 + i1] = (-t[n12*i3 + n1*i2 + (i1+2)] + 8*t[n12*i3 + n1*i2 + (i1+1)]-8*t[n12*i3 + n1*i2 + (i1-1)] + t[n12*i3+n1*i2+(i1-2)])/(12*d1);

                } // for i1

            dtdz[n12*i3+n1*i2+0] = (-3*t[n12*i3+n1*i2 + 4] + 16* t[n12*i3+n1*i2 + 3] -36* t[n12*i3+n1*i2 + 2] +48*t[n12*i3+n1*i2 + 1] -25*t[n12*i3+n1*i2 + 0])/(12*d1);

            dtdz[n12*i3+n1*i2+1] = (t[n12*i3+n1*i2 + 4] - 6* t[n12*i3+n1*i2 + 3] +18* t[n12*i3+n1*i2 + 2]-10*t[n12*i3+n1*i2 + 1] -3*t[n12*i3+n1*i2 + 0])/(12*d1);

            dtdz[n12*i3+n1*i2+(n1-1)] = (25*t[n12*i3+n1*i2 + n1-1]- 48*t[n12*i3+n1*i2 + n1-2] + 36*t[n12*i3+n1*i2 + n1-3] - 16*t[n12*i3+n1*i2 + n1-4] + 3*t[n12*i3+n1*i2 + n1-5])/(12*d1);

            dtdz[n12*i3+n1*i2+(n1-2)] = (3*t[n12*i3+n1*i2 + n1-1] + 10*t[n12*i3+n1*i2 + n1-2]- 18*t[n12*i3+n1*i2 + n1-3] + 6*t[n12*i3+n1*i2 + n1-4] - t[n12*i3+n1*i2 + n1-5])/(12*d1);
  
            }// for i2

        } // for i3

        for(i3 = 0; i3 < n3; i3++){
            for(i1 = 0; i1 < n1; i1++){
                for(i2 = 2; i2 < n2-2; i2++){

                dtdy[n12*i3+n1*i2 + i1] = (-t[n12*i3+n1*(i2+2)+i1] + 8* t[n12*i3+n1*(i2+1)+i1] -8*t[n12*i3+n1*(i2-1)+i1] + t[n12*i3+n1*(i2-2)+i1])/(12*d2);

                }// for i2

            dtdy[n12*i3+n1*0+i1] = (-3*t[n12*i3+n1*4 + i1] +16*t[n12*i3+n1*3 + i1] -36*t[n12*i3+n1*2 + i1] +48*t[n12*i3+n1*1 + i1] -25*t[n12*i3+n1*0 + i1])/(12*d2);

            dtdy[n12*i3+n1*1+i1] = (t[n12*i3+n1*4 + i1] - 6*t[n12*i3+n1*3 + i1] + 18*t[n12*i3+n1*2 + i1] -10*t[n12*i3+n1*1 + i1] -3*t[n12*i3+n1*0 + i1])/(12*d2);


            dtdy[n12*i3+n1*(n2-1)+i1] = (25*t[n12*i3+n1*(n2-1) + i1]- 48*t[n12*i3+n1*(n2-2) + i1] + 36*t[n12*i3+n1*(n2-3) + i1]-16*t[n12*i3+n1*(n2-4) + i1]+3*t[n12*i3+n1*(n2-5) + i1])/(12*d2);

            dtdy[n12*i3+n1*(n2-2)+i1] = (3*t[n12*i3+n1*(n2-1) + i1] + 10*t[n12*i3+n1*(n2-2) + i1]- 18*t[n12*i3+n1*(n2-3) + i1] + 6*t[n12*i3+n1*(n2-4) + i1]- t[n12*i3+n1*(n2-5) + i1])/(12*d2);

            } // for i1

        } // for i3

        for(i2 = 0; i2 < n2; i2++){
            for(i1 = 0; i1 < n1; i1++){
                for(i3 = 2; i3 < n3-2; i3++){

                dtdx[n12*i3+n1*i2 + i1] = (-t[n12*(i3+2)+n1*i2+i1] + 8* t[n12*(i3+1)+n1*i2+i1] -8*t[n12*(i3-1)+n1*i2+i1] + t[n12*(i3-2)+n1*i2+i1])/(12*d3);

                } // for i3

            dtdx[n12*0+n1*i2+i1] = (-3*t[n12*4+n1*i2 + i1] +16*t[n12*3+n1*i2 + i1] -36*t[n12*2+n1*i2 + i1] +48*t[n12*1+n1*i2 + i1] -25*t[n12*0+n1*i2 + i1])/(12*d3);

            dtdx[n12*1+n1*i2+i1] = (t[n12*4+n1*i2 + i1] - 6*t[n12*3+n1*i2 + i1] + 18*t[n12*2+n1*i2 + i1] -10*t[n12*1+n1*i2 + i1] -3*t[n12*0+n1*i2 + i1])/(12*d3);


            dtdx[n12*(n3-1)+n1*i2+i1] = (25*t[n12*(n3-1)+n1*i2 + i1]- 48*t[n12*(n3-2)+n1*i2 + i1] + 36*t[n12*(n3-3)+n1*i2 + i1]-16*t[n12*(n3-4)+n1*i2 + i1]+3*t[n12*(n3-5)+n1*i2 + i1])/(12*d3);

            dtdx[n12*(n3-2)+n1*i2+i1] = (3*t[n12*(n3-1)+n1*i2 + i1] + 10*t[n12*(n3-2)+n1*i2 + i1]- 18*t[n12*(n3-3)+n1*i2 + i1] + 6*t[n12*(n3-4)+n1*i2 + i1]- t[n12*(n3-5)+n1*i2 + i1])/(12*d3);

            } // for i1

        } // for i2

    } // if accuracy==4
   
} // sf_derivative_3D


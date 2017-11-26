/* Generate grid disk with gaussain amplitude for shooting portion of a wavefield.
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

#include <math.h>
#include <rsf.h>


void eurotate(float *obj      /* what to rotate*/,
              float *ref      /* reference direction*/,
              float angle    /* angle to rotate*/) {
/*<rotate obj with respect to a reference axis>*/

    float u,v,w,t,mat11,mat12,mat13,mat21,mat22,mat23,mat31,mat32,mat33;

    u = ref[0]; v = ref[1]; w = ref[2];
    t = angle*SF_PI/180;
/*
    mat11 = (u*u + (v*v+w*w)*cos(t))/(u*u+v*v+w*w);
    mat12 = (u*v*(1-cos(t))-w*sqrt(u*u+v*v+w*w)*sin(t))/(u*u+v*v+w*w);
    mat13 = (u*w*(1-cos(t))+v*sqrt(u*u+v*v+w*w)*sin(t))/(u*u+v*v+w*w);
    mat21 = (u*v*(1-cos(t))+w*sqrt(u*u+v*v+w*w)*sin(t))/(u*u+v*v+w*w);
    mat22 = (v*v + (u*u+w*w)*cos(t))/(u*u+v*v+w*w);
    mat23 = (v*w*(1-cos(t))-u*sqrt(u*u+v*v+w*w)*sin(t))/(u*u+v*v+w*w);
    mat31 = (u*w*(1-cos(t))-v*sqrt(u*u+v*v+w*w)*sin(t))/(u*u+v*v+w*w);
    mat32 = (v*w*(1-cos(t))+u*sqrt(u*u+v*v+w*w)*sin(t))/(u*u+v*v+w*w);
    mat33 = (w*w + (v*v+u*u)*cos(t))/(u*u+v*v+w*w);
    */
    
    mat11 = u*u*(1-cos(t))+ cos(t);
    mat12 = u*v*(1-cos(t))-w*sin(t);
    mat13 = u*w*(1-cos(t))+v*sin(t);
    mat21 = u*v*(1-cos(t))+w*sin(t);
    mat22 = v*v*(1-cos(t))+ cos(t);
    mat23 = v*w*(1-cos(t))-u*sin(t);
    mat31 = u*w*(1-cos(t))-v*sin(t);
    mat32 = v*w*(1-cos(t))+u*sin(t);
    mat33 = w*w*(1-cos(t))+ cos(t);
    
    obj[0] = mat11*obj[0] + mat12*obj[1] + mat13*obj[2]; 
    obj[1] = mat21*obj[0] + mat22*obj[1] + mat23*obj[2]; 
    obj[2] = mat31*obj[0] + mat32*obj[1] + mat33*obj[2]; 

    return;
}


int main(int argc, char* argv[])
{
    int n[3], na, nr, nt, i, j, k;
    float p[3], pnext[3], portho[3], d[3], o[3], c[3];
    float sd, radius, length, dt, ot, normp, normportho;
    float *amp,***wlt, *inpwlt,**disk;
    sf_file  gridin, sourcegrid, wavelet,waveletout ;

    sf_init (argc,argv);
    gridin = sf_input("in");
    wavelet = sf_input("wavelet");

    /* get 3-D grid parameters z,x,y */
    if (!sf_histint(gridin,"n1",n))     sf_error("No n1= in input"); // z
    if (!sf_histint(gridin,"n2",n+1))   sf_error("No n2= in input"); // x
    if (!sf_histint(gridin,"n3",n+2))   sf_error("No n3= in input"); // y
    if (!sf_histfloat(gridin,"d1",d))   sf_error("No d1= in input");
    if (!sf_histfloat(gridin,"d2",d+1)) sf_error("No d2= in input");
    if (!sf_histfloat(gridin,"d3",d+2)) sf_error("No d3= in input");
    if (!sf_histfloat(gridin,"o1",o))   o[0]=0.;
    if (!sf_histfloat(gridin,"o2",o+1)) o[1]=0.;
    if (!sf_histfloat(gridin,"o3",o+2)) o[2]=0.;
    
    if (!sf_histint(wavelet,"n1",&nt))     sf_error("No nt= in input");
    if (!sf_histfloat(wavelet,"d1",&dt))     sf_error("No dt= in input");
    if (!sf_histfloat(wavelet,"o1",&ot))     sf_error("No ot= in input");

    /* additional parameters */
    if(!sf_getint("na",&na)) na=180;
    /* number of angle */
    
    if(!sf_getint("nr",&nr)) nr=20;
    /* number of samples per radius */
    
    if(!sf_getfloat("radius",&radius)) radius=0.1;
    /* radius length */
    
    if(!sf_getfloat("sd",&sd)) sd=0.1;
    /* SD for gaussian amplitude */
    
    if (!sf_getfloats("direction",p,3)) sf_error("Need direction=");
    /* Direction of source  (x,y,z) */
    
    if (!sf_getfloats("center",c,3)) sf_error("Need center=");
    /* Position of source */


    /* specify output dimensions */
    sourcegrid = sf_output("out");
    waveletout = sf_output("waveletout");
    sf_putint (sourcegrid, "n1", 3);
    sf_putfloat (sourcegrid, "d1", 1.0);
    sf_putfloat (sourcegrid, "o1", 0.0);
    sf_putstring (sourcegrid, "label1", "Components");
    sf_putstring (sourcegrid, "unit1", "");
    sf_putint (sourcegrid, "n2", na*nr);
    sf_putfloat (sourcegrid, "o2", 0.0);
    sf_putfloat (sourcegrid, "d2", 1.0);
    sf_putstring (sourcegrid, "label2", "Number of sources");
    sf_putstring (sourcegrid, "unit2", "");
    sf_putint (sourcegrid, "n3", 1);
    sf_putfloat (sourcegrid, "o3", 0.0);
    sf_putfloat (sourcegrid, "d3", 1.0);
    sf_putstring (sourcegrid, "label3", "");
    sf_putstring (sourcegrid, "unit3", "");
    
    sf_putint (waveletout, "n1", 3);
    sf_putfloat (waveletout, "d1", 1.0);
    sf_putfloat (waveletout, "o1", 0.0);
    sf_putstring (waveletout, "label1", "Components");
    sf_putstring (waveletout, "unit1", "");
    sf_putint (waveletout, "n2", nt);
    sf_putfloat (waveletout, "d2", dt);
    sf_putfloat (waveletout, "o2", ot);
    sf_putstring (waveletout, "label2", "Time");
    sf_putstring (waveletout, "unit2", "s");
    sf_putint (waveletout, "n3", na*nr);
    sf_putfloat (waveletout, "o3", 0.0);
    sf_putfloat (waveletout, "d3", 1.0);
    sf_putstring (waveletout, "label3", "Number of sources");
    sf_putstring (waveletout, "unit3", "");

    /* Allocate space*/
    disk = sf_floatalloc2(3,na*nr); // different from whats need but to optimize speed !
    amp = sf_floatalloc(na*nr);
    wlt = sf_floatalloc3(3,nt,na*nr);
    inpwlt = sf_floatalloc(nt);
    
    /* Read input*/
    sf_floatread(inpwlt,nt,wavelet);

    /* Computation of disk-----------------------------------------------------------------------------------------------*/
    normp = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]); 
    p[0]/=normp; p[1]/=normp; p[2]/=normp;
    
    portho[0] = 0.0; // orthogonal to p (p cross x axis)
    portho[1] = p[2];
    portho[2] = -p[1];

    normportho = sqrt(portho[0]*portho[0]+portho[1]*portho[1]+portho[2]*portho[2]);
    if (normportho<0.001) { // x axis
        portho[0] = 0.0;
        portho[1] = 1.0;
        portho[2] = 0.0;
        normportho = 1.0;
    }

    pnext[0] = portho[0]/normportho; pnext[1] = portho[1]/normportho; pnext[2] = portho[2]/normportho; 
     
    for(i=0;i<na;i++) { // Number of angle rotations
        
        eurotate(pnext,p,360.0/na);
        for(j=0;j<nr;j++) { // Number of radius steps
        
            length = (j+1)*radius/nr;
            
            disk[j+i*nr][0] = pnext[0]*length + c[0]-o[1]; // x
            disk[j+i*nr][1] = pnext[1]*length + c[1]-o[2]; // y
            disk[j+i*nr][2] = pnext[2]*length + c[2]-o[0]; // z
            
            amp[j+i*nr] = (1/(sd*sqrt(2*SF_PI)))*exp(-length*length/(2*sd*sd)); // Gaussian amp
            
            for (k=0;k<nt;k++) { // Loop over time for wavelet modification
                wlt[j+i*nr][k][0] = amp[j+i*nr]*inpwlt[k];
            }
        }
    }


sf_floatwrite(disk[0],3*na*nr,sourcegrid);
sf_floatwrite(wlt[0][0],3*na*nr*nt,waveletout);

    exit (0);
}


/*         $Id$         */

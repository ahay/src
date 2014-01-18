/* Judgement of chaos  
   Input  - Complex;
   Output - Float
*/
/*
  Copyright (C) 2014 Jilin University
  
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

int main(int argc,char *argv[])
{
    int i,j,k,p,a,b,lx,ly,n1,n2,n3,fix,*w;
    float *x,*y,o2,d2,delta,max1,max2,min1,min2,gx,*output,*outmask;
    bool verb,fixgrid,ma;
    char none[]= " ";
    sf_complex *input;
    sf_file in, out , mask=NULL;
    
    sf_init (argc, argv); 
    in = sf_input("in");
    out = sf_output("out");
    if(NULL != sf_getstring ("mask"))
    { 
	mask = sf_output("mask");
	sf_setformat(mask,"float");
    }else{
	mask = NULL ;
    }
    
    sf_setformat(out,"float");
    
    
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1;
    if (!sf_histfloat(in,"o2",&o2)) o2=1;
    if (!sf_histfloat(in,"d2",&d2)) d2=1;
    n3 = sf_leftsize(in,2);
    
    /*judgement of input dimension*/
    if(!sf_getfloat("gx",&gx)) gx=2.0;
    /*Total Size of fixed grid*/
    if (!sf_getfloat("delta",&delta)) delta=0.01;
    /*The cell size of grid*/
    if(!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */ 
    if(!sf_getbool("fixgrid",&fixgrid)) fixgrid = false;
    /* if y ,the total size of grid determined by gx */
    if(!sf_getbool("ma",&ma)) ma = false;
    /* if y ,output auxilily file = mask*/ 
    
    sf_putint(out, "n1", n2);
    sf_putfloat(out, "d1", d2);
    sf_putfloat(out, "o1", o2);
    sf_putint(out, "n2", n3);
    sf_putint(out, "n3", 1);
    fix = (fabs(gx)*2)/delta;
    
    input = sf_complexalloc(n1*n2);
    x = sf_floatalloc(n1*n2);
    y = sf_floatalloc(n1*n2);
    output = sf_floatalloc(n2);
    w = sf_intalloc(fix*fix);

    /*get the realpart&imagepart of cmplex data*/
    
    for(k = 0; k < n3 ;k++)
    {
	if(fixgrid){
	    if(ma)
	    {
		outmask = sf_floatalloc(fix*fix);
		sf_putint(mask, "n1", fix);
		sf_putfloat(mask, "d1", delta);
		sf_putfloat(mask, "o1", -gx);
		sf_putint(mask, "n2", fix);
		sf_putfloat(mask, "d2", delta);
		sf_putfloat(mask, "o2", -gx);
		sf_putint(mask, "n3",n2);
		sf_putint(mask, "n4",n3);
		sf_putstring(mask,"label1",none);
		sf_putstring(mask,"label2",none);
		sf_putstring(mask,"unit1",none);
		sf_putstring(mask,"unit2",none);
	    }
	    sf_complexread(input,n1*n2,in);
	    
	    for(i=0;i<n1*n2;i++){
		x[i] = crealf(input[i]);
		y[i] = cimagf(input[i]);
	    }
	    
	    for(j = 0; j < n2;j++){
		for(i = 0; i < fix*fix; i++){
		    w[i]=0;
		    if(ma)outmask[i]=0;
		}
		for(i = 0; i < n1; i++){
		    a = (x[i+j*n1]+gx)/delta;
		    b = (y[i+j*n1]+gx)/delta;
		    w[b*fix+a]=1;
		    if(ma)outmask[a*fix+b]=1;
		}
		
		p = 0;
		for(i=0; i < (fix*fix); i++){
		    p=p+w[i];
		}
		output[j] = p;
		if(ma)sf_floatwrite(outmask,fix*fix,mask);
	    }	
	}else{
	    sf_complexread(input,n1*n2,in);
	    for(i=0;i<n1*n2;i++){
		x[i] = crealf(input[i]);
		y[i] = cimagf(input[i]);
	    }
	    
	    max1=-100,min1=100,max2=-100,min2=100;
	    for(i=0;i<n1*n2;i++){
		if(x[i]>max1)max1=x[i];
		if(x[i]<min1)min1=x[i];
		if(y[i]>max2)max2=y[i];
		if(y[i]<min2)min2=y[i];
	    }
	    /* determine the size of gird*/
	    lx = (max1-min1)/delta;
	    ly = (max2-min2)/delta;  
	    
	    w = sf_intalloc((lx+2)*(ly+2));
	    if(ma)
	    {
		outmask = sf_floatalloc((lx+2)*(ly+2));
		sf_putint(mask, "n1", ly+2);
		sf_putfloat(mask, "d1", delta);
		sf_putfloat(mask, "o1", min2);
		sf_putint(mask, "n2", lx+2);
		sf_putfloat(mask, "d2", delta);
		sf_putfloat(mask, "o2", min1);
		sf_putint(mask, "n3",n2);
		sf_putint(mask, "n4",n3);
		sf_putstring(mask,"label1",none);
		sf_putstring(mask,"label2",none);
		sf_putstring(mask,"unit1",none);
		sf_putstring(mask,"unit2",none);
	    }
	    for(j = 0; j < n2;j++){
		for(i=0;i<(lx+2)*(ly+2);i++){
		    w[i]=0;
		    if(ma)outmask[i]=0;
		}
		
		for(i=0;i<n1;i++){
		    a = (x[i+j*n1]-min1)/delta;
		    b = (y[i+j*n1]-min2)/delta;
		    w[(b+1)*(lx+2)+(a+1)]=1;
		    if(ma)outmask[(a+1)*(ly+2)+(b+1)]=1;
		}
		
		p = 0;
		for(i=0;i<((lx+2)*(ly+2));i++){
		    p=p+w[i];
		}
		output[j] = p;
		if(ma)sf_floatwrite(outmask,(lx+2)*(ly+2),mask);
	    }
	}
	
	sf_floatwrite(output,n2,out);
    }
    exit(0);
}

/* 	$Id$	 */

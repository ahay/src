/* Boundary for median filters*/
/*
  Copyright (C) 2007 University of Texas at Austin

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

#include "boundary.h"

void bound1(float* tempt,float* extendt,int nfw,int n1,int n2,bool boundary)
/*<extend seismic data>*/
{
    int m=(nfw-1)/2;
    int i,j;
    
    for(i=0;i<(n1+2*m)*(n2);i++){
	extendt[i]=0.0;
    }
    /*extend the number of samples*/
    for(i=0;i<n2;i++){
	for(j=0;j<m;j++){
	    if (boundary){
		extendt[(n1+2*m)*i+j]=tempt[n1*i+0];
	    }
	    else{
		extendt[(n1+2*m)*i+j]=0.0;
	    }
	}
    }
    for(i=0;i<n2;i++){
	for(j=0;j<n1;j++){
	    extendt[(n1+2*m)*i+j+m]=tempt[n1*i+j];
	}
    }
    for(i=0;i<n2;i++){
	for(j=0;j<m;j++){
	    if (boundary){
		extendt[(n1+2*m)*i+j+n1+m]=tempt[n1*i+n1-1];
	    }
	    else{
		extendt[(n1+2*m)*i+j+n1+m]=0.0;
	    }
	}
    }
    
}

void bound2(float* temp2,float* temp3,int n1,int tempnfw,int j,bool boundary)
/*<extend temporary seismic data>*/
{
    int k;
    /*extend trace*/
    if((j-tempnfw/2)>=0&&(j+tempnfw/2)<n1){
	for(k=0;k<tempnfw;k++){
	    temp3[k]=temp2[(j-tempnfw/2+k)];
	}
    }
    else if((j-tempnfw/2)<0&&(j+tempnfw/2)<n1){
	for(k=0;k<(abs(j-tempnfw/2));k++){
	    if (boundary){
		temp3[k]=temp2[0];
	    }
	    else{
		temp3[k]=0.0;
	    }
	}
	for(k=(abs(j-tempnfw/2));k<tempnfw;k++){
	    temp3[k]=temp2[k-abs(j-tempnfw/2)];
	}
    }
    else if((j-tempnfw/2)>=0&&(j+tempnfw/2)>=n1){
	for(k=0;k<(tempnfw-abs(j+tempnfw/2-n1+1));k++){
	    temp3[k]=temp2[j-tempnfw/2+k];
	}
	for(k=(tempnfw-abs(j+tempnfw/2-n1+1));k<tempnfw;k++){
	    if (boundary){
		temp3[k]=temp2[n1];
	    }
	    else{
		temp3[k]=0.0;
	    }
	}
    }
    else{
	for(k=0;k<(abs(j-tempnfw/2));k++){
	    if (boundary){
		temp3[k]=temp2[0];
	    }
	    else{
		temp3[k]=0.0;
	    }
	}
	for(k=(abs(j-tempnfw/2));k<(tempnfw-abs(j+tempnfw/2-n1+1));k++){
	    temp3[k]=temp2[k-abs(j-tempnfw/2)];
	}
	for(k=(tempnfw-abs(j+tempnfw/2-n1+1));k<tempnfw;k++){
	    if (boundary){
		temp3[k]=temp2[n1];
	    }
	    else{
		temp3[k]=0.0;
	    }
	}
    }	
}

void bound3(float* tempt,
	    float* extendt,
	    int nfw1        /* Sample direction*/,
	    int nfw2        /* Trace direction */,
	    int n1,int n2,
	    bool boundary)
/*<extend seismic data>*/
{
    int m1=(nfw1-1)/2;
    int m2=(nfw2-1)/2;
    int i,j;
    
    for(i=0;i<(n1+2*m1)*(n2+2*m2);i++){
	extendt[i]=0.0;
    }
    /*extend trace*/
    for(i=0;i<m2;i++){
	for(j=0;j<n1;j++){
	    if (boundary){
		extendt[(n1+2*m1)*i+j+m1]=tempt[n1*0+j];
	    }
	    else{
		extendt[(n1+2*m1)*i+j+m1]=0.0;
	    }
	}
    }
    for(i=0;i<n2;i++){
	for(j=0;j<n1;j++){
	    extendt[(n1+2*m1)*(i+m2)+j+m1]=tempt[n1*i+j];
	}
    }
    for(i=0;i<m2;i++){
	for(j=0;j<n1;j++){
	    if (boundary){
		extendt[(n1+2*m1)*(i+m2+n2)+j+m1]=tempt[n1*(n2-1)+j];
	    }
	    else{
		extendt[(n1+2*m1)*(i+m2+n2)+j+m1]=0.0;
	    }
	}
    }
    /*extend the number of samples*/
    for(i=0;i<(n2+2*m2);i++){
	for(j=0;j<m1;j++){
	    if (boundary){
		extendt[(n1+2*m1)*i+j]=extendt[(n1+2*m1)*i+m1];
	    }
	    else{
		extendt[(n1+2*m1)*i+j]=0.0;
	    }
	}
    }
    for(i=0;i<(n2+2*m2);i++){
	for(j=0;j<m1;j++){
	    if (boundary){
		extendt[(n1+2*m1)*i+j+n1+m1]=extendt[(n1+2*m1)*i+n1+m1-1];
	    }
	    else{
		extendt[(n1+2*m1)*i+j+n1+m1]=0.0;
	    }
	}
    }
    
}
/* 	$Id$	 */

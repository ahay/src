/* Normalize the data. 
Normalize by dividing the data set by its absolute maximum value. */
/*
  Copyright (C) 2013 University of Texas at Austin
  
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
#include <stdio.h>
#include <math.h>
#include <rsf.h>

int main(int argc, char* argv[])
{
  int  i, i1,n23,n1,n23a,n1a;
  float  max=0, *trace,*tracea;
  sf_file inp, outp,apply;
  sf_warning("input done");
  sf_init(argc,argv);
  outp = sf_output("out");

  if(!sf_getfloat("max",&max))
    {
      inp = sf_input("in");
      apply= sf_input("apply");
      sf_fileflush(outp,apply);
      if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
      if (!sf_histint(apply,"n1",&n1a)) sf_error("No n1= in apply");
      
      trace = sf_floatalloc(n1);
      tracea= sf_floatalloc(n1a);
      n23=sf_leftsize(inp,1);
      n23a=sf_leftsize(apply,1);
      max=0;
      for (i=0;i<n23;i++)
	{
	  sf_floatread(trace,n1,inp);
	  for (i1=0;i1<n1;i1++)
	    if(fabs(trace[i1])>max)max=fabs(trace[i1]);
	}    
      if(max==0){sf_error("At least one number is not equal to 0");}
	else{sf_warning("The absolute maximum is %f",max);}
    }
  else
    {
      apply=sf_input("apply");
      sf_fileflush(outp,apply);
      if (!sf_histint(apply,"n1",&n1a)) sf_error("No n1= in apply");
      tracea= sf_floatalloc(n1a);
      n23a=sf_leftsize(apply,1);
      
    }
  /*Apply normalization*/
  for (i=0;i<n23a;i++)
    {
      sf_floatread(tracea,n1a,apply);
      for (i1=0;i1<n1a;i1++)
	tracea[i1]=tracea[i1]/max;
      sf_floatwrite(tracea,n1a,outp);
    }
  exit(0);
}




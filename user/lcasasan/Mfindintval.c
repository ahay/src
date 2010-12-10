/* Find a certain integer value position in an array [n1]*/
/*
  Copyright (C) 2010 University of Texas at Austin

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
#include <string.h>
#include <rsf.h>

void find_val (const int* a/*in*/, 
				int n, 
				int val, 
				int val2,
				int* val_ind_pt /*out*/,
				const char* type /*comparison type*/);

int main (int argc, char* argv[])
{

	char *type;
	int n1,val,val2;	
	int *column=NULL; 
	int *index=NULL;	
	sf_file in=NULL,out=NULL;
	
	sf_init (argc,argv);
	
	in=sf_input("in");
	out=sf_output("out");
	sf_settype (out, SF_INT);
	
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
	if (!sf_getint("val",&val)) sf_error("No val= in input");


	type=sf_getstring("type"); /* type of comparison eq (=) leq (<=) geq(>=) */
	if (type==NULL) type="eq";
	
	if (!sf_getint("val2",&val2) ) val2=val;
	else type="ran"; 	

	if (val>val2) sf_error("val2=%d must be greater that val=%d",val,val2);
	
	column = sf_intalloc(n1);

    index = sf_intalloc(n1);

	sf_intread(column,n1,in);

	find_val (column, n1,val,val2, index, type);
	
	sf_intwrite(index,n1,out);	
	
	
    exit (0);
}

void find_val (const int* a/*in*/, 
				int n, 
				int val,
				int val2, 
				int* val_ind_pt /*out*/,
				const char* type /*comparison type*/) {

	int i;
        
	if (strcmp(type,"eq")==0) {
		sf_warning("eq");
			for (i = 1; i < n; i++)  {
                if (a[i] == val)                        
                	val_ind_pt[i]  = 1;
      			else 
				 	val_ind_pt[i]  = 0;
					
			}
	}

	if (strcmp(type,"leq")==0) {
		sf_warning("leq");
			for (i = 1; i < n; i++)  {
                if (a[i] <= val)                        
                	val_ind_pt[i]  = 1;
      			else 
				 	val_ind_pt[i]  = 0;
					
			}
	}

	if (strcmp(type,"geq")==0) {
		sf_warning("geq");
			for (i = 1; i < n; i++)  {
                if (a[i] >= val)                        
                	val_ind_pt[i]  = 1;
      			else 
				 	val_ind_pt[i]  = 0;
					
			}
	}

	if (strcmp(type,"ran")==0) {
		sf_warning("ran");
			for (i = 1; i < n; i++)  {
                if ((a[i] >= val)  && (a[i] <= val2))                      
                	val_ind_pt[i]  = 1;
      			else 
				 	val_ind_pt[i]  = 0;
					
			}
	}
        	

		
}

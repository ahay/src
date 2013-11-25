//  Compare the difference of two rsf data sets with the same size. 
//
//  Copyright (C) 2013 University of Texas at Austin
//  
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//  
//  This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//  
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//

#include <stdio.h>
#include <math.h>

#include <rsf.hh>
//#include "vecmatop.hh"

using namespace std;

int main(int argc, char* argv[])
{
    sf_init(argc,argv);

	iRSF par(0);			/* input parameters */
	iRSF inp1;				/* input 1 */
	iRSF inp2("match");		/* input 2 */
	oRSF dif;				/* output */

    int  i1, i2, n1, n2, n12;
    float  s=0;

	inp1.get("n1",n1);
	inp1.get("n2",n2);
	n12=n1*n2;

	std::valarray<float> pp1(n12);
	std::valarray<float> pp2(n12);	

	inp1 >> pp1;
	inp2 >> pp2;


    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    s = s + pow((pp1[i2*n1+i1]-pp2[i2*n1+i1]),2);
	}
    }

	dif.put("o2",0);
	dif.put("o1",0);
	dif.put("n1",1);
	dif.put("n2",1);	
    
    sf_warning("The difference is %f", s );
    dif << s;
    return 0;
}




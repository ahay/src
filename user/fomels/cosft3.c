/* 3-D cosine Fourier transform */
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
#include <rsf.h>

float cosft_dk(int n, float d)
/*< wavenumber sampling for cosine FT >*/
{
    if (1==n) return d;

    return 1./(2*kiss_fft_next_fast_size(n-1)*d);
}

void cosft3(bool inv,               /* forward or inverse */ 
	    int n1, int n2, int n3, /* dimensions */
	    float ***data           /* data [n3][n2][n1] */)
/*< 3-D transform (in place) >*/
{
    int i1, i2, i3;

    if (inv) {
	if (n3 > 1) {
	    sf_cosft_init(n3);
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    sf_cosft_inv(data[0][0],i1+i2*n1,n1*n2);
		}
	    }
	    sf_cosft_close();
	}
	
	if (n2 > 1) {
	    sf_cosft_init(n2);
	    for (i3=0; i3 < n3; i3++) {
		for (i1=0; i1 < n1; i1++) {
		    sf_cosft_inv(data[i3][0],i1,n1);
		}
	    }
	    sf_cosft_close();
	}
	
	if (n1 > 1) {
	    sf_cosft_init(n1);	
	    for (i3=0; i3 < n3; i3++) {
		for (i2=0; i2 < n2; i2++) {
		    sf_cosft_inv(data[i3][i2],0,1);
		}
	    }	
	    sf_cosft_close();
	}
    } else {
	if (n1 > 1) {
	    sf_cosft_init(n1);	
	    for (i3=0; i3 < n3; i3++) {
		for (i2=0; i2 < n2; i2++) {
		    sf_cosft_frw(data[i3][i2],0,1);
		}
	    }	
	    sf_cosft_close();
	}
	
	if (n2 > 1) {
	    sf_cosft_init(n2);
	    for (i3=0; i3 < n3; i3++) {
		for (i1=0; i1 < n1; i1++) {
		    sf_cosft_frw(data[i3][0],i1,n1);
		}
	    }
	    sf_cosft_close();
	}
	
	if (n3 > 1) {
	    sf_cosft_init(n3);
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    sf_cosft_frw(data[0][0],i1+i2*n1,n1*n2);
		}
	    }
	    sf_cosft_close();
	}
    }
}

void cosft12(bool inv,               /* forward or inverse */ 
	     int n1, int n2, int n3, /* dimensions */
	     float ***data           /* data [n3][n2][n1] */)
/*< transform first two coordinates (in place) >*/
{
    int i1, i2, i3;

    if (inv) {
	if (n3 > 1) {
	    sf_cosft_init(n3);
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    sf_cosft_inv(data[0][0],i1+i2*n1,n1*n2);
		}
	    }
	    sf_cosft_close();
	}
	
	if (n2 > 1) {
	    sf_cosft_init(n2);
	    for (i3=0; i3 < n3; i3++) {
		for (i1=0; i1 < n1; i1++) {
		    sf_cosft_inv(data[i3][0],i1,n1);
		}
	    }
	    sf_cosft_close();
	}
    } else {
	if (n2 > 1) {
	    sf_cosft_init(n2);
	    for (i3=0; i3 < n3; i3++) {
		for (i1=0; i1 < n1; i1++) {
		    sf_cosft_frw(data[i3][0],i1,n1);
		}
	    }
	    sf_cosft_close();
	}
	
	if (n3 > 1) {
	    sf_cosft_init(n3);
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    sf_cosft_frw(data[0][0],i1+i2*n1,n1*n2);
		}
	    }
	    sf_cosft_close();
	}
    }
}

void cosft2(bool inv,       /* forward or inverse */ 
	    int n1, int n2, /* dimensions */
	    float **data    /* data [n2][n1] */)
/*< 2-D transform (in place) >*/
{
    int i1, i2;

    if (inv) {
	if (n2 > 1) {
	    sf_cosft_init(n2);
	    for (i1=0; i1 < n1; i1++) {
		sf_cosft_inv(data[0],i1,n1);
	    }
	    sf_cosft_close();
	}	    

	if (n1 > 1) {
	    sf_cosft_init(n1);	
	    for (i2=0; i2 < n2; i2++) {
		sf_cosft_inv(data[i2],0,1);
	    }
	    sf_cosft_close();
	}
    } else {
	if (n1 > 1) {
	    sf_cosft_init(n1);	
	    for (i2=0; i2 < n2; i2++) {
		sf_cosft_frw(data[i2],0,1);
	    }
	    sf_cosft_close();
	}	

	if (n2 > 1) {
	    sf_cosft_init(n2);
	    for (i1=0; i1 < n1; i1++) {
		sf_cosft_frw(data[0],i1,n1);
	    }
	    sf_cosft_close();
	}
    }
}

void cosft1(bool inv,    /* forward or inverse */ 
	    int n1,      /* dimensions */
	    float *data  /* data [n3][n2][n1] */)
/*< 1-D transform (in place) >*/
{
    if (n1 > 1) {
	sf_cosft_init(n1);	
	if (inv) {
	    sf_cosft_inv(data,0,1);
	} else {
	    sf_cosft_frw(data,0,1);
	}
	sf_cosft_close();
    }
}

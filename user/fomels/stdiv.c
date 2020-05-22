/* Streaming division */
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

static int n;
static float lam;

void stdiv_init(int nd       /* data size */, 
		float lambda /* smoothing parameter */)
/*< initialize >*/
{
    n = nd;
    lam = lambda*lambda;
}

void stdiv (const float* num, const float* den,  float* rat)
/*< smoothly divide rat=num/den >*/
{
    int i;
    float prev;

    prev=0.0f;
    for (i=0; i < n; i++) {
	rat[i] = (num[i]*den[i]+lam*prev)/(den[i]*den[i]+lam);
	prev=rat[i];
    }    
}

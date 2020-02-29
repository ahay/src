/* 3D acoustic variable-velocity variable-density time-domain FD modeling 

The code uses a standard second-order stencil in time.
The coefficients of the spatial stencil are computed
by matching the transfer function of the 6-point discretized
first-derivative operator to the ideal response.

The code implements the linearized operator obtained from the
system of first-order PDEs parametrized in incompressibility and density

dv/dt = - 1./rho * grad(p)
dp/dt = - K * div(v)

where
  rho  : density
  K    : incompressibility
  div  : divergence operator
  grad : gradient  operator
  p,v    : pressure and particle velocity wavefields

The models supplied by the user are wave speed and density, the code performs
the conversion internally to buoyancy (inverse density) and incompressibility.

Author: Francesco Perrone
Date: November 2020
*/

/*
	Copyright (C) 2013 Colorado School of Mines

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
#include "prep_utils.h"
#include "kernels.h"
#include <time.h>

int main(int argc, char* argv[])
{

    exit (0);
}

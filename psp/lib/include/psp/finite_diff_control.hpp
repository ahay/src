/*
   Parallel Sweeping Preconditioner (PSP): a distributed-memory implementation
   of a sweeping preconditioner for 3d Helmholtz equations.

   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and
   The University of Texas at Austin

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef PSP_FINITE_DIFF_CONTROL_HPP
#define PSP_FINITE_DIFF_CONTROL_HPP 1

#include "elemental.hpp"

namespace psp {

enum BoundaryCondition { PML, DIRICHLET };
enum Stencil { SEVEN_POINT }; // TWENTY_SEVEN_POINT not yet supported

// Our control structure for defining the basic parameters for the problem 
// to be solved. The domain is assumed to be of the form:
//
//                 _______________ (wx,wy,0)
//                /              /|
//            x  /              / |
//  sweep dir   /              /  |
//     /\      /______________/   |
//     ||      |              |   |
//     ||      |              |   / (wx,wy,wz)
//     ||    z |              |  /  
//     ||      |              | /  
//             |______________|/
//          (0,0,wz)    y    (0,wy,wz)
//
// PML must be enforced at least on the bottom face, and the remaining faces 
// may be set as either PML or a zero Dirichlet boundary condition.
//
template<typename R>
struct FiniteDiffControl
{
    typedef typename elem::Complex<R> C;

    Stencil stencil; // only 7-point supported so far
    int nx, ny, nz;  // number of grid points in each direction
    R wx, wy, wz;    // width of the PML-padded box in each direction
    R omega;         // frequency of problem [rad/sec]

    R Cx, Cy, Cz;   // coefficient for PML in each direction
    int bx, by, bz; // number of grid points of PML in each direction

    R imagShift;           // stabilizing shift for preconditioner
    int cutoff;            // minimum acceptable leaf size (when depth=1)
    int numPlanesPerPanel; // number of xy planes to put into each panel

    BoundaryCondition frontBC;
    BoundaryCondition rightBC;
    BoundaryCondition backBC;
    BoundaryCondition leftBC;
    BoundaryCondition topBC;
    // The bottom boundary condition must be PML since we are sweeping from it.
};

} // namespace psp

#endif // PSP_FINITE_DIFF_CONTROL_HPP

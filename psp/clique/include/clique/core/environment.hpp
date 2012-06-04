/*
   Clique: a scalable implementation of the multifrontal algorithm

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
#ifndef CLIQUE_CORE_ENVIRONMENT_HPP
#define CLIQUE_CORE_ENVIRONMENT_HPP 1

#include "elemental.hpp"
#include <algorithm>
#include <map>

#include "clique/config.h"

namespace cliq {

typedef unsigned char byte;
 
typedef elem::Complex<float> scomplex;
typedef elem::Complex<double> dcomplex;

bool Initialized();
void Initialize( int& argc, char**& argv );
void Finalize();

#ifndef RELEASE
void PushCallStack( std::string s );
void PopCallStack();
void DumpCallStack();
#endif

// Pull in some of Elemental's imported libraries
namespace blas = elem::blas;
namespace lapack = elem::lapack;
namespace mpi = elem::mpi;

// Pull in a number of useful enums from Elemental
using namespace elem::unit_or_non_unit_wrapper;
using namespace elem::distribution_wrapper;
using namespace elem::orientation_wrapper;
using namespace elem::upper_or_lower_wrapper;
using namespace elem::left_or_right_wrapper;

// Pull in a few classes from Elemental
using elem::Complex;
using elem::Matrix;
using elem::Grid;
using elem::DistMatrix;

} // namespace cliq

#endif /* CLIQUE_CORE_ENVIRONMENT_HPP */

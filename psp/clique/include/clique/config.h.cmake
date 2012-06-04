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
#ifndef CLIQUE_CONFIG_HPP
#define CLIQUE_CONFIG_HPP 1

#define Clique_VERSION_MAJOR @Clique_VERSION_MAJOR@
#define Clique_VERSION_MINOR @Clique_VERSION_MINOR@

#cmakedefine USE_CUSTOM_ALLTOALLV_FOR_FACT
#cmakedefine USE_CUSTOM_ALLTOALLV_FOR_MULT
#cmakedefine USE_CUSTOM_ALLTOALLV_FOR_SOLVE

#endif /* CLIQUE_CONFIG_HPP */

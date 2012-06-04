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
#ifndef CLIQUE_NUMERIC_SYMM_FRONT_TREE_HPP
#define CLIQUE_NUMERIC_SYMM_FRONT_TREE_HPP 1

namespace cliq {

enum SolveMode { NORMAL_1D, NORMAL_2D, FAST_1D_LDL, FAST_2D_LDL };

inline bool 
ModeIs1d( SolveMode mode )
{
    if( mode == NORMAL_1D || mode == FAST_1D_LDL )
        return true;
    else
        return false;
}

namespace numeric {

// Only keep track of the left and bottom-right piece of the fronts
// (with the bottom-right piece stored in workspace) since only the left side
// needs to be kept after the factorization is complete.

template<typename F>
struct LocalSymmFront
{
    Matrix<F> frontL;
    mutable Matrix<F> work;
};

template<typename F>
struct LocalSymmFrontTree
{
    std::vector<LocalSymmFront<F> > fronts;
};

template<typename F>
struct DistSymmFront
{
    // The 'SolveMode' member variable of the parent 'DistSymmFrontTree' 
    // determines which of the following fronts is active.
    //   {NORMAL_1D,FAST_1D_LDL} -> front1d
    //   {NORMAL_2D,FAST_2D_LDL} -> front2d
    //
    // Split each front into a left and right piece such that the right piece
    // is not needed after the factorization (and can be freed).

    // TODO: Think about the fact that almost everything is now mutable...

    mutable DistMatrix<F,VC,STAR> front1dL;
    mutable DistMatrix<F,VC,STAR> work1d;

    mutable DistMatrix<F> front2dL;
    mutable DistMatrix<F> work2d;

    DistMatrix<F,VC,STAR> diag;
};

template<typename F>
struct DistSymmFrontTree
{
    SolveMode mode;
    std::vector<DistSymmFront<F> > fronts;
};

template<typename F>
struct SymmFrontTree
{
    LocalSymmFrontTree<F> local;
    DistSymmFrontTree<F> dist;
};

template<typename F>
void SetSolveMode( SymmFrontTree<F>& L, SolveMode solveMode );

} // namespace numeric
} // namespace cliq

#endif // CLIQUE_SYMM_FRONT_TREE_HPP

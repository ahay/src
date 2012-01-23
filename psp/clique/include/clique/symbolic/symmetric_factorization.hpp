/*
   Clique: a scalable implementation of the multifrontal algorithm

   Copyright (C) 2010-2011 Jack Poulson <jack.poulson@gmail.com>
   Copyright (C) 2011 Jack Poulson, Lexing Ying, and 
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
#ifndef CLIQUE_SYMBOLIC_SYMM_FACT_HPP
#define CLIQUE_SYMBOLIC_SYMM_FACT_HPP 1

namespace clique {
namespace symbolic {

struct LocalSymmOrigSupernode
{
    int size, offset;
    int parent; // -1 if root separator
    std::vector<int> children;
    std::vector<int> lowerStruct;
};

struct LocalSymmOrig
{
    std::vector<LocalSymmOrigSupernode> supernodes;
};

struct DistSymmOrigSupernode
{
    int size, offset;
    std::vector<int> lowerStruct;
};

struct DistSymmOrig
{
    mpi::Comm comm;
    std::vector<DistSymmOrigSupernode> supernodes;
};

struct SymmOrig
{
    LocalSymmOrig local;
    DistSymmOrig dist;
};

struct LocalSymmFactSupernode
{
    bool isLeftChild;
    int size, offset, myOffset;
    int parent; // -1 if root separator
    std::vector<int> children;

    std::vector<int> lowerStruct, origLowerStruct;
    std::vector<int> origLowerRelIndices;
    std::vector<int> leftChildRelIndices, rightChildRelIndices;
};

struct LocalSymmFact
{
    std::vector<LocalSymmFactSupernode> supernodes;
};

struct DistSymmFactSupernode
{
    mpi::Comm comm;
    Grid* grid;

    int size, offset, myOffset, leftChildSize, rightChildSize;
    std::vector<int> lowerStruct, origLowerStruct;

    // Useful for expanding sparse matrices into this frontal matrix
    std::vector<int> origLowerRelIndices;

    // The relative indices of the left and right children
    // (maps from the child update indices to our frontal indices).
    // These could be replaced with just the relative indices of our local 
    // submatrices of the child updates.
    std::vector<int> leftChildRelIndices, rightChildRelIndices;

    //
    // Helpers for the factorization
    //
    std::vector<int> numChildFactSendIndices;
    std::deque<int> leftChildFactColIndices, leftChildFactRowIndices,
                    rightChildFactColIndices, rightChildFactRowIndices;
    // This information does not necessarily have to be kept and can be
    // computed from the above information (albeit somewhat expensively).
    mutable std::vector<std::deque<int> > childFactRecvIndices;

    //
    // Helpers for 1d solves (few right-hand sides)
    //
    int localSize1d, localOffset1d;
    std::vector<int> numChildSolveSendIndices;
    std::deque<int> leftChildSolveIndices, rightChildSolveIndices;
    std::vector<std::deque<int> > childSolveRecvIndices;

    //
    // Helpers for 2d solves (many right-hand sides)
    //
    // TODO
};

struct DistSymmFact
{
    std::vector<DistSymmFactSupernode> supernodes;
};

struct SymmFact
{
    LocalSymmFact local;
    DistSymmFact dist;
};

void SymmetricFactorization
( const SymmOrig& orig, SymmFact& fact, bool storeFactRecvIndices=true );

void LocalSymmetricFactorization( const SymmOrig& orig, SymmFact& fact );

void DistSymmetricFactorization
( const SymmOrig& orig, SymmFact& fact, bool storeFactRecvIndices=true );

void ComputeFactRecvIndices
( const DistSymmFactSupernode& factSN,
  const DistSymmFactSupernode& factChildSN );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

inline void SymmetricFactorization
( const SymmOrig& orig, SymmFact& fact, bool storeFactRecvIndices )
{
#ifndef RELEASE
    PushCallStack("symbolic::SymmetricFactorization");
#endif
    LocalSymmetricFactorization( orig, fact );
    DistSymmetricFactorization( orig, fact, storeFactRecvIndices );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace symbolic
} // namespace clique

#endif /* CLIQUE_SYMBOLIC_SYMM_FACT_HPP */


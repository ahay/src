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
#include "clique.hpp"

void clique::symbolic::LocalSymmetricFactorization
( const SymmOrig& orig, SymmFact& fact )
{
#ifndef RELEASE
    PushCallStack("symbolic::LocalSymmetricFactorization");
#endif
    const int numSupernodes = orig.local.supernodes.size();
    fact.local.supernodes.resize( numSupernodes );

    // Perform the symbolic factorization
    int myOffset = 0;
    std::vector<int>::iterator it;
    std::vector<int> childrenStruct, partialStruct, fullStruct,
                     supernodeIndices;
    for( int s=0; s<numSupernodes; ++s )
    {
        const LocalSymmOrigSupernode& origSN = orig.local.supernodes[s];
        LocalSymmFactSupernode& factSN = fact.local.supernodes[s];
        factSN.size = origSN.size;
        factSN.offset = origSN.offset;
        factSN.myOffset = myOffset;
        factSN.parent = origSN.parent;
        factSN.children = origSN.children;
        factSN.origLowerStruct = origSN.lowerStruct;

        const int numChildren = origSN.children.size();
        const int numOrigLowerIndices = origSN.lowerStruct.size();
#ifndef RELEASE
        if( numChildren != 0 && numChildren != 2 )
            throw std::logic_error("Tree must be built from bisections");
#endif
        if( numChildren == 2 )
        {
            const int left = origSN.children[0];
            const int right = origSN.children[1];
            LocalSymmFactSupernode& leftChild = fact.local.supernodes[left];
            LocalSymmFactSupernode& rightChild = fact.local.supernodes[right];
            leftChild.isLeftChild = true;
            rightChild.isLeftChild = false;

            // Union the child lower structs
            const int numLeftIndices = leftChild.lowerStruct.size();
            const int numRightIndices = rightChild.lowerStruct.size();
#ifndef RELEASE
            for( int i=1; i<numLeftIndices; ++i )
            {
                if( leftChild.lowerStruct[i] <= leftChild.lowerStruct[i-1] )
                {
                    std::ostringstream msg;
                    msg << "Left child struct was not sorted for s=" << s 
                        << "\n";
                    throw std::logic_error( msg.str().c_str() );
                }
            }
            for( int i=1; i<numRightIndices; ++i )
            {
                if( rightChild.lowerStruct[i] <= rightChild.lowerStruct[i-1] )
                {
                    std::ostringstream msg;
                    msg << "Right child struct was not sorted for s=" << s
                        << "\n";
                    throw std::logic_error( msg.str().c_str() );
                }
            }
#endif
            const int childrenStructMaxSize = numLeftIndices + numRightIndices;
            childrenStruct.resize( childrenStructMaxSize );
            it = std::set_union
            ( leftChild.lowerStruct.begin(), leftChild.lowerStruct.end(),
              rightChild.lowerStruct.begin(), rightChild.lowerStruct.end(),
              childrenStruct.begin() );
            const int childrenStructSize = int(it-childrenStruct.begin());
            childrenStruct.resize( childrenStructSize );

            // Union the lower structure of this supernode

#ifndef RELEASE
            for( int i=1; i<numOrigLowerIndices; ++i )
            {
                if( origSN.lowerStruct[i] <= origSN.lowerStruct[i-1] )
                {
                    std::ostringstream msg;
                    msg << "Original struct was not sorted for s=" << s
                        << "\n";
                    throw std::logic_error( msg.str().c_str() );
                }
            }
#endif
            const int partialStructMaxSize = 
                childrenStructSize + numOrigLowerIndices;
            partialStruct.resize( partialStructMaxSize );
            it = std::set_union
            ( origSN.lowerStruct.begin(), origSN.lowerStruct.end(),
              childrenStruct.begin(), childrenStruct.end(),
              partialStruct.begin() );
            const int partialStructSize = int(it-partialStruct.begin());
            partialStruct.resize( partialStructSize );

            // Union again with the supernode indices
            supernodeIndices.resize( origSN.size );
            for( int i=0; i<origSN.size; ++i )
                supernodeIndices[i] = origSN.offset + i;
            fullStruct.resize( origSN.size + partialStructSize );
            it = std::set_union
            ( partialStruct.begin(), partialStruct.end(),
              supernodeIndices.begin(), supernodeIndices.end(),
              fullStruct.begin() );
            fullStruct.resize( int(it-fullStruct.begin()) );

            // Construct the relative indices of the original lower structure
            it = fullStruct.begin();
            factSN.origLowerRelIndices.resize( numOrigLowerIndices );
            for( int i=0; i<numOrigLowerIndices; ++i )
            {
                const int index = origSN.lowerStruct[i];
                it = std::lower_bound ( it, fullStruct.end(), index );
                factSN.origLowerRelIndices[i] = int(it-fullStruct.begin());
            }

            // Construct the relative indices of the children
            factSN.leftChildRelIndices.resize( numLeftIndices );
            it = fullStruct.begin();
            for( int i=0; i<numLeftIndices; ++i )
            {
                const int index = leftChild.lowerStruct[i];
                it = std::lower_bound( it, fullStruct.end(), index );
                factSN.leftChildRelIndices[i] = int(it-fullStruct.begin());
            }
            factSN.rightChildRelIndices.resize( numRightIndices );
            it = fullStruct.begin();
            for( int i=0; i<numRightIndices; ++i )
            {
                const int index = rightChild.lowerStruct[i];
                it = std::lower_bound( it, fullStruct.end(), index );
                factSN.rightChildRelIndices[i] = int(it-fullStruct.begin());
            }

            // Form lower struct of this supernode by removing supernode indices
            // (which take up the first origSN.size indices of fullStruct)
            const int lowerStructSize = fullStruct.size()-origSN.size;
            factSN.lowerStruct.resize( lowerStructSize );
            for( int i=0; i<lowerStructSize; ++i )
                factSN.lowerStruct[i] = fullStruct[origSN.size+i];
        }
        else // numChildren == 0, so this is a leaf supernode 
        {
            factSN.lowerStruct = origSN.lowerStruct;
            
            // Construct the trivial relative indices of the original structure
            factSN.origLowerRelIndices.resize( numOrigLowerIndices );
            for( int i=0; i<numOrigLowerIndices; ++i )
                factSN.origLowerRelIndices[i] = i + factSN.size;
        }

        myOffset += factSN.size;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}


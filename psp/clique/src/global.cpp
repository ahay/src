/*
   Clique: a scalable implementation of the multifrontal algorithm

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

namespace { 
bool cliqueInitializedElemental; 
bool initializedClique = false;
}

bool clique::Initialized()
{ return ::initializedClique; }

void clique::Initialize( int& argc, char**& argv )
{
    // If Clique has already been initialized, this is a no-op
    if( ::initializedClique )
        return;

    const bool mustInitElemental = !elemental::Initialized();
    if( mustInitElemental )
    {
        elemental::Initialize( argc, argv );
        ::cliqueInitializedElemental = true;
    }
    else
    {
        ::cliqueInitializedElemental = false;
    }
    ::initializedClique = true;
}

void clique::Finalize()
{
    // If Clique is not currently initialized, then this is a no-op
    if( !::initializedClique )
        return;
    
    if( ::cliqueInitializedElemental )
        elemental::Finalize();

    ::initializedClique = false;
}

#ifndef RELEASE
namespace { std::stack<std::string> callStack; }

void clique::PushCallStack( std::string s )
{ ::callStack.push( s ); }

void clique::PopCallStack()
{ ::callStack.pop(); }

void clique::DumpCallStack()
{
    std::ostringstream msg;
    while( !::callStack.empty() )
    {
        msg << "[" << ::callStack.size() << "]: " << ::callStack.top() << "\n";
        ::callStack.pop();
    }
    std::cerr << msg.str() << std::endl;
}
#endif // ifndef RELEASE

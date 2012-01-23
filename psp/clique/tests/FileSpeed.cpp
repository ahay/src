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
using namespace elemental;

void Usage()
{
    std::cout << "FileSpeed <num MB> <num pieces> <num files> <base name>\n"
              << "  <num MB>:     number of megabytes to write to file\n"
              << "  <num pieces>: number of pieces to break it into\n"
              << "  <num files>:  number of files to read/write\n"
              << "  <base name>:  base file name to write/read\n"
              << std::endl;
}

int
main( int argc, char* argv[] )
{
    clique::Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    if( argc < 5 )
    {
        if( commRank == 0 )        
            Usage();
        clique::Finalize();
        return 0;
    }

    int argNum = 1;
    const unsigned numMB = atoi( argv[argNum++] );
    const unsigned numPieces = atoi( argv[argNum++] );
    const unsigned numFiles = atoi( argv[argNum++] );
    const char* baseName = argv[argNum++];
    if( numPieces == 0 && commRank == 0 )
    {
        std::cout << "Number of pieces must be positive." << std::endl;
        clique::Finalize();
        return 0;
    }
    if( numFiles == 0 && commRank == 0 )
    {
        std::cout << "Number of files must be positive." << std::endl;
        clique::Finalize();
        return 0;
    }
    if( commRank == 0 )
        std::cout << "numMB:     " << numMB << "\n"
                  << "numPieces: " << numPieces << "\n"
                  << "numFiles:  " << numFiles << "\n"
                  << "baseName:  " << baseName << "\n"
                  << std::endl;

    try
    {
        const std::size_t bufferSize = numMB<<20;
        std::vector<char> buffer( bufferSize );

        const std::size_t normalPieceSize = bufferSize / numPieces;
        const std::size_t lastPieceSize = 
            normalPieceSize + (bufferSize%numPieces);

        // Write the files
        for( unsigned f=0; f<numFiles; ++f )
        {
            std::ostringstream filename;
            filename << baseName << "-" << commRank << "-" << f << ".dat";

            mpi::Barrier( comm );
            const double writeStartTime = mpi::Time();
            std::ofstream outFile
            ( filename.str().c_str(), std::ios::out|std::ios::binary );
            for( unsigned i=0; i<numPieces-1; ++i )
                outFile.write( &buffer[i*normalPieceSize], normalPieceSize );
            outFile.write
            ( &buffer[(numPieces-1)*normalPieceSize], lastPieceSize );
            outFile.close();
            mpi::Barrier( comm );
            const double writeStopTime = mpi::Time();
            if( commRank == 0 )
                std::cout << "Write time: " << writeStopTime-writeStartTime
                          << " secs." << std::endl;
        }

        // Read the files
        for( unsigned f=0; f<numFiles; ++f )
        {
            std::ostringstream filename;
            filename << baseName << "-" << commRank << "-" << f << ".dat";

            mpi::Barrier( comm );
            const double readStartTime = mpi::Time();
            std::ifstream inFile
            ( filename.str().c_str(), std::ios::in|std::ios::binary );
            for( unsigned i=0; i<numPieces-1; ++i )
            inFile.read( &buffer[i*normalPieceSize], normalPieceSize );
            inFile.read
            ( &buffer[(numPieces-1)*normalPieceSize], lastPieceSize );
            inFile.close();
            mpi::Barrier( comm );
            const double readStopTime = mpi::Time();
            if( commRank == 0 )
                std::cout << "Read time: " << readStopTime-readStartTime
                          << " secs." << std::endl;
        }
    }
    catch( std::exception& e )
    {
#ifndef RELEASE
        elemental::DumpCallStack();
        clique::DumpCallStack();
#endif
        std::ostringstream msg;
        msg << "Process " << commRank << " caught message:\n"
            << e.what() << "\n";
        std::cerr << msg.str() << std::endl;
    }

    clique::Finalize();
    return 0;
}

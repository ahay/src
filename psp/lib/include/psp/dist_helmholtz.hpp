/*
   Parallel Sweeping Preconditioner (PSP): a distributed-memory implementation
   of a sweeping preconditioner for 3d Helmholtz equations.

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
#ifndef PSP_DIST_HELMHOLTZ_HPP
#define PSP_DIST_HELMHOLTZ_HPP 1

namespace psp {

template<typename R>
class DistHelmholtz
{
public:
    typedef std::complex<R> C;

    DistHelmholtz
    ( const FiniteDiffControl<R>& control, elemental::mpi::Comm comm );

    ~DistHelmholtz();

    // Build the sparse matrix and the preconditioner
    void Initialize( const GridData<R>& velocity, bool accelerate=true );

    // Solves an O(1) set of linear systems with the sweeping preconditioner
    void SolveWithSQMR( GridData<C>& B, R bcgRelTol=1e-4 ) const;
    void SolveWithGMRES( GridData<C>& B, int m=20, R relTol=1e-3 ) const;

    // Destroy the sparse matrix and the preconditioner
    void Finalize();

private:
    elemental::mpi::Comm comm_;
    unsigned log2CommSize_;
    const FiniteDiffControl<R> control_;
    const R hx_, hy_, hz_; // grid spacings

    // Whether or not we have used velocity information to set up a 
    // preconditioner
    bool initialized_; 

    //
    // Solve helpers
    //

    void PullRightHandSides
    ( const GridData<C>& gridB, elemental::Matrix<C>& B ) const;
    void PushRightHandSides
    ( GridData<C>& gridB, const elemental::Matrix<C>& B ) const;

    void InternalSolveWithGMRES
    ( elemental::Matrix<C>& B, int m, R relTol ) const;
    void InternalSolveWithSQMR
    ( elemental::Matrix<C>& B, R bcgRelTol ) const;

    bool CheckForNaN( R alpha ) const;
    bool CheckForNaN( C alpha ) const;
    bool CheckForNaN( const std::vector<R>& alpha ) const;
    bool CheckForNaN( const std::vector<C>& alpha ) const;
    bool CheckForZero( const std::vector<R>& alpha ) const;
    bool CheckForZero( const std::vector<C>& alpha ) const;

    void Norms
    ( const elemental::Matrix<C>& xList, std::vector<R>& normList ) const;
    void InnerProducts
    ( const elemental::Matrix<C>& xList, const elemental::Matrix<C>& yList,
      std::vector<C>& alphaList ) const;

    void PseudoNorms
    ( const elemental::Matrix<C>& xList, std::vector<C>& pseudoNormList ) const;
    void PseudoInnerProducts
    ( const elemental::Matrix<C>& xList, const elemental::Matrix<C>& yList,
      std::vector<C>& alphaList ) const;

    void MultiplyColumns
    ( elemental::Matrix<C>& xList, const std::vector<C>& deltaList ) const;
    void DivideColumns
    ( elemental::Matrix<C>& xList, const std::vector<R>& deltaList ) const;
    void AddScaledColumns
    ( const std::vector<C>& deltaList,
      const elemental::Matrix<C>& xList, elemental::Matrix<C>& yList ) const;
    void SubtractScaledColumns
    ( const std::vector<C>& deltaList,
      const elemental::Matrix<C>& xList, elemental::Matrix<C>& yList ) const;

    void Multiply( elemental::Matrix<C>& B ) const;
    void Precondition( elemental::Matrix<C>& B ) const;

    // B_i := T_i B_i
    void SolvePanel( elemental::Matrix<C>& B, int i ) const;

    // B_{i+1} := B_{i+1} - A_{i+1,i} B_i
    void SubdiagonalUpdate( elemental::Matrix<C>& B, int i ) const;

    // Z := B_i, B_i := 0
    void ExtractPanel
    ( elemental::Matrix<C>& B, int i, elemental::Matrix<C>& Z ) const;

    // B_i := -A_{i,i+1} B_{i+1}
    void MultiplySuperdiagonal( elemental::Matrix<C>& B, int i ) const;

    // B_i := B_i + Z
    void UpdatePanel
    ( elemental::Matrix<C>& B, int i, const elemental::Matrix<C>& Z ) const;

    //
    // Information related to the decomposition of the domain into panels
    //

    // Information about the bottom (and full inner) panels
    int bottomDepth_;          // including the original PML
    int localBottomHeight_;    // local height of the bottom panel
    int localFullInnerHeight_; // local height of full inner panels
    int innerDepth_;           // total physical inner depth
    int numFullInnerPanels_;   // number of full inner panels we need
    int numPanels_;            // total number of panels

    // Information about the leftover inner panel
    bool haveLeftover_;
    int leftoverInnerDepth_;
    int localLeftoverInnerHeight_;

    // Information about the top panel
    int topOrigDepth_;
    int localTopHeight_;

    // Symbolic factorizations of each class of panels
    clique::symbolic::SymmFact 
        bottomSymbolicFact_, 
        leftoverInnerSymbolicFact_, 
        topSymbolicFact_;

    // Numeric factorizations of each panel
    clique::numeric::SymmFrontTree<C> bottomFact_;
    std::vector<clique::numeric::SymmFrontTree<C>*> fullInnerFacts_;
    clique::numeric::SymmFrontTree<C> leftoverInnerFact_;
    clique::numeric::SymmFrontTree<C> topFact_;

    //
    // Information related to the global sparse matrix
    //

    // Sparse matrix storage
    int localHeight_;
    std::vector<int> localToNaturalMap_;
    std::vector<int> localRowOffsets_;
    std::vector<int> owningProcesses_;
    std::vector<C> localEntries_;

    // Global sparse matrix communication information
    std::vector<int> globalSendCounts_, globalRecvCounts_; 
    std::vector<int> globalSendDispls_, globalRecvDispls_; 
    std::vector<int> globalSendIndices_; 

    // For redistributing pieces of B_i to form A_{i+1,i} B_i
    std::vector<int> subdiagSendCounts_, subdiagRecvCounts_;
    std::vector<int> subdiagSendDispls_, subdiagRecvDispls_;
    std::vector<int> subdiagPanelSendCounts_, subdiagPanelRecvCounts_;
    std::vector<int> subdiagPanelSendDispls_, subdiagPanelRecvDispls_;
    std::vector<int> subdiagSendIndices_;
    std::vector<int> subdiagRecvLocalIndices_, 
                     subdiagRecvLocalRows_;

    // For redistributing pieces of B_{i+1} to form A_{i,i+1} B_{i+1}
    std::vector<int> supdiagSendCounts_, supdiagRecvCounts_;
    std::vector<int> supdiagSendDispls_, supdiagRecvDispls_;
    std::vector<int> supdiagPanelSendCounts_, supdiagPanelRecvCounts_;
    std::vector<int> supdiagPanelSendDispls_, supdiagPanelRecvDispls_;
    std::vector<int> supdiagSendIndices_;
    std::vector<int> supdiagRecvLocalIndices_,
                     supdiagRecvLocalRows_;

    //
    // General helper routines 
    //

    C s1Inv( int x ) const;
    C s2Inv( int y ) const;
    C s3Inv( int v ) const;
    C s3InvArtificial( int v, int vOffset ) const;

    int NumLocalSupernodes( unsigned commRank ) const;
    static void NumLocalSupernodesRecursion
    ( int xSize, int ySize, int cutoff, 
      unsigned commRank, unsigned depthTilSerial, int& numLocal );

    //
    // Global sparse helper routines
    //

    clique::numeric::SymmFrontTree<C>& 
    PanelNumericFactorization( int whichPanel );

    const clique::numeric::SymmFrontTree<C>&
    PanelNumericFactorization( int whichPanel ) const;

    clique::symbolic::SymmFact&
    PanelSymbolicFactorization( int whichPanel );
    const clique::symbolic::SymmFact&
    PanelSymbolicFactorization( int whichPanel ) const;

    int PanelPadding( int whichPanel ) const;
    int PanelDepth( int whichPanel ) const;
    int WhichPanel( int v ) const;
    int PanelV( int whichPanel ) const;
    int LocalV( int v ) const;
    int LocalPanelOffset( int whichPanel ) const;
    int LocalPanelHeight( int whichPanel ) const;

    void GetGlobalVelocity
    ( const GridData<R>& velocity,
      std::vector<R>& myGlobalVelocity,
      std::vector<int>& offsets ) const;

    void FormGlobalRow( R alpha, int x, int y, int v, int row );

    int LocalPanelHeight( int vSize, int vPadding, unsigned commRank ) const;
    static void LocalPanelHeightRecursion
    ( int xSize, int ySize, int vSize, int vPadding, int cutoff, 
      unsigned teamRank, unsigned depthTilSerial, int& localHeight );

    // For building localToNaturalMap_, which takes our local index in the 
    // global sparse matrix and returns the original 'natural' index.
    void MapLocalPanelIndices( unsigned commRank, int whichPanel );
    static void MapLocalPanelIndicesRecursion
    ( int nx, int ny, int nz, int xSize, int ySize, int vSize, int vPadding,
      int xOffset, int yOffset, int vOffset, int cutoff,
      unsigned teamRank, unsigned depthTilSerial,
      std::vector<int>& localToNaturalMap, std::vector<int>& localRowOffsets,
      int& localOffset );

    void MapLocalConnectionIndices
    ( unsigned commRank, 
      std::vector<int>& localConnections, int whichPanel ) const;
    static void MapLocalConnectionIndicesRecursion
    ( int nx, int ny, int nz, int xSize, int ySize, int vSize, int vPadding,
      int xOffset, int yOffset, int vOffset, int cutoff,
      unsigned teamRank, unsigned depthTilSerial,
      std::vector<int>& localConnections, int& localOffset );

    int OwningProcess( int naturalIndex ) const;
    int OwningProcess( int x, int y, int v ) const;
    static void OwningProcessRecursion
    ( int x, int y, int vLocal, int xSize, int ySize, 
      unsigned depthTilSerial, int& process );

    //
    // Helpers for the PML-padded sparse-direct portion
    //

    void GetPanelVelocity
    ( int vOffset, int vSize, 
      const clique::symbolic::SymmFact& fact,
      const GridData<R>& velocity,
      std::vector<R>& myPanelVelocity,
      std::vector<int>& offsets,
      std::map<int,int>& panelNestedToNatural, 
      std::map<int,int>& panelNaturalToNested ) const;

    void FillPanelFronts
    ( int vOffset, int vSize,
      const clique::symbolic::SymmFact& symbFact,
            clique::numeric::SymmFrontTree<C>& fact,
      const GridData<R>& velocity,
      const std::vector<R>& myPanelVelocity,
            std::vector<int>& offsets,
            std::map<int,int>& panelNestedToNatural,
            std::map<int,int>& panelNaturalToNested ) const;

    void FormLowerColumnOfSupernode
    ( R alpha, int x, int y, int v, int vOffset, int vSize,
      int offset, int size, int j, 
      const std::vector<int>& origLowerStruct, 
      const std::vector<int>& origLowerRelIndices,
      std::map<int,int>& panelNaturalToNested, 
      std::vector<int>& frontIndices, std::vector<C>& values ) const;

    void LocalReordering( std::map<int,int>& reordering, int vSize ) const;
    static void LocalReorderingRecursion
    ( std::map<int,int>& reordering, int offset, 
      int xOffset, int yOffset, int xSize, int ySize, int vSize, int nx, int ny,
      int depthTilSerial, int cutoff, unsigned commRank );

    int ReorderedIndex( int x, int y, int vLocal, int vSize ) const;
    static int ReorderedIndexRecursion
    ( int x, int y, int vLocal, int xSize, int ySize, int vSize,
      int depthTilSerial, int cutoff, int offset );

    void FillOrigPanelStruct( int vSize, clique::symbolic::SymmOrig& S ) const;
    void FillDistOrigPanelStruct
    ( int vSize, int& nxSub, int& nySub, int& xOffset, int& yOffset, 
      clique::symbolic::SymmOrig& S ) const;
    void FillLocalOrigPanelStruct
    ( int vSize, int& nxSub, int& nySub, int& xOffset, int& yOffset, 
      clique::symbolic::SymmOrig& S ) const;
   
    // For use in FillLocalOrigPanelStruct
    struct Box
    {
        int parentIndex, nx, ny, xOffset, yOffset;
        bool leftChild;
    };
};

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename R>
inline 
DistHelmholtz<R>::~DistHelmholtz() 
{ }

} // namespace psp

#endif // PSP_DIST_HELMHOLTZ_HPP

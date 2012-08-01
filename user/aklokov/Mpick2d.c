/* 2D picking

Input:
	dataFile_.rsf - parameter spectrum (semblance) panel

Output:
	outFile_.rsf - picked optimal values
*/

/*
  Copyright (C) 2012 University of Texas at Austin
  
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

// VARIABLES

// files
sf_file dataFile_;          // input file - parameter spectrum file
sf_file sembFile_;          // auxiliary file - semblance
sf_file paramFile_;         // auxiliary file - max param values
sf_file outFile_;           // output file - picked values

// data
float*  panel_;             // parameter spectrum
float*  maxTrace_;       	// picked max param values

// param panel dimensions:
int     tNum_;                 
float   tStart_;
float   tStep_;

int     paramNum_;
float   paramStart_;
float   paramStep_;

int     xNum_;
float   xStart_;
float   xStep_;

// regularization params
int     xApp_;
float   eps_;

// FUNCTIONS

void TridiagonalSolve (const float *Ma, const float *Mb, float *Mc, float *Md, float *x, unsigned int n) {

    Mc[0] /= Mb[0];	
    Md[0] /= Mb[0];		
    for (int i = 1; i < n; ++i) {
        double id = (Mb[i] - Mc[i - 1] * Ma[i]);
		Mc[i] /= id;		
		Md[i] = (Md[i] - Md[i - 1] * Ma[i]) / id;
    }
 
    // back substitute
    x[n - 1] = Md[n - 1];
    for (int i = n - 2; i >= 0; --i)
        x[i] = Md[i] - Mc[i] * x[i + 1];

	return;	
}

int main (int argc, char* argv[]) {
   
// Initialize RSF 
    sf_init (argc, argv);

// Input files
    dataFile_   = sf_input("in");
// check that the input is float 
    if ( SF_FLOAT != sf_gettype (dataFile_) )   sf_error ("Need float input: dip-angle gathers");
    /* dip-angle gathers - stacks in the scattering-angle direction */
// Output auxiliary files
	char* paramTag = "maxForPicking.rsf";
    paramFile_ = sf_output (paramTag);
	char* sembTag = "sembForPicking.rsf";
	sembFile_ = sf_output (sembTag);
    outFile_ = sf_output ("out");

// Depth/time axis 
    if ( !sf_histint   (dataFile_, "n1", &tNum_) )   sf_error ("Need n1= in input");
    if ( !sf_histfloat (dataFile_, "d1", &tStep_) )  sf_error ("Need d1= in input");
    if ( !sf_histfloat (dataFile_, "o1", &tStart_) ) sf_error ("Need o1= in input");
// Dip angle axis 
    if ( !sf_histint   (dataFile_, "n2", &paramNum_) )   sf_error ("Need n2= in input");
    if ( !sf_histfloat (dataFile_, "d2", &paramStep_) )  sf_error ("Need d2= in input");
    if ( !sf_histfloat (dataFile_, "o2", &paramStart_) ) sf_error ("Need o2= in input");
// x axis 
    if ( !sf_histint   (dataFile_, "n3", &xNum_) )     sf_error ("Need n3= in input");
    if ( !sf_histfloat (dataFile_, "d3", &xStep_) )    sf_error ("Need d3= in input");
    if ( !sf_histfloat (dataFile_, "o3", &xStart_) )   sf_error ("Need o3= in input");

    if ( !sf_getfloat ("eps", &eps_) ) eps_ = 0;
    /* smoothness measure */
    if ( !sf_getint   ("xApp", &xApp_) ) xApp_ = 1;
    /* x-aperture */

// OUTPUT FILE
    sf_putint    (sembFile_, "n1", tNum_);     	 sf_putint    (sembFile_, "n2", xNum_); 
    sf_putfloat  (sembFile_, "d1", tStep_); 	 sf_putfloat  (sembFile_, "d2", xStep_); 
    sf_putfloat  (sembFile_, "o1", tStart_);	 sf_putfloat  (sembFile_, "o2", xStart_); 
    sf_putstring (sembFile_, "label1", "time");  sf_putstring (sembFile_, "label2", "inline");
    sf_putstring (sembFile_, "unit1", "s");      sf_putstring (sembFile_, "unit2", "m"); 

    sf_putint    (paramFile_, "n1", tNum_);      sf_putint    (paramFile_, "n2", xNum_); 
    sf_putfloat  (paramFile_, "d1", tStep_);   	 sf_putfloat  (paramFile_, "d2", xStep_); 
    sf_putfloat  (paramFile_, "o1", tStart_);	 sf_putfloat  (paramFile_, "o2", xStart_); 
    sf_putstring (paramFile_, "label1", "time"); sf_putstring (paramFile_, "label2", "inline");
    sf_putstring (paramFile_, "unit1", "s");     sf_putstring (paramFile_, "unit2", "m"); 

    sf_putint    (outFile_, "n1", tNum_);     	 sf_putint    (outFile_, "n2", xNum_); 
    sf_putfloat  (outFile_, "d1", tStep_);   	 sf_putfloat  (outFile_, "d2", xStep_); 
    sf_putfloat  (outFile_, "o1", tStart_);	     sf_putfloat  (outFile_, "o2", xStart_); 
    sf_putstring (outFile_, "label1", "time");   sf_putstring (outFile_, "label2", "inline");
    sf_putstring (outFile_, "unit1", "s");       sf_putstring (outFile_, "unit2", "m"); 

	const int panelSize = tNum_ * paramNum_;	
	panel_ = sf_floatalloc (panelSize);
	memset ( panel_, 0, panelSize * sizeof (float) );   

	float* maxTrace = sf_floatalloc (tNum_);
	float* paramTrace = sf_floatalloc (tNum_);

	// I. PICK MAX VALUES

	for (int ix = 0; ix < xNum_; ++ix) {
		sf_warning ("max picking: CIG %d of %d;", ix + 1, xNum_);

		const size_t startPos = (size_t) ix * panelSize * sizeof (float);
		sf_seek (dataFile_, startPos, SEEK_SET);
		sf_floatread (panel_, panelSize, dataFile_);

		memset ( maxTrace, 0, tNum_ * sizeof (float) );   
	
		float* pPanel = panel_;

		for (int ip = 0; ip < paramNum_; ++ip) {
			float* pMax   = maxTrace;			
			float* pParam = paramTrace;			
			const float curParam = paramStart_ + ip * paramStep_;
			for (int it = 0; it < tNum_; ++it, ++pMax, ++pPanel, ++pParam) {
				if (*pMax < *pPanel) { *pMax = *pPanel; *pParam = curParam; }							
			}
		}
		sf_floatwrite (maxTrace,   tNum_, sembFile_);
		sf_floatwrite (paramTrace, tNum_, paramFile_);
	}

	sf_warning (".");

	free (panel_);
	free (maxTrace);		
	free (paramTrace);		

	sf_fileclose (dataFile_);
	sf_fileclose (paramFile_);
	sf_fileclose (sembFile_);

	// II. REGULARIZATION

	// constants
	const int   halfXApp = xApp_ / 2;
	const float epsSq = eps_ * eps_;

	// input files
	sembFile_  = sf_input (sembTag);
	paramFile_ = sf_input (paramTag);

	for (int ix = 0; ix < xNum_; ++ix) {
		sf_warning ("regularization: CIG %d of %d;", ix + 1, xNum_);

        int xNum = xApp_;       
        size_t startInd = ix - halfXApp;

		// boundary checking
        int temp = ix - halfXApp;
		if (temp < 0) { xNum += temp; startInd -= temp; }
		temp = xNum_ - (ix + halfXApp) - 1;
		if (temp < 0) xNum += temp;

		// memory allocation
		const int curPanelSize = tNum_ * xNum; 

	    float* data = sf_floatalloc (curPanelSize);
	    float* semb = sf_floatalloc (curPanelSize);
		float*  res = sf_floatalloc (tNum_);
		float*   Ma = sf_floatalloc (tNum_);
		float*   Mb = sf_floatalloc (tNum_);    
		float*   Mc = sf_floatalloc (tNum_);    
		float* rightPart = sf_floatalloc (tNum_);

		// read data
		int offset = tNum_ * startInd * sizeof (float);
		sf_seek (paramFile_, offset, SEEK_SET);
		sf_seek (sembFile_ , offset, SEEK_SET);

		sf_floatread (data, curPanelSize, paramFile_);	
		sf_floatread (semb, curPanelSize, sembFile_);	

		// fix if semblance value equals zero
		float* pSemb = semb;
		for (int ip = 0; ip < curPanelSize; ++ip, ++pSemb) {
			if (*pSemb < 1e-6) *pSemb = 1e-6;
		}
	
		// fill out the equation matrices
		// - left part
		float* pMa = Ma;
		float* pMb = Mb;
		float* pMc = Mc;
		// --- basic - D^T D
		for (int it = 0; it < tNum_; ++it, ++pMa, ++pMb, ++pMc) {
			*pMa = -epsSq;
			*pMb = 2 * epsSq;
			*pMc = -epsSq;
		}
		// --- add weigths in the main diagonal
		pMb = Mb;		
		float* pR = rightPart;
		for (int it = 0; it < tNum_; ++it, ++pMb, ++pR) {
			float* pSemb = semb + it;			
			float sembStack = 0.f;		    
			for (int ix = 0; ix < xNum; ++ix, pSemb += tNum_)
				sembStack += *pSemb;
			*pMb += sembStack * sembStack;	
			*pR = sembStack;	
		}
		// - complete right part
		pR = rightPart;		
		for (int it = 0; it < tNum_; ++it, ++pR) {
			float* pSemb = semb + it;			
			float* pData = data + it;						
			float horStack = 0.f;		    
		    for (int ix = 0; ix < xNum; ++ix, pSemb += tNum_, pData += tNum_)
				horStack += *pData * (*pSemb);
			*pR *= horStack;
		}	
		// - correct end values
		Ma[0] = 0.f;
		const int lSamp = tNum_ - 1;
		Mb[lSamp] -= epsSq;
		Mc[lSamp] = 0.f;

		// solve the equation  
		TridiagonalSolve (Ma, Mb, Mc, rightPart, res, tNum_);

		// write down the result
		sf_floatwrite (res, tNum_, outFile_);
	
		// free memory
		free (Ma);
		free (Mb);
		free (Mc);

		free (rightPart);

		free (data);
		free (semb);
		free (res);
	}

	sf_warning (".");

	sf_fileclose (paramFile_);
	sf_fileclose (sembFile_);
	sf_fileclose (outFile_);

	return 0;
}

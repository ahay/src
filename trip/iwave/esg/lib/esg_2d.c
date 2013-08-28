/**
   2-4 stagger finite difference in 2D for isotropic elastic wave
 **/

// Version by Mario Bencomo (summer 2013)

#include "esgsteps_mod.h"
#include "esgn.h"

#ifdef _OPENMP
#include <omp.h>
#endif


//1/24 constant.
#define C24 4.166666666666666666666667e-2

// #define MARIO_VERBOSE
// #define SUPER_MARIO_VERBOSE

//	H = halo
//	B = boundary
// 		++++++++++++++++++++++++++++++++++++++++++++++++++
// 		+						 +
// 		+		        HII			 +
// 		+			 			 +
// 	+++++++++====================== BII =====================+++++++++
// 	+	||	|				|	||	 +
// 	+	||	|		II		|	||	 +
// 	+	||	|				|	||	 +
// 	+	||----------------------------------------------||	 +
// 	+	||	|				|	||	 +
// 	+	||	|				|	||	 +
// 	+	||	|				|	||	 +
// 	+ HIII BIII III	|	     physical		|  IV  BIV  HIV	 +
// 	+	||	|				|	||	 +
// 	+	||	|				|	||	 +
// 	+ 	||	|				|	||	 +
// 	+	||	|				|	||	 +
// 	+	------------------------------------------------||	 +
// 	+	||	|				|	||	 +
// 	+ 	||	|	       I		|	||	 +
// 	+	||	|				|	||	 +
// 	+++++++++==================== BI ========================+++++++++
// 		+						 +
// 		+		      HI			 +
// 		+						 +
// 		++++++++++++++++++++++++++++++++++++++++++++++++++
//    
//  	      y	^
//		^
//		^ >>>>
//		     x
//NOTE: px does not have stencil points in HI and HII
//      py does not have stencil points in HIII and HIV
//	sxy has stencil points in all halo regions
//	Also, vx,vy will not be updated at halo regions. 

//Declaration of auxilary functions
int esg_2d_4_ns( float ** _px,   float ** _py, 
		 float ** _mp00, float ** _mp01, 
	         float ** _vx,   float ** _vy,
		 float * _epx,   float * _epy,
		 float ** _px_x_I, float ** _px_x_II, float ** _px_x_III, float ** _px_x_IV,
		 float ** _px_y_I, float ** _px_y_II, float ** _px_y_III, float ** _px_y_IV,
		 float ** _py_x_I, float ** _py_x_II, float ** _py_x_III, float ** _py_x_IV,
		 float ** _py_y_I, float ** _py_y_II, float ** _py_y_III, float ** _py_y_IV,
		 int pml_I_width_y, int pml_II_width_y, int pml_III_width_x, int pml_IV_width_x,
		 int hlo_I_width_y, int hlo_II_width_y, int hlo_III_width_x, int hlo_IV_width_x,
		 int * n_,  int * n_V0,  int * n_V1,
		 int * gs_, int * gs_V0, int * gs_V1,
		 float la_x, float la_y, float dt      );

int esg_2d_4_ss0( float ** _sxy,  float ** _ms0, 
	          float ** _vx,   float ** _vy,
		  float * _epx,   float * _epy,
		  float ** _sxy_x_I, float ** _sxy_x_II, float ** _sxy_x_III, float ** _sxy_x_IV,
		  float ** _sxy_y_I, float ** _sxy_y_II, float ** _sxy_y_III, float ** _sxy_y_IV,
		  int pml_I_width_y, int pml_II_width_y, int pml_III_width_x, int pml_IV_width_x,
		  int hlo_I_width_y, int hlo_II_width_y, int hlo_III_width_x, int hlo_IV_width_x,
		  int * n_,  int * n_V0,  int * n_V1,
		  int * gs_, int * gs_V0, int * gs_V1,
		  float la_x, float la_y, float dt     );

int esg_2d_4_v0( float ** _vx,  float ** _mv0, 
	         float ** _px,  float ** _sxy,
		 float * _epx,  float * _epy,
		 float ** _vx_x_I, float ** _vx_x_II, float ** _vx_x_III, float ** _vx_x_IV,
		 float ** _vx_y_I, float ** _vx_y_II, float ** _vx_y_III, float ** _vx_y_IV,
		 int pml_I_width_y, int pml_II_width_y, int pml_III_width_x, int pml_IV_width_x,
		 int * n_,  int * n_P0,  int * n_S0,
		 int * gs_, int * gs_P0, int * gs_S0,
		 float la_x, float la_y, float dt     );

int esg_2d_4_v1( float ** _vy,  float ** _mv1, 
	         float ** _sxy, float ** _py,
		 float * _epx,  float * _epy,
		 float ** _vy_x_I, float ** _vy_x_II, float ** _vy_x_III, float ** _vy_x_IV,
		 float ** _vy_y_I, float ** _vy_y_II, float ** _vy_y_III, float ** _vy_y_IV,
		 int pml_I_width_y, int pml_II_width_y, int pml_III_width_x, int pml_IV_width_x,
		 int * n_,  int * n_S0,  int * n_P1,
		 int * gs_, int * gs_S0, int * gs_P1,
		 float la_x, float la_y, float dt     );



/*----------------------------------------------------------------------------*/
int esg_2d(RDOM *dom, int iarr, void *pars) {
/*----------------------------------------------------------------------------*/
#ifdef MARIO_VERBOSE
	fprintf(stderr,">> Inside esg_2d_4\n");
#endif

 	ESGN_TS_PARS * esgnpars = (ESGN_TS_PARS *) pars;
	if (esgnpars->k!=2){
		fprintf(stderr,"ERROR: inside esg_2d, order 2-%d not implemented for 2D\n",esgnpars->k*2);
		return 1;
	}

  	RARR * s       = dom->_s;
	RDOM * ld_pml  = ((ESGN_TS_PARS*)pars)->ld_pml;
	RARR * s_pml;

	register float ** restrict _px,  ** restrict _py;
	register float ** restrict _vx,	 ** restrict _vy;
	register float ** restrict _sxy;

	register float ** restrict _f1_x_I, ** restrict _f1_x_II, ** restrict _f1_x_III, ** restrict _f1_x_IV;
	register float ** restrict _f1_y_I, ** restrict _f1_y_II, ** restrict _f1_y_III, ** restrict _f1_y_IV;
	register float ** restrict _f2_x_I, ** restrict _f2_x_II, ** restrict _f2_x_III, ** restrict _f2_x_IV;
	register float ** restrict _f2_y_I, ** restrict _f2_y_II, ** restrict _f2_y_III, ** restrict _f2_y_IV;

	register float ** restrict _multi1, ** restrict _multi2;
	register float  * restrict _epx,     * restrict _epy;
	
	register int pml_I_width_y, pml_II_width_y, pml_III_width_x, pml_IV_width_x;
	register int hlo_I_width_y, hlo_II_width_y, hlo_III_width_x, hlo_IV_width_x;
	register float la_x, la_y, dt;
	IPNT n_,   n_0,   n_1;
	IPNT gs_,  gs_0,  gs_1;
	IPNT temp;
	int empty;

	la_x = ((ESGN_TS_PARS*)pars)->lam[0] * C24;
	la_y = ((ESGN_TS_PARS*)pars)->lam[1] * C24;
	dt = ((ESGN_TS_PARS*)pars)->dt;
	

	/*~~~ updating pressures, px and py ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  	if ( iarr == D_P1 )  return 0;
  	if ( iarr == D_P0 ){
		/* Initializing fields */
		_px     = s[D_P0]._s2;
		_py     = s[D_P1]._s2; 
		_vx     = s[D_V0]._s2;
		_vy     = s[D_V1]._s2;

		_multi1 = s[D_MP00]._s2;
		_multi2 = s[D_MP01]._s2;

		_epx    = s[D_EP[0]]._s1;	/* 1D */
		_epy    = s[D_EP[1]]._s1;  	/* 1D */

		/* pml fields */
		s_pml   = ld_pml[0]._s;
		_f1_x_I = s_pml[D_P0]._s2; 
		_f2_x_I = s_pml[D_P1]._s2;
		s_pml   = ld_pml[1]._s;
		_f1_y_I = s_pml[D_P0]._s2;
		_f2_y_I = s_pml[D_P1]._s2;
	
		s_pml    = ld_pml[2]._s;
		_f1_x_II = s_pml[D_P0]._s2; 
		_f2_x_II = s_pml[D_P1]._s2;
		s_pml    = ld_pml[3]._s;
		_f1_y_II = s_pml[D_P0]._s2;
		_f2_y_II = s_pml[D_P1]._s2;
	
		s_pml     = ld_pml[4]._s;
		_f1_x_III = s_pml[D_P0]._s2; 
		_f2_x_III = s_pml[D_P1]._s2;
		s_pml     = ld_pml[5]._s;
		_f1_y_III = s_pml[D_P0]._s2;
		_f2_y_III = s_pml[D_P1]._s2;
		
		s_pml = ld_pml[6]._s;
		_f1_x_IV = s_pml[D_P0]._s2;
		_f2_x_IV = s_pml[D_P1]._s2;
		s_pml = ld_pml[7]._s;
		_f1_y_IV = s_pml[D_P0]._s2;
		_f2_y_IV = s_pml[D_P1]._s2;

		/* getting starting indexes and size of fields */
		ra_gse( &s[D_P0], gs_ , temp);
		ra_gse( &s[D_V0], gs_0, temp);
		ra_gse( &s[D_V1], gs_1, temp);
		ra_size( &s[D_P0], n_  );
		ra_size( &s[D_V0], n_0 );
		ra_size( &s[D_V1], n_1 );

		//computing widths of pml regions
		/* pml region I */
		rd_empty( ld_pml+0, D_P0, &empty );
		if (empty) {
			pml_I_width_y = 0;
		}
		else {
			s_pml = ld_pml[0]._s;
			pml_I_width_y = s_pml[D_P0]._dims[1].n;
		}
	
		/* pml region II */
		rd_empty( ld_pml+2, D_P0, &empty );
		if (empty) {
			pml_II_width_y = 0;
		}
		else {
			s_pml = ld_pml[2]._s;;
			pml_II_width_y = s_pml[D_P0]._dims[1].n;
		}
	
		/* pml region III */
		rd_empty( ld_pml+4, D_P0, &empty );
		if (empty) {
			pml_III_width_x = 0;
		}
		else {
			s_pml = ld_pml[4]._s;
			pml_III_width_x = s_pml[D_P0]._dims[0].n;
		}
	
		/* pml region IV */
		rd_empty( ld_pml+6, D_P0, &empty );
		if (empty) {
			pml_IV_width_x = 0;
		}
		else {
			s_pml = ld_pml[6]._s;
			pml_IV_width_x = s_pml[D_P0]._dims[0].n;
		}
		
		// computing widths of halo regions
		hlo_I_width_y   = s[D_P1]._dims[1].gs - s[D_P1]._dims[1].gs0;
		hlo_II_width_y  = s[D_P1]._dims[1].ge0 - s[D_P1]._dims[1].ge;
		hlo_III_width_x = s[D_P0]._dims[0].gs - s[D_P0]._dims[0].gs0;
		hlo_IV_width_x  = s[D_P0]._dims[0].ge0 - s[D_P0]._dims[0].ge;
	
		return esg_2d_4_ns ( _px,     _py, 
		 		     _multi1, _multi2, 
	         		     _vx,     _vy,
		 		     _epx,    _epy,
		 		     _f1_x_I, _f1_x_II, _f1_x_III, _f1_x_IV,
		 		     _f1_y_I, _f1_y_II, _f1_y_III, _f1_y_IV,
		 		     _f2_x_I, _f2_x_II, _f2_x_III, _f2_x_IV,
		 		     _f2_y_I, _f2_y_II, _f2_y_III, _f2_y_IV,
		 		     pml_I_width_y, pml_II_width_y, pml_III_width_x, pml_IV_width_x,
				     hlo_I_width_y, hlo_II_width_y, hlo_III_width_x, hlo_IV_width_x,
		 		     n_,  n_0,  n_1,
		 		     gs_, gs_0, gs_1,
		 		     la_x, la_y, dt    );
	}
	/*~~~ updating shear stress sxy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  	if ( iarr == D_S0 ){
		/* Initializing fields */
		_sxy    = s[D_S0]._s2;
		_vx     = s[D_V0]._s2;
		_vy     = s[D_V1]._s2;

		_multi1 = s[D_MS0]._s2;

		_epx    = s[D_EV[0]]._s1;	/* 1D */
		_epy    = s[D_EV[1]]._s1;  	/* 1D */

		/* pml fields */
		s_pml   = ld_pml[0]._s;
		_f1_x_I = s_pml[D_S0]._s2; 
		s_pml   = ld_pml[1]._s;
		_f1_y_I = s_pml[D_S0]._s2;
	
		s_pml    = ld_pml[2]._s;
		_f1_x_II = s_pml[D_S0]._s2; 
		s_pml    = ld_pml[3]._s;
		_f1_y_II = s_pml[D_S0]._s2;
	
		s_pml     = ld_pml[4]._s;
		_f1_x_III = s_pml[D_S0]._s2; 
		s_pml     = ld_pml[5]._s;
		_f1_y_III = s_pml[D_S0]._s2;
		
		s_pml = ld_pml[6]._s;
		_f1_x_IV = s_pml[D_S0]._s2;
		s_pml = ld_pml[7]._s;
		_f1_y_IV = s_pml[D_S0]._s2;

		/* getting starting indexes and size of fields */
		ra_gse( &s[D_S0], gs_ , temp);
		ra_gse( &s[D_V0], gs_0, temp);
		ra_gse( &s[D_V1], gs_1, temp);
		ra_size( &s[D_S0], n_  );
		ra_size( &s[D_V0], n_0 );
		ra_size( &s[D_V1], n_1 );

		//computing widths of pml regions
		/* pml region I */
		rd_empty( ld_pml+0, D_S0, &empty );
		if (empty) {
			pml_I_width_y = 0;
		}
		else {
			s_pml = ld_pml[0]._s;
			pml_I_width_y = s_pml[D_S0]._dims[1].n;
		}
	
		/* pml region II */
		rd_empty( ld_pml+2, D_S0, &empty );
		if (empty) {
			pml_II_width_y = 0;
		}
		else {
			s_pml = ld_pml[2]._s;;
			pml_II_width_y = s_pml[D_S0]._dims[1].n;
		}
	
		/* pml region III */
		rd_empty( ld_pml+4, D_S0, &empty );
		if (empty) {
			pml_III_width_x = 0;
		}
		else {
			s_pml = ld_pml[4]._s;
			pml_III_width_x = s_pml[D_S0]._dims[0].n;
		}
	
		/* pml region IV */
		rd_empty( ld_pml+6, D_S0, &empty );
		if (empty) {
			pml_IV_width_x = 0;
		}
		else {
			s_pml = ld_pml[6]._s;
			pml_IV_width_x = s_pml[D_S0]._dims[0].n;
		}

		// computing widths of halo regions
		hlo_I_width_y   = s[D_S0]._dims[1].gs - s[D_S0]._dims[1].gs0;
		hlo_II_width_y  = s[D_S0]._dims[1].ge0 - s[D_S0]._dims[1].ge;
		hlo_III_width_x = s[D_S0]._dims[0].gs - s[D_S0]._dims[0].gs0;
		hlo_IV_width_x  = s[D_S0]._dims[0].ge0 - s[D_S0]._dims[0].ge;

		return esg_2d_4_ss0 ( _sxy, _multi1, 
	         		      _vx,   _vy,
		 		      _epx,  _epy,
		 		      _f1_x_I, _f1_x_II, _f1_x_III, _f1_x_IV,
		 		      _f1_y_I, _f1_y_II, _f1_y_III, _f1_y_IV,
		 		      pml_I_width_y, pml_II_width_y, pml_III_width_x, pml_IV_width_x,
		 		      hlo_I_width_y, hlo_II_width_y, hlo_III_width_x, hlo_IV_width_x,
		 		      n_,  n_0,  n_1,
		 		      gs_, gs_0, gs_1,
		 		      la_x, la_y, dt     );
	}
	/*~~~ updating velocity vx ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  	if ( iarr == D_V0 ){
		/* Initializing fields */
		_vx     = s[D_V0]._s2;
		_px     = s[D_P0]._s2;
		_sxy    = s[D_S0]._s2;

		_multi1 = s[D_MV0]._s2;

		_epx    = s[D_EV[0]]._s1;	/* 1D */
		_epy    = s[D_EP[1]]._s1;  	/* 1D */

		/* pml fields */
		s_pml   = ld_pml[0]._s;
		_f1_x_I = s_pml[D_V0]._s2; 
		s_pml   = ld_pml[1]._s;
		_f1_y_I = s_pml[D_V0]._s2;
	
		s_pml    = ld_pml[2]._s;
		_f1_x_II = s_pml[D_V0]._s2; 
		s_pml    = ld_pml[3]._s;
		_f1_y_II = s_pml[D_V0]._s2;
	
		s_pml     = ld_pml[4]._s;
		_f1_x_III = s_pml[D_V0]._s2; 
		s_pml     = ld_pml[5]._s;
		_f1_y_III = s_pml[D_V0]._s2;
		
		s_pml = ld_pml[6]._s;
		_f1_x_IV = s_pml[D_V0]._s2;
		s_pml = ld_pml[7]._s;
		_f1_y_IV = s_pml[D_V0]._s2;

		/* getting starting indexes and size of fields */
		ra_gse( &s[D_V0], gs_ , temp);
		ra_gse( &s[D_P0], gs_0, temp);
		ra_gse( &s[D_S0], gs_1, temp);
		ra_size( &s[D_V0], n_  );
		ra_size( &s[D_P0], n_0 );
		ra_size( &s[D_S0], n_1 );

		//computing widths of pml regions
		/* pml region I */
		rd_empty( ld_pml+0, D_V0, &empty );
		if (empty) {
			pml_I_width_y = 0;
		}
		else {
			s_pml = ld_pml[0]._s;
			pml_I_width_y = s_pml[D_V0]._dims[1].n;
		}
	
		/* pml region II */
		rd_empty( ld_pml+2, D_V0, &empty );
		if (empty) {
			pml_II_width_y = 0;
		}
		else {
			s_pml = ld_pml[2]._s;;
			pml_II_width_y = s_pml[D_V0]._dims[1].n;
		}
	
		/* pml region III */
		rd_empty( ld_pml+4, D_V0, &empty );
		if (empty) {
			pml_III_width_x = 0;
		}
		else {
			s_pml = ld_pml[4]._s;
			pml_III_width_x = s_pml[D_V0]._dims[0].n;
		}
	
		/* pml region IV */
		rd_empty( ld_pml+6, D_V0, &empty );
		if (empty) {
			pml_IV_width_x = 0;
		}
		else {
			s_pml = ld_pml[6]._s;
			pml_IV_width_x = s_pml[D_V0]._dims[0].n;
		}

		return esg_2d_4_v0 ( _vx,  _multi1, 
	         		     _px,  _sxy,
		 		     _epx, _epy,
		 		     _f1_x_I, _f1_x_II, _f1_x_III, _f1_x_IV,
		 		     _f1_y_I, _f1_y_II, _f1_y_III, _f1_y_IV,
		 	             pml_I_width_y, pml_II_width_y, pml_III_width_x, pml_IV_width_x,
		 		     n_,  n_0,  n_1,
		 		     gs_, gs_0, gs_1,
		 		     la_x, la_y, dt     );
	}
	/*~~~ updating velocity vy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  	if ( iarr == D_V1 ){
		/* Initializing fields */
		_vy     = s[D_V1]._s2;
		_sxy    = s[D_S0]._s2;
		_py     = s[D_P1]._s2;

		_multi1 = s[D_MV1]._s2;

		_epx    = s[D_EP[0]]._s1;	/* 1D */
		_epy    = s[D_EV[1]]._s1;  	/* 1D */

		/* pml fields */
		s_pml   = ld_pml[0]._s;
		_f1_x_I = s_pml[D_V1]._s2; 
		s_pml   = ld_pml[1]._s;
		_f1_y_I = s_pml[D_V1]._s2;
	
		s_pml    = ld_pml[2]._s;
		_f1_x_II = s_pml[D_V1]._s2; 
		s_pml    = ld_pml[3]._s;
		_f1_y_II = s_pml[D_V1]._s2;
	
		s_pml     = ld_pml[4]._s;
		_f1_x_III = s_pml[D_V1]._s2; 
		s_pml     = ld_pml[5]._s;
		_f1_y_III = s_pml[D_V1]._s2;
		
		s_pml = ld_pml[6]._s;
		_f1_x_IV = s_pml[D_V1]._s2;
		s_pml = ld_pml[7]._s;
		_f1_y_IV = s_pml[D_V1]._s2;

		/* getting starting indexes and size of fields */
		ra_gse( &s[D_V1], gs_ , temp);
		ra_gse( &s[D_S0], gs_0, temp);
		ra_gse( &s[D_P1], gs_1, temp);
		ra_size( &s[D_V1], n_  );
		ra_size( &s[D_S0], n_0 );
		ra_size( &s[D_P1], n_1 );

		//computing widths of pml regions
		/* pml region I */
		rd_empty( ld_pml+0, D_V1, &empty );
		if (empty) {
			pml_I_width_y = 0;
		}
		else {
			s_pml = ld_pml[0]._s;
			pml_I_width_y = s_pml[D_V1]._dims[1].n;
		}
	
		/* pml region II */
		rd_empty( ld_pml+2, D_V1, &empty );
		if (empty) {
			pml_II_width_y = 0;
		}
		else {
			s_pml = ld_pml[2]._s;;
			pml_II_width_y = s_pml[D_V1]._dims[1].n;
		}
	
		/* pml region III */
		rd_empty( ld_pml+4, D_V1, &empty );
		if (empty) {
			pml_III_width_x = 0;
		}
		else {
			s_pml = ld_pml[4]._s;
			pml_III_width_x = s_pml[D_V1]._dims[0].n;
		}
	
		/* pml region IV */
		rd_empty( ld_pml+6, D_V1, &empty );
		if (empty) {
			pml_IV_width_x = 0;
		}
		else {
			s_pml = ld_pml[6]._s;
			pml_IV_width_x = s_pml[D_V1]._dims[0].n;
		}

		return esg_2d_4_v1 ( _vy,  _multi1, 
	         		     _sxy, _py,
		 		     _epx, _epy,
		 		     _f1_x_I, _f1_x_II, _f1_x_III, _f1_x_IV,
		 		     _f1_y_I, _f1_y_II, _f1_y_III, _f1_y_IV,
		 	             pml_I_width_y, pml_II_width_y, pml_III_width_x, pml_IV_width_x,
		 		     n_,  n_0,  n_1,
		 		     gs_, gs_0, gs_1,
		 		     la_x, la_y, dt     );
	}
	return E_NOTIMESTEP;
}

//Definition of auxilary functions
/*----------------------------------------------------------------------------*/
int esg_2d_4_ns( float ** _px,   float ** _py, 
		 float ** _mp00, float ** _mp01, 
	         float ** _vx,   float ** _vy,
		 float * _epx,   float * _epy,
		 float ** _px_x_I, float ** _px_x_II, float ** _px_x_III, float ** _px_x_IV,
		 float ** _px_y_I, float ** _px_y_II, float ** _px_y_III, float ** _px_y_IV,
		 float ** _py_x_I, float ** _py_x_II, float ** _py_x_III, float ** _py_x_IV,
		 float ** _py_y_I, float ** _py_y_II, float ** _py_y_III, float ** _py_y_IV,
		 int pml_I_width_y, int pml_II_width_y, int pml_III_width_x, int pml_IV_width_x,
		 int hlo_I_width_y, int hlo_II_width_y, int hlo_III_width_x, int hlo_IV_width_x,
		 int * n_,  int * n_V0,  int * n_V1,
		 int * gs_, int * gs_V0, int * gs_V1,
		 float la_x, float la_y, float dt ){
/*----------------------------------------------------------------------------*/
#ifdef MARIO_VERBOSE
	fprintf(stderr,"> Inside esg_2d_4_ns\n");
#endif

	if ( n_[0] * n_[1] == 0 ) return 0;

	int tsz, tid;

	int width[] = {0,0};

	int offset_[]   = {0,0}; 
	int offset_V0[] = {0,0};
	int offset_V1[] = {0,0};

	int offset_III[]    = {0,0};
	int offset_V0_III[] = {0,0};
	int offset_V1_III[] = {0,0};

	int offset_IV[]    = {0,0};
	int offset_V0_IV[] = {0,0};
	int offset_V1_IV[] = {0,0};

	int index[] = {0,0};
	int i_[]    = {0,0};
	int i_V0[]  = {0,0};
	int i_V1[]  = {0,0};

	register float vx3, vx2, vx1, vx0;
	register float vy3, vy2, vy1, vy0; 
	register float dfdx, dfdy, etaxdt, etaydt;
	register float aux;
	register float dt2 = dt / 2.0;
	
#ifdef SUPER_MARIO_VERBOSE
	fprintf(stderr,"gs_  [0] = %d, gs_  [1] = %d\n",gs_[0],gs_[1]);
	fprintf(stderr,"gs_V0[0] = %d, gs_V0[1] = %d\n",gs_V0[0],gs_V0[1]);
	fprintf(stderr,"gs_V1[0] = %d, gs_V1[1] = %d\n\n",gs_V1[0],gs_V1[1]);

	fprintf(stderr,"n_  [0] = %d, n_  [1] = %d\n",n_[0],n_[1]);
	fprintf(stderr,"n_V0[0] = %d, n_V0[1] = %d\n",n_V0[0],n_V0[1]);
	fprintf(stderr,"n_V1[0] = %d, n_V1[1] = %d\n\n",n_V1[0],n_V1[1]);

	fprintf(stderr,"pml_I_width_y   = %d\n",pml_I_width_y);
	fprintf(stderr,"pml_II_width_y  = %d\n",pml_II_width_y);
	fprintf(stderr,"pml_III_width_x = %d\n",pml_III_width_x);
	fprintf(stderr,"pml_IV_width_x  = %d\n\n",pml_IV_width_x);
#endif

#pragma omp parallel private(tsz,tid,index,i_,i_V0,i_V1,offset_,offset_V0_III,offset_V1_III,offset_V0_IV,offset_V1_IV,_vx,_vy,_mp00,_mp01,_epx,_epy,_vx,_vy,vx3,vx2,vx1,vx0,vy3,vy2,vy1,vy0,_px_x_I,_px_y_I,_py_x_I,_py_y_I,_px_x_II,_px_y_II,_py_x_II,_py_y_II,_px_x_III,_px_y_III,_py_x_III,_py_y_III,_px_x_IV,_px_y_IV,_py_x_IV,_py_y_IV,dfdx,dfdy,etaxdt,etaydt)
	{
#ifdef _OPENMP
	tsz = omp_get_num_threads();
	tid = omp_get_thread_num();
#else
	tsz = 1;
	tid = 0;
#endif
	
#pragma omp single

#pragma omp barrier


#ifdef SUPER_MARIO_VERBOSE
	fprintf(stderr,"going into pml\n");
#endif
	//##################################################################################//
	//# PML + PHYSICAL REGIONS #########################################################//
	//##################################################################################//

	/* pml region I *********************************************************************/
	if ( pml_I_width_y>1 ){

		//computing width of loops and offsets
		width[0]  = n_[0]-1;
		width[1]  = pml_I_width_y;
		
		offset_[0]   = gs_[0];
		offset_V0[0] = gs_V0[0];
		offset_V1[0] = gs_V0[0];
		
		offset_[1]   = gs_[1] + tid;
		offset_V0[1] = gs_V0[1] + tid;
		offset_V1[1] = gs_V1[1] + tid;

		//loop in y-direction
		for ( index[1] = 1+tid; index[1] < width[1]; index[1]++ ){
			i_[1]   = index[1] + offset_[1];
			i_V0[1] = index[1] + offset_V0[1];
			i_V1[1] = index[1] + offset_V1[1];

			etaydt = _epy[i_[1]] * dt2;

			//loop in x-direction
			for ( index[0] = 1; index[0] < width[0]; index[0]++ ){
				i_[0]   = index[0] + offset_[0];
				i_V0[0] = index[0] + offset_V0[0];
				i_V1[0] = index[0] + offset_V1[0];

				etaxdt = _epx[i_[0]] * dt2;
				
				vx0 = _vx[ i_V0[1] ][ i_V0[0]-2 ];
				vx1 = _vx[ i_V0[1] ][ i_V0[0]-1 ];
				vx2 = _vx[ i_V0[1] ][ i_V0[0]   ];
				vx3 = _vx[ i_V0[1] ][ i_V0[0]+1 ];

				vy0 = _vy[ i_V1[1]-2 ][ i_V1[0] ];
				vy1 = _vy[ i_V1[1]-1 ][ i_V1[0] ];
				vy2 = _vy[ i_V1[1]   ][ i_V1[0] ];
				vy3 = _vy[ i_V1[1]+1 ][ i_V1[0] ];
				
				dfdx = ( (vx0 - vx3) + (vx2 - vx1) * 27.0 ) * la_x; 
				dfdy = ( (vy0 - vy3) + (vy2 - vy1) * 27.0 ) * la_y;

				_px_x_I[ i_[1] ][ i_[0] ] = ( _px_x_I[ i_[1] ][ i_[0] ]*(1.0f - etaxdt) 
                                                                                      + dfdx*_mp00[ i_[1] ][ i_[0] ] ) / (1.0f + etaxdt);
				_py_x_I[ i_[1] ][ i_[0] ] = ( _py_x_I[ i_[1] ][ i_[0] ]*(1.0f - etaxdt) 
                                                                                      + dfdx*_mp01[ i_[1] ][ i_[0] ] ) / (1.0f + etaxdt);
				_px_y_I[ i_[1] ][ i_[0] ] = ( _px_y_I[ i_[1] ][ i_[0] ]*(1.0f - etaydt) 
                                                                                      + dfdy*_mp01[ i_[1] ][ i_[0] ] ) / (1.0f + etaydt);
				_py_y_I[ i_[1] ][ i_[0] ] = ( _py_y_I[ i_[1] ][ i_[0] ]*(1.0f - etaydt) 
                                                                                      + dfdy*_mp00[ i_[1] ][ i_[0] ] ) / (1.0f + etaydt);

				_px[ i_[1] ][ i_[0] ] = _px_x_I[ i_[1] ][ i_[0] ] + _px_y_I[ i_[1] ][ i_[0] ];
				_py[ i_[1] ][ i_[0] ] = _py_x_I[ i_[1] ][ i_[0] ] + _py_y_I[ i_[1] ][ i_[0] ];
			}
		}
	}

#ifdef SUPER_MARIO_VERBOSE
	fprintf(stderr,"out of pml_I loop\n");
#endif

	/* pml region III, IV and physical region *******************************************/
	
	//computing width of loops and offsets
	width[0] = n_[0] - pml_III_width_x - pml_IV_width_x;
	width[1] = n_[1] - pml_I_width_y   - pml_II_width_y;

	offset_[0]       = pml_III_width_x + gs_[0];
	offset_V0[0]     = pml_III_width_x + gs_V0[0];
	offset_V1[0]     = pml_III_width_x + gs_V1[0];

	offset_[1]       = pml_I_width_y + gs_[1]  + tid;
	offset_V0[1]     = pml_I_width_y + gs_V0[1] + tid;
	offset_V1[1]     = pml_I_width_y + gs_V1[1] + tid;

	offset_III[0]    = gs_[0];
	offset_V0_III[0] = gs_V0[0];
	offset_V1_III[0] = gs_V1[0];

	offset_IV[0]     = n_[0] - pml_IV_width_x + gs_[0];
	offset_V0_IV[0]  = n_[0] - pml_IV_width_x + gs_V0[0];
	offset_V1_IV[0]  = n_[0] - pml_IV_width_x + gs_V1[0];

	//loop in y-direction
	for ( index[1] = tid; index[1] < width[1]; index[1]++ ){
		i_[1]   = index[1] + offset_[1];
		i_V0[1] = index[1] + offset_V0[1];
		i_V1[1] = index[1] + offset_V1[1];

		etaydt = _epy[i_[1]] * dt2;

		/* pml region III */
		for ( index[0] = 1; index[0] < pml_III_width_x; index[0]++ ){
			i_[0]   = index[0] + offset_III[0];
			i_V0[0] = index[0] + offset_V0_III[0];
			i_V1[0] = index[0] + offset_V1_III[0];

			etaxdt = _epx[i_[0]] * dt2;
			
			vx0 = _vx[ i_V0[1] ][ i_V0[0]-2 ];
			vx1 = _vx[ i_V0[1] ][ i_V0[0]-1 ];
			vx2 = _vx[ i_V0[1] ][ i_V0[0]   ];
			vx3 = _vx[ i_V0[1] ][ i_V0[0]+1 ];

			vy0 = _vy[ i_V1[1]-2 ][ i_V1[0] ];
			vy1 = _vy[ i_V1[1]-1 ][ i_V1[0] ];
			vy2 = _vy[ i_V1[1]   ][ i_V1[0] ];
			vy3 = _vy[ i_V1[1]+1 ][ i_V1[0] ];

			dfdx = ( (vx0 - vx3) + (vx2 - vx1) * 27.0 ) * la_x; 
			dfdy = ( (vy0 - vy3) + (vy2 - vy1) * 27.0 ) * la_y;
			
			_px_x_III[ i_[1] ][ i_[0] ] = ( _px_x_III[ i_[1] ][ i_[0] ]*(1.0f - etaxdt) 
										  + dfdx*_mp00[ i_[1] ][ i_[0] ] ) / (1.0f + etaxdt);
			_py_x_III[ i_[1] ][ i_[0] ] = ( _py_x_III[ i_[1] ][ i_[0] ]*(1.0f - etaxdt) 
							                          + dfdx*_mp01[ i_[1] ][ i_[0] ] ) / (1.0f + etaxdt);
			_px_y_III[ i_[1] ][ i_[0] ] =   _px_y_III[ i_[1] ][ i_[0] ] + dfdy*_mp01[ i_[1] ][ i_[0] ];
			_py_y_III[ i_[1] ][ i_[0] ] =   _py_y_III[ i_[1] ][ i_[0] ] + dfdy*_mp00[ i_[1] ][ i_[0] ];
		
			_px[ i_[1] ][ i_[0] ] = _px_x_III[ i_[1] ][ i_[0] ] + _px_y_III[ i_[1] ][ i_[0] ];
			_py[ i_[1] ][ i_[0] ] = _py_x_III[ i_[1] ][ i_[0] ] + _py_y_III[ i_[1] ][ i_[0] ];
		}
	
		/* physical region */
		for ( index[0] = 0; index[0] < width[0]; index[0]++ ){
			i_[0]   = index[0] + offset_[0];
			i_V0[0] = index[0] + offset_V0[0];
			i_V1[0] = index[0] + offset_V1[0];

			vx0 = _vx[ i_V0[1] ][ i_V0[0]-2 ];
			vx1 = _vx[ i_V0[1] ][ i_V0[0]-1 ];
			vx2 = _vx[ i_V0[1] ][ i_V0[0]   ];
			vx3 = _vx[ i_V0[1] ][ i_V0[0]+1 ];

			vy0 = _vy[ i_V1[1]-2 ][ i_V1[0] ];
			vy1 = _vy[ i_V1[1]-1 ][ i_V1[0] ];
			vy2 = _vy[ i_V1[1]   ][ i_V1[0] ];
			vy3 = _vy[ i_V1[1]+1 ][ i_V1[0] ];
			
			dfdx = ( (vx0 - vx3) + (vx2 - vx1) * 27.0 ) * la_x; 
			dfdy = ( (vy0 - vy3) + (vy2 - vy1) * 27.0 ) * la_y;
			
			_px[ i_[1] ][ i_[0] ] = _px[ i_[1] ][ i_[0] ] + dfdx*_mp00[ i_[1] ][ i_[0] ] + dfdy*_mp01[ i_[1] ][ i_[0] ];
			_py[ i_[1] ][ i_[0] ] = _py[ i_[1] ][ i_[0] ] + dfdx*_mp01[ i_[1] ][ i_[0] ] + dfdy*_mp00[ i_[1] ][ i_[0] ];
		}
	
		/* pml region IV */
		for ( index[0] = 0; index[0] < pml_IV_width_x-1; index[0]++ ){
			i_[0]   = index[0] + offset_IV[0];
			i_V0[0] = index[0] + offset_V0_IV[0];
			i_V1[0] = index[0] + offset_V1_IV[0];

			etaxdt = _epx[i_[0]] * dt2;

			vx0 = _vx[ i_V0[1] ][ i_V0[0]-2 ];
			vx1 = _vx[ i_V0[1] ][ i_V0[0]-1 ];
			vx2 = _vx[ i_V0[1] ][ i_V0[0]   ];
			vx3 = _vx[ i_V0[1] ][ i_V0[0]+1 ];

			vy0 = _vy[ i_V1[1]-2 ][ i_V1[0] ];
			vy1 = _vy[ i_V1[1]-1 ][ i_V1[0] ];
			vy2 = _vy[ i_V1[1]   ][ i_V1[0] ];
			vy3 = _vy[ i_V1[1]+1 ][ i_V1[0] ];

			dfdx = ( (vx0 - vx3) + (vx2 - vx1) * 27.0 ) * la_x; 
			dfdy = ( (vy0 - vy3) + (vy2 - vy1) * 27.0 ) * la_y;
			
			_px_x_IV[ i_[1] ][ i_[0] ] = ( _px_x_IV[ i_[1] ][ i_[0] ] * (1.0f - etaxdt) 
										+ dfdx*_mp00[ i_[1] ][ i_[0] ] ) / (1.0f + etaxdt);
			_py_x_IV[ i_[1] ][ i_[0] ] = ( _py_x_IV[ i_[1] ][ i_[0] ] * (1.0f - etaxdt) 
										+ dfdx*_mp01[ i_[1] ][ i_[0] ] ) / (1.0f + etaxdt);
			_px_y_IV[ i_[1] ][ i_[0] ] =   _px_y_IV[ i_[1] ][ i_[0] ] + dfdy*_mp01[ i_[1] ][ i_[0] ];
			_py_y_IV[ i_[1] ][ i_[0] ] =   _py_y_IV[ i_[1] ][ i_[0] ] + dfdy*_mp00[ i_[1] ][ i_[0] ];
		
			_px[ i_[1] ][ i_[0] ] = _px_x_IV[ i_[1] ][ i_[0] ] + _px_y_IV[ i_[1] ][ i_[0] ];
			_py[ i_[1] ][ i_[0] ] = _py_x_IV[ i_[1] ][ i_[0] ] + _py_y_IV[ i_[1] ][ i_[0] ];
		}
	}

#ifdef SUPER_MARIO_VERBOSE
	fprintf(stderr,"out of pml_III + pml_IV + physical loops\n");
#endif

    	/* pml region II ********************************************************************/
	if ( pml_II_width_y>1 ){

		//computing loop widths and offsets
		width[0] = n_[0]-1;
		width[1] = pml_II_width_y-1;

		offset_[0]   = gs_[0];
		offset_V0[0] = gs_V0[0];
		offset_V1[0] = gs_V1[0];

		offset_[1]   = n_[1] - pml_II_width_y + gs_[1]   + tid;
		offset_V0[1] = n_[1] - pml_II_width_y + gs_V0[1] + tid;
		offset_V1[1] = n_[1] - pml_II_width_y + gs_V1[1] + tid;

		//main loop in y-direction
		for ( index[1] = 0+tid; index[1] < width[1]; index[1]++ ){
			i_[1]   = index[1] + offset_[1];
			i_V0[1] = index[1] + offset_V0[1];
			i_V1[1] = index[1] + offset_V1[1];

			etaydt = _epy[i_[1]] * dt2;
	
			//loop in x-direction
			for ( index[0] = 1; index[0] < width[0]; index[0]++ ){
				i_[0]   = index[0] + offset_[0];
				i_V0[0] = index[0] + offset_V0[0];
				i_V1[0] = index[0] + offset_V1[0];

				etaxdt = _epx[i_[0]] * dt2;
				
				vx0 = _vx[ i_V0[1] ][ i_V0[0]-2 ];
				vx1 = _vx[ i_V0[1] ][ i_V0[0]-1 ];
				vx2 = _vx[ i_V0[1] ][ i_V0[0]   ];
				vx3 = _vx[ i_V0[1] ][ i_V0[0]+1 ];

				vy0 = _vy[ i_V1[1]-2 ][ i_V1[0] ];
				vy1 = _vy[ i_V1[1]-1 ][ i_V1[0] ];
				vy2 = _vy[ i_V1[1]   ][ i_V1[0] ];
				vy3 = _vy[ i_V1[1]+1 ][ i_V1[0] ];
				
				dfdx = ( (vx0 - vx3) + (vx2 - vx1) * 27.0 ) * la_x; 
				dfdy = ( (vy0 - vy3) + (vy2 - vy1) * 27.0 ) * la_y;

				_px_x_II[ i_[1] ][ i_[0] ] = ( _px_x_II[ i_[1] ][ i_[0] ]*(1.0f - etaxdt) 
                                                                                      + dfdx*_mp00[ i_[1] ][ i_[0] ] ) / (1.0f + etaxdt);
				_py_x_II[ i_[1] ][ i_[0] ] = ( _py_x_II[ i_[1] ][ i_[0] ]*(1.0f - etaxdt) 
                                                                                      + dfdx*_mp01[ i_[1] ][ i_[0] ] ) / (1.0f + etaxdt);
				_px_y_II[ i_[1] ][ i_[0] ] = ( _px_y_II[ i_[1] ][ i_[0] ]*(1.0f - etaydt) 
                                                                                      + dfdy*_mp01[ i_[1] ][ i_[0] ] ) / (1.0f + etaydt);
				_py_y_II[ i_[1] ][ i_[0] ] = ( _py_y_II[ i_[1] ][ i_[0] ]*(1.0f - etaydt) 
                                                                                      + dfdy*_mp00[ i_[1] ][ i_[0] ] ) / (1.0f + etaydt);

				_px[ i_[1] ][ i_[0] ] = _px_x_II[ i_[1] ][ i_[0] ] + _px_y_II[ i_[1] ][ i_[0] ];
				_py[ i_[1] ][ i_[0] ] = _py_x_II[ i_[1] ][ i_[0] ] + _py_y_II[ i_[1] ][ i_[0] ];
			}
		}
	}

#ifdef SUPER_MARIO_VERBOSE
	fprintf(stderr,"out of pml_II loops\n");
#endif

	//##################################################################################//
	//# BOUNDARY REGIONS ###############################################################//
	//##################################################################################//

    	/* boundary region I ***************************************************************/
	//computing loop widths and offsets
	width[0]  = n_[0] - 1; //excluding corners!

	offset_[0]   = gs_[0];
	offset_V0[0] = gs_V0[0];
	offset_V1[0] = gs_V1[0];

	offset_[1]   = gs_[1];
	offset_V0[1] = gs_V0[1];
	offset_V1[1] = gs_V1[1];

	i_[1]   = 0 + offset_[1]; 
	i_V0[1] = 0 + offset_V0[1]; 
	i_V1[1] = 0 + offset_V1[1]; 

	//case where boundary is in pml region
	if ( pml_I_width_y > 0 ){
		//loop in x-direction
		for (index[0] = 1; index[0] < width[0]; index[0]++){
			i_[0]   = index[0] + offset_[0]; 
			i_V0[0] = index[0] + offset_V0[0]; 
			i_V1[0] = index[0] + offset_V1[0]; 

			etaxdt = _epx[i_[0]] * dt2;

			/* Only updating _px field since _py = 0. Note that x-derivative is only used. */
			vx0 = _vx[ i_V0[1] ][ i_V0[0]-2 ];
			vx1 = _vx[ i_V0[1] ][ i_V0[0]-1 ];
			vx2 = _vx[ i_V0[1] ][ i_V0[0]   ];
			vx3 = _vx[ i_V0[1] ][ i_V0[0]+1 ];

			dfdx = ( (vx0 - vx3) + (vx2 - vx1) * 27.0 ) * la_x; 
			// aux = lambda + 2*mu + lambda^2/(lambda+2*mu)
			aux  = _mp00[i_[1]][i_[0]] + _mp01[i_[1]][i_[0]]*_mp01[i_[1]][i_[0]] / _mp00[i_[1]][i_[0]]; 

			_px[ i_[1] ][ i_[0] ] = ( _px[ i_[1] ][ i_[0] ]*(1.0f - etaxdt) + dfdx*aux ) / (1.0f + etaxdt);	
		}
	}
	//case where boundary is in physical region
	else if ( pml_I_width_y == 0 ){
		//loop in x-direction
		for (index[0] = 1; index[0] < width[0]; index[0]++){
			i_[0]   = index[0] + offset_[0]; 
			i_V0[0] = index[0] + offset_V0[0]; 
			i_V1[0] = index[0] + offset_V1[0];

			/* Only updating _px field since _py = 0. Note that x-derivative is only used. */
			vx0 = _vx[ i_V0[1] ][ i_V0[0]-2 ];
			vx1 = _vx[ i_V0[1] ][ i_V0[0]-1 ];
			vx2 = _vx[ i_V0[1] ][ i_V0[0]   ];
			vx3 = _vx[ i_V0[1] ][ i_V0[0]+1 ];

			dfdx = ( (vx0 - vx3) + (vx2 - vx1) * 27.0 ) * la_x; 
			// aux = lambda + 2*mu + lambda^2/(lambda+2*mu)
			aux  = _mp00[i_[1]][i_[0]] + _mp01[i_[1]][i_[0]]*_mp01[i_[1]][i_[0]] / _mp00[i_[1]][i_[0]]; 

			_px[ i_[1] ][ i_[0] ] = _px[ i_[1] ][ i_[0] ] + dfdx*aux;	
		}
	}

#ifdef SUPER_MARIO_VERBOSE
	fprintf(stderr,"out of boundary I loops\n");
#endif

    	/* boundary region II **************************************************************/
	//computing loop widths and offsets
	width[0]  = n_[0] - 1;//excluding corners!

	offset_[0]   = gs_[0];
	offset_V0[0] = gs_V0[0];
	offset_V1[0] = gs_V1[0];

	offset_[1]   = gs_[1];
	offset_V0[1] = gs_V0[1];
	offset_V1[1] = gs_V1[1];

	i_[1]   = n_[1]-1 + offset_[1];
	i_V0[1] = n_[1]-1 + offset_V0[1];
	i_V1[1] = n_[1]-1 + offset_V1[1];

	//case where boundary is in pml region
	if ( pml_II_width_y > 0 ){
		//loop in x-direction
		for (index[0] = 1; index[0] < width[0]; index[0]++){
			i_[0]   = index[0] + offset_[0]; 
			i_V0[0] = index[0] + offset_V0[0]; 
			i_V1[0] = index[0] + offset_V1[0]; 

			etaxdt = _epx[i_[0]] * dt2;

			/* Only updating _px field since _py = 0. Note that x-derivative is only used. */
			vx0 = _vx[ i_V0[1] ][ i_V0[0]-2 ];
			vx1 = _vx[ i_V0[1] ][ i_V0[0]-1 ];
			vx2 = _vx[ i_V0[1] ][ i_V0[0]   ];
			vx3 = _vx[ i_V0[1] ][ i_V0[0]+1 ];

			dfdx = ( (vx0 - vx3) + (vx2 - vx1) * 27.0 ) * la_x; 
			// aux = lambda + 2*mu + lambda^2/(lambda+2*mu)
			aux  = _mp00[i_[1]][i_[0]] + _mp01[i_[1]][i_[0]]*_mp01[i_[1]][i_[0]] / _mp00[i_[1]][i_[0]]; 

			_px[ i_[1] ][ i_[0] ] = ( _px[ i_[1] ][ i_[0] ]*(1.0f - etaxdt) + dfdx*aux ) / (1.0f + etaxdt);	
		}
	}
	//case where boundary is in physical region
	else if ( pml_II_width_y == 0 ){
		//loop in x-direction
		for (index[0] = 1; index[0] < width[0]; index[0]++){
			i_[0]   = index[0] + offset_[0]; 
			i_V0[0] = index[0] + offset_V0[0]; 
			i_V1[0] = index[0] + offset_V1[0];

			/* Only updating _px field since _py = 0. Note that x-derivative is only used. */
			vx0 = _vx[ i_V0[1] ][ i_V0[0]-2 ];
			vx1 = _vx[ i_V0[1] ][ i_V0[0]-1 ];
			vx2 = _vx[ i_V0[1] ][ i_V0[0]   ];
			vx3 = _vx[ i_V0[1] ][ i_V0[0]+1 ];

			dfdx = ( (vx0 - vx3) + (vx2 - vx1) * 27.0 ) * la_x; 
			// aux = lambda + 2*mu + lambda^2/(lambda+2*mu)
			aux  = _mp00[i_[1]][i_[0]] + _mp01[i_[1]][i_[0]]*_mp01[i_[1]][i_[0]] / _mp00[i_[1]][i_[0]]; 

			_px[ i_[1] ][ i_[0] ] = _px[ i_[1] ][ i_[0] ] + dfdx*aux;	
		}
	}

#ifdef SUPER_MARIO_VERBOSE
	fprintf(stderr,"out of boundary II loops\n");
#endif

    	/* boundary region III **************************************************************/
	//computing loop widths and offsets
	width[1]  = n_[1] - 1;//excluding corners!

	offset_[0]   = gs_[0];
	offset_V0[0] = gs_V0[0];
	offset_V1[0] = gs_V1[0];

	offset_[1]   = gs_[1];
	offset_V0[1] = gs_V0[1];
	offset_V1[1] = gs_V1[1];

	i_[0]    = 0 + offset_[0];
	i_V0[0]  = 0 + offset_V0[0];
	i_V1[0]  = 0 + offset_V1[0];


	//case where boundary is in pml region
	if ( pml_III_width_x > 0 ){
		//loop in y-direction
		for (index[1] = 1; index[1] < width[1]; index[1]++){
			i_[1]   = index[1] + offset_[1];
			i_V0[1] = index[1] + offset_V0[1];
			i_V1[1] = index[1] + offset_V1[1];

			etaydt = _epy[i_[1]] * dt2;

			/* Only updating _py field since _px = 0. Note that y-derivative is only used. */
			vy0 = _vy[ i_V1[1]-2 ][ i_V1[0] ];
			vy1 = _vy[ i_V1[1]-1 ][ i_V1[0] ];
			vy2 = _vy[ i_V1[1]   ][ i_V1[0] ];
			vy3 = _vy[ i_V1[1]+1 ][ i_V1[0] ];

			dfdy = ( (vy0 - vy3) + (vy2 - vy1) * 27.0 ) * la_y; 
			// aux = lambda + 2*mu + lambda^2/(lambda+2*mu)
			aux  = _mp00[i_[1]][i_[0]] + _mp01[i_[1]][i_[0]]*_mp01[i_[1]][i_[0]] / _mp00[i_[1]][i_[0]]; 

			_py[ i_[1] ][ i_[0] ] = ( _py[ i_[1] ][ i_[0] ]*(1.0f - etaydt) + dfdy*aux ) / (1.0f + etaydt);	
		}
	}
	//case where boundary is in physical region
	else if ( pml_III_width_x == 0 ){
		//loop in y-direction
		for (index[1] = 1; index[1] < width[1]; index[1]++){
			i_[1]   = index[1] + offset_[1];
			i_V0[1] = index[1] + offset_V0[1];
			i_V1[1] = index[1] + offset_V1[1];

			/* Only updating _py field since _px = 0. Note that y-derivative is only used. */
			vy0 = _vy[ i_V1[1]-2 ][ i_V1[0] ];
			vy1 = _vy[ i_V1[1]-1 ][ i_V1[0] ];
			vy2 = _vy[ i_V1[1]   ][ i_V1[0] ];
			vy3 = _vy[ i_V1[1]+1 ][ i_V1[0] ];

			dfdy = ( (vy0 - vy3) + (vy2 - vy1) * 27.0 ) * la_y; 
			// aux = lambda + 2*mu + lambda^2/(lambda+2*mu)
			aux  = _mp00[i_[1]][i_[0]] + _mp01[i_[1]][i_[0]]*_mp01[i_[1]][i_[0]] / _mp00[i_[1]][i_[0]];

			_py[ i_[1] ][ i_[0] ] = _py[ i_[1] ][ i_[0] ] + dfdy*aux;	
		}
	}

#ifdef SUPER_MARIO_VERBOSE
	fprintf(stderr,"out of boundary III loops\n");
#endif

    	/* boundary region IV ***************************************************************/
	//computing loop widths and offsets
	width[1] = n_[1] - 1;//excluding corners!

	offset_[0]   = gs_[0];
	offset_V0[0] = gs_V0[0];
	offset_V1[0] = gs_V1[0];

	offset_[1]   = gs_[1];
	offset_V0[1] = gs_V0[1];
	offset_V1[1] = gs_V1[1];

	i_[0]    = n_[0]-1 + offset_[0];
	i_V0[0]  = n_[0]-1 + offset_V0[0];
	i_V1[0]  = n_[0]-1 + offset_V1[0];

	//case where boundary is in pml region
	if ( pml_IV_width_x > 0 ){
		//loop in y-direction
		for (index[1] = 1; index[1] < width[1]; index[1]++){
			i_[1]   = index[1] + offset_[1];
			i_V0[1] = index[1] + offset_V0[1];
			i_V1[1] = index[1] + offset_V1[1];

			etaydt = _epy[i_[1]] * dt2;

			/* Only updating _py field since _px = 0. Note that y-derivative is only used. */
			vy0 = _vy[ i_V1[1]-2 ][ i_V1[0] ];
			vy1 = _vy[ i_V1[1]-1 ][ i_V1[0] ];
			vy2 = _vy[ i_V1[1]   ][ i_V1[0] ];
			vy3 = _vy[ i_V1[1]+1 ][ i_V1[0] ];

			dfdy = ( (vy0 - vy3) + (vy2 - vy1) * 27.0 ) * la_y; 
			// aux = lambda + 2*mu + lambda^2/(lambda+2*mu)
			aux  = _mp00[i_[1]][i_[0]] + _mp01[i_[1]][i_[0]]*_mp01[i_[1]][i_[0]] / _mp00[i_[1]][i_[0]]; 

			_py[ i_[1] ][ i_[0] ] = ( _py[ i_[1] ][ i_[0] ]*(1.0f - etaydt) + dfdy*aux ) / (1.0f + etaydt);	
		}
	}
	//case where boundary is in physical region
	else if ( pml_IV_width_x == 0 ){
		//loop in y-direction
		for (index[1] = 1; index[1] < width[1]; index[1]++){
			i_[1]   = index[1] + offset_[1];
			i_V0[1] = index[1] + offset_V0[1];
			i_V1[1] = index[1] + offset_V1[1];

			/* Only updating _py field since _px = 0. Note that y-derivative is only used. */
			vy0 = _vy[ i_V1[1]-2 ][ i_V1[0] ];
			vy1 = _vy[ i_V1[1]-1 ][ i_V1[0] ];
			vy2 = _vy[ i_V1[1]   ][ i_V1[0] ];
			vy3 = _vy[ i_V1[1]+1 ][ i_V1[0] ];

			dfdy = ( (vy0 - vy3) + (vy2 - vy1) * 27.0 ) * la_y; 
			// aux = lambda + 2*mu + lambda^2/(lambda+2*mu)
			aux  = _mp00[i_[1]][i_[0]] + _mp01[i_[1]][i_[0]]*_mp01[i_[1]][i_[0]] / _mp00[i_[1]][i_[0]];

			_py[ i_[1] ][ i_[0] ] = _py[ i_[1] ][ i_[0] ] + dfdy*aux;	
		}
	}

#ifdef SUPER_MARIO_VERBOSE
	fprintf(stderr,"out of boundary IV loops\n");
#endif

	//##################################################################################//
	//# HALO REGIONS ###################################################################//
	//##################################################################################//

    	/* Halo region I *********************************************************************/
	//computing loop widths and offsets
	width[0] = n_[0];
	width[1] = hlo_I_width_y; //halo y-width

	offset_[0] = gs_[0];
	offset_[1] = gs_[1];

	//loop in y-direction
	for (index[1] = 1; index[1] <= width[1]; index[1]++){
		i_[1] = offset_[1] - index[1];

		//loop in x-direction
		for (index[0] = 0; index[0] < width[0]; index[0]++){
			i_[0] = index[0] + offset_[0];

			_py[ i_[1] ][ i_[0] ] = - _py[ offset_[1]+index[1] ][ i_[0] ]; //odd extension
		}
	}

    	/* Halo region II ********************************************************************/
	//computing loop widths and offsets
	width[0] = n_[0];
	width[1] = hlo_II_width_y; //halo y-width

	offset_[0] = gs_[0];
	offset_[1] = n_[1] - pml_II_width_y;

	//loop in y-direction
	for (index[1] = 1; index[1] <= width[1]; index[1]++){
		i_[1] = offset_[1] + index[1];

		//loop in x-direction
		for (index[0] = 0; index[0] < width[0]; index[0]++){
			i_[0] = index[0] + offset_[0];

			_py[ i_[1] ][ i_[0] ] = - _py[ offset_[1]-index[1] ][ i_[0] ]; //odd extension
		}
	}

    	/* Halo region III *******************************************************************/
	//computing loop widths and offsets
	width[0] = hlo_III_width_x; //halo x-width
	width[1] = n_[1];

	offset_[0] = gs_[0];
	offset_[1] = gs_[1];

	//loop in y-direction
	for (index[1] = 0; index[1] < width[1]; index[1]++){
		i_[1] = index[1] + offset_[1];

		//loop in x-direction
		for (index[0] = 1; index[0] <= width[0]; index[0]++){
			i_[0] = offset_[0] - index[0];

			_px[ i_[1] ][ i_[0] ] = - _px[ i_[1] ][ offset_[0]+index[0] ]; //odd extension
		}
	}

    	/* Halo region IV ********************************************************************/
	//computing loop widths and offsets
	width[0] = hlo_IV_width_x;; //halo x-width
	width[1] = n_[1];

	offset_[0] = n_[0] - pml_IV_width_x;
	offset_[1] = gs_[1];

	//loop in y-direction
	for (index[1] = 0; index[1] < width[1]; index[1]++){
		i_[1] = index[1] + offset_[1];

		//loop in x-direction
		for (index[0] = 1; index[0] <= width[0]; index[0]++){
			i_[0] = offset_[0] + index[0];

			_px[ i_[1] ][ i_[0] ] = - _px[ i_[1] ][ offset_[0]-index[0] ]; //odd extension
		}
	}

	} /* omp parallel */

  	return 0;
}

/*----------------------------------------------------------------------------*/
int esg_2d_4_ss0( float ** _sxy,  float ** _ms0, 
	          float ** _vx,   float ** _vy,
		  float * _epx,   float * _epy,
		  float ** _sxy_x_I, float ** _sxy_x_II, float ** _sxy_x_III, float ** _sxy_x_IV,
		  float ** _sxy_y_I, float ** _sxy_y_II, float ** _sxy_y_III, float ** _sxy_y_IV,
		  int pml_I_width_y, int pml_II_width_y, int pml_III_width_x, int pml_IV_width_x,
		  int hlo_I_width_y, int hlo_II_width_y, int hlo_III_width_x, int hlo_IV_width_x,
		  int * n_,  int * n_V0,  int * n_V1,
		  int * gs_, int * gs_V0, int * gs_V1,
		  float la_x, float la_y, float dt     ){
/*----------------------------------------------------------------------------*/
#ifdef MARIO_VERBOSE
	fprintf(stderr,"> esg_2d_4_ss0\n");
#endif
	if ( n_[0] * n_[1] == 0 ) return 0;

	int tsz, tid;

	int width[] = {0,0};

	int offset_[]   = {0,0}; 
	int offset_V0[] = {0,0};
	int offset_V1[] = {0,0};

	int offset_III[]    = {0,0};
	int offset_V0_III[] = {0,0};
	int offset_V1_III[] = {0,0};

	int offset_IV[]    = {0,0};
	int offset_V0_IV[] = {0,0};
	int offset_V1_IV[] = {0,0};

	int index[] = {0,0};
	int i_[]    = {0,0};
	int i_V0[]  = {0,0};
	int i_V1[]  = {0,0};
	
	register float vx3, vx2, vx1, vx0;
	register float vy3, vy2, vy1, vy0; 
	register float dfdx, dfdy, etaxdt, etaydt;
	register float aux;
	register float dt2 = dt / 2.0;

	
#pragma omp parallel private(tsz,tid,index,i_,i_V0,i_V1,offset_,offset_V0_III,offset_V1_III,offset_V0_IV,offset_V1_IV,_sxy,_ms0,_epx,_epy,_vx,_vy,vx3,vx2,vx1,vx0,vy3,vy2,vy1,vy0,_sxy_x_I,_sxy_y_I,_sxy_x_II,_sxy_y_II,_sxy_x_III,_sxy_y_III,_sxy_x_IV,_sxy_y_IV,dfdx,dfdy,etaxdt,etaydt)
	{
#ifdef _OPENMP
	tsz = omp_get_num_threads();
	tid = omp_get_thread_num();
#else
	tsz = 1;
	tid = 0;
#endif
	
#pragma omp single

#pragma omp barrier


#ifdef SUPER_MARIO_VERBOSE
	fprintf(stderr,"going into pml\n");
#endif

	//##################################################################################//
	//# PML + PHYSICAL REGIONS #########################################################//
	//##################################################################################//

	/* pml region I *********************************************************************/
	if ( pml_I_width_y>0 ){

		//computing width of loops and offsets
		width[0]  = n_[0];
		width[1]  = pml_I_width_y;
		
		offset_[0]   = gs_[0];
		offset_V0[0] = gs_V0[0];
		offset_V1[0] = gs_V0[0];
		
		offset_[1]   = gs_[1]   + tid;
		offset_V0[1] = gs_V0[1] + tid;
		offset_V1[1] = gs_V1[1] + tid;

		//loop in y-direction
		for ( index[1] = tid; index[1] < width[1]; index[1]++ ){
			i_[1]   = index[1] + offset_[1];
			i_V0[1] = index[1] + offset_V0[1];
			i_V1[1] = index[1] + offset_V1[1];

			etaydt = _epy[i_[1]] * dt2;

			//loop in x-direction
			for ( index[0] = 0; index[0] < width[0]; index[0]++ ){
				i_[0]   = index[0] + offset_[0];
				i_V0[0] = index[0] + offset_V0[0];
				i_V1[0] = index[0] + offset_V1[0];

				etaxdt = _epx[i_[0]] * dt2;
				
				vy0 = _vy[ i_V0[1] ][ i_V0[0]-1 ];
				vy1 = _vy[ i_V0[1] ][ i_V0[0]   ];
				vy2 = _vy[ i_V0[1] ][ i_V0[0]+1 ];
				vy3 = _vy[ i_V0[1] ][ i_V0[0]+2 ];

				vx0 = _vx[ i_V1[1]-1 ][ i_V1[0] ];
				vx1 = _vx[ i_V1[1]   ][ i_V1[0] ];
				vx2 = _vx[ i_V1[1]+1 ][ i_V1[0] ];
				vx3 = _vx[ i_V1[1]+2 ][ i_V1[0] ];
				
				dfdx = ( (vy0 - vy3) + (vy2 - vy1) * 27.0 ) * la_x; 
				dfdy = ( (vx0 - vx3) + (vx2 - vx1) * 27.0 ) * la_y;

				_sxy_x_I[ i_[1] ][ i_[0] ] = ( _sxy_x_I[ i_[1] ][ i_[0] ]*(1.0f - etaxdt) 
                                                                                      + dfdx*_ms0[ i_[1] ][ i_[0] ] ) / (1.0f + etaxdt);
				_sxy_y_I[ i_[1] ][ i_[0] ] = ( _sxy_y_I[ i_[1] ][ i_[0] ]*(1.0f - etaydt) 
                                                                                      + dfdy*_ms0[ i_[1] ][ i_[0] ] ) / (1.0f + etaydt);

				_sxy[ i_[1] ][ i_[0] ] = _sxy_x_I[ i_[1] ][ i_[0] ] + _sxy_y_I[ i_[1] ][ i_[0] ];
			}
		}
	}

#ifdef SUPER_MARIO_VERBOSE
	fprintf(stderr,"out of pml_I loop\n");
#endif

	/* pml region III, IV and physical region *******************************************/
	
	//computing width of loops and offsets
	width[0] = n_[0] - pml_III_width_x - pml_IV_width_x;
	width[1] = n_[1] - pml_I_width_y   - pml_II_width_y;

	offset_[0]       = pml_III_width_x + gs_[0];
	offset_V0[0]     = pml_III_width_x + gs_V0[0];
	offset_V1[0]     = pml_III_width_x + gs_V1[0];

	offset_[1]       = pml_I_width_y + gs_[1]   + tid;
	offset_V0[1]     = pml_I_width_y + gs_V0[1] + tid;
	offset_V1[1]     = pml_I_width_y + gs_V1[1] + tid;

	offset_III[0]    = gs_[0];
	offset_V0_III[0] = gs_V0[0];
	offset_V1_III[0] = gs_V1[0];

	offset_IV[0]     = n_[0] - pml_IV_width_x + gs_[0];
	offset_V0_IV[0]  = n_[0] - pml_IV_width_x + gs_V0[0];
	offset_V1_IV[0]  = n_[0] - pml_IV_width_x + gs_V1[0];

	//loop in y-direction
	for ( index[1] = tid; index[1] < width[1]; index[1]++ ){
		i_[1]   = index[1] + offset_[1];
		i_V0[1] = index[1] + offset_V0[1];
		i_V1[1] = index[1] + offset_V1[1];

		etaydt = _epy[i_[1]] * dt2;

		/* pml region III */
		for ( index[0] = 0; index[0] < pml_III_width_x; index[0]++ ){
			i_[0]   = index[0] + offset_III[0];
			i_V0[0] = index[0] + offset_V0_III[0];
			i_V1[0] = index[0] + offset_V1_III[0];

			etaxdt = _epx[i_[0]] * dt2;
			
			vy0 = _vy[ i_V0[1] ][ i_V0[0]-1 ];
			vy1 = _vy[ i_V0[1] ][ i_V0[0]   ];
			vy2 = _vy[ i_V0[1] ][ i_V0[0]+1 ];
			vy3 = _vy[ i_V0[1] ][ i_V0[0]+2 ];

			vx0 = _vx[ i_V1[1]-1 ][ i_V1[0] ];
			vx1 = _vx[ i_V1[1]   ][ i_V1[0] ];
			vx2 = _vx[ i_V1[1]+1 ][ i_V1[0] ];
			vx3 = _vx[ i_V1[1]+2 ][ i_V1[0] ];

			dfdx = ( (vy0 - vy3) + (vy2 - vy1) * 27.0 ) * la_x; 
			dfdy = ( (vx0 - vx3) + (vx2 - vx1) * 27.0 ) * la_y;
			
			_sxy_x_III[ i_[1] ][ i_[0] ] = ( _sxy_x_III[ i_[1] ][ i_[0] ]*(1.0f - etaxdt) 
										  + dfdx*_ms0[ i_[1] ][ i_[0] ] ) / (1.0f + etaxdt);
			
			_sxy_y_III[ i_[1] ][ i_[0] ] =   _sxy_y_III[ i_[1] ][ i_[0] ] + dfdy*_ms0[ i_[1] ][ i_[0] ];
		
			_sxy[ i_[1] ][ i_[0] ] = _sxy_x_III[ i_[1] ][ i_[0] ] + _sxy_y_III[ i_[1] ][ i_[0] ];
		}
	
		/* physical region */
		for ( index[0] = 0; index[0] < width[0]; index[0]++ ){
			i_[0]   = index[0] + offset_[0];
			i_V0[0] = index[0] + offset_V0[0];
			i_V1[0] = index[0] + offset_V1[0];

			vy0 = _vy[ i_V0[1] ][ i_V0[0]-1 ];
			vy1 = _vy[ i_V0[1] ][ i_V0[0]   ];
			vy2 = _vy[ i_V0[1] ][ i_V0[0]+1 ];
			vy3 = _vy[ i_V0[1] ][ i_V0[0]+2 ];

			vx0 = _vx[ i_V1[1]-1 ][ i_V1[0] ];
			vx1 = _vx[ i_V1[1]   ][ i_V1[0] ];
			vx2 = _vx[ i_V1[1]+1 ][ i_V1[0] ];
			vx3 = _vx[ i_V1[1]+2 ][ i_V1[0] ];
			
			dfdx = ( (vy0 - vy3) + (vy2 - vy1) * 27.0 ) * la_x; 
			dfdy = ( (vx0 - vx3) + (vx2 - vx1) * 27.0 ) * la_y;
			
			_sxy[ i_[1] ][ i_[0] ] = _sxy[ i_[1] ][ i_[0] ] + dfdx*_ms0[ i_[1] ][ i_[0] ] + dfdy*_ms0[ i_[1] ][ i_[0] ];
		}
	
		/* pml region IV */
		for ( index[0] = 0; index[0] < pml_IV_width_x; index[0]++ ){
			i_[0]   = index[0] + offset_IV[0];
			i_V0[0] = index[0] + offset_V0_IV[0];
			i_V1[0] = index[0] + offset_V1_IV[0];

			etaxdt = _epx[i_[0]] * dt2;

			vy0 = _vy[ i_V0[1] ][ i_V0[0]-1 ];
			vy1 = _vy[ i_V0[1] ][ i_V0[0]   ];
			vy2 = _vy[ i_V0[1] ][ i_V0[0]+1 ];
			vy3 = _vy[ i_V0[1] ][ i_V0[0]+2 ];

			vx0 = _vx[ i_V1[1]-1 ][ i_V1[0] ];
			vx1 = _vx[ i_V1[1]   ][ i_V1[0] ];
			vx2 = _vx[ i_V1[1]+1 ][ i_V1[0] ];
			vx3 = _vx[ i_V1[1]+2 ][ i_V1[0] ];

			dfdx = ( (vy0 - vy3) + (vy2 - vy1) * 27.0 ) * la_x; 
			dfdy = ( (vx0 - vx3) + (vx2 - vx1) * 27.0 ) * la_y;
			
			_sxy_x_IV[ i_[1] ][ i_[0] ] = ( _sxy_x_IV[ i_[1] ][ i_[0] ] * (1.0f - etaxdt) 
										+ dfdx*_ms0[ i_[1] ][ i_[0] ] ) / (1.0f + etaxdt);
			_sxy_y_IV[ i_[1] ][ i_[0] ] =   _sxy_y_IV[ i_[1] ][ i_[0] ] + dfdy*_ms0[ i_[1] ][ i_[0] ];
		
			_sxy[ i_[1] ][ i_[0] ] = _sxy_x_IV[ i_[1] ][ i_[0] ] + _sxy_y_IV[ i_[1] ][ i_[0] ];
		}
	}
    	/*************************************************************************************/
#ifdef SUPER_MARIO_VERBOSE
	fprintf(stderr,"out of pml_III + pml_IV + physical loops\n");
#endif

    	/* pml region II ********************************************************************/
	if ( pml_II_width_y>0 ){

		//computing loop widths and offsets
		width[0] = n_[0];
		width[1] = pml_II_width_y;

		offset_[0]   = gs_[0];
		offset_V0[0] = gs_V0[0];
		offset_V1[0] = gs_V1[0];

		offset_[1]   = n_[1] - pml_II_width_y + gs_[1]   + tid;
		offset_V0[1] = n_[1] - pml_II_width_y + gs_V0[1] + tid;
		offset_V1[1] = n_[1] - pml_II_width_y + gs_V1[1] + tid;

		//main loop in y-direction
		for ( index[1] = tid; index[1] < width[1]; index[1]++ ){
			i_[1]   = index[1] + offset_[1];
			i_V0[1] = index[1] + offset_V0[1];
			i_V1[1] = index[1] + offset_V1[1];

			etaydt = _epy[i_[1]] * dt2;
	
			//loop in x-direction
			for ( index[0] = 0; index[0] < width[0]; index[0]++ ){
				i_[0]   = index[0] + offset_[0];
				i_V0[0] = index[0] + offset_V0[0];
				i_V1[0] = index[0] + offset_V1[0];

				etaxdt = _epx[i_[0]] * dt2;
				
				vy0 = _vy[ i_V0[1] ][ i_V0[0]-1 ];
				vy1 = _vy[ i_V0[1] ][ i_V0[0]   ];
				vy2 = _vy[ i_V0[1] ][ i_V0[0]+1 ];
				vy3 = _vy[ i_V0[1] ][ i_V0[0]+2 ];

				vx0 = _vx[ i_V1[1]-1 ][ i_V1[0] ];
				vx1 = _vx[ i_V1[1]   ][ i_V1[0] ];
				vx2 = _vx[ i_V1[1]+1 ][ i_V1[0] ];
				vx3 = _vx[ i_V1[1]+2 ][ i_V1[0] ];
				
				dfdx = ( (vy0 - vy3) + (vy2 - vy1) * 27.0 ) * la_x; 
				dfdy = ( (vx0 - vx3) + (vx2 - vx1) * 27.0 ) * la_y;

				_sxy_x_II[ i_[1] ][ i_[0] ] = ( _sxy_x_II[ i_[1] ][ i_[0] ]*(1.0f - etaxdt) 
                                                                                      + dfdx*_ms0[ i_[1] ][ i_[0] ] ) / (1.0f + etaxdt);
				_sxy_y_II[ i_[1] ][ i_[0] ] = ( _sxy_y_II[ i_[1] ][ i_[0] ]*(1.0f - etaydt) 
                                                                                      + dfdy*_ms0[ i_[1] ][ i_[0] ] ) / (1.0f + etaydt);

				_sxy[ i_[1] ][ i_[0] ] = _sxy_x_II[ i_[1] ][ i_[0] ] + _sxy_y_II[ i_[1] ][ i_[0] ];
			}
		}
	}
    	/*************************************************************************************/
#ifdef SUPER_MARIO_VERBOSE
	fprintf(stderr,"out of pml_II loops\n");
#endif

	//##################################################################################//
	//# HALO REGIONS ###################################################################//
	//##################################################################################//

    	/* Halo region I *********************************************************************/
	//computing loop widths and offsets
	width[0] = n_[0];
	width[1] = hlo_I_width_y; //halo y-width

	offset_[0] = gs_[0];
	offset_[1] = gs_[1];

	//loop in y-direction
	for (index[1] = 1; index[1] <= width[1]; index[1]++){
		i_[1] = offset_[1] - index[1];

		//loop in x-direction
		for (index[0] = 0; index[0] < width[0]; index[0]++){
			i_[0] = index[0] + offset_[0];

			_sxy[ i_[1] ][ i_[0] ] = - _sxy[ offset_[1]+index[1]-1 ][ i_[0] ]; //odd extension
		}
	}

    	/* Halo region II ********************************************************************/
	//computing loop widths and offsets
	width[0] = n_[0];
	width[1] = hlo_II_width_y; //halo y-width

	offset_[0] = gs_[0];
	offset_[1] = n_[1] - pml_II_width_y - 1;

	//loop in y-direction
	for (index[1] = 1; index[1] <= width[1]; index[1]++){
		i_[1] = offset_[1] + index[1];

		//loop in x-direction
		for (index[0] = 0; index[0] < width[0]; index[0]++){
			i_[0] = index[0] + offset_[0];

			_sxy[ i_[1] ][ i_[0] ] = - _sxy[ offset_[1]-index[1]+1 ][ i_[0] ]; //odd extension
		}
	}

    	/* Halo region III *******************************************************************/
	//computing loop widths and offsets
	width[0] = hlo_III_width_x; //halo x-width
	width[1] = n_[1];

	offset_[0] = gs_[0];
	offset_[1] = gs_[1];

	//loop in y-direction
	for (index[1] = 0; index[1] < width[1]; index[1]++){
		i_[1] = index[1] + offset_[1];

		//loop in x-direction
		for (index[0] = 1; index[0] <= width[0]; index[0]++){
			i_[0] = offset_[0] - index[0];

			_sxy[ i_[1] ][ i_[0] ] = - _sxy[ i_[1] ][ offset_[0]+index[0]-1 ]; //odd extension
		}
	}

    	/* Halo region IV ********************************************************************/
	//computing loop widths and offsets
	width[0] = hlo_IV_width_x; //halo x-width
	width[1] = n_[1];

	offset_[0] = n_[0] - pml_IV_width_x - 1;
	offset_[1] = gs_[1];

	//loop in y-direction
	for (index[1] = 0; index[1] < width[1]; index[1]++){
		i_[1] = index[1] + offset_[1];

		//loop in x-direction
		for (index[0] = 1; index[0] <= width[0]; index[0]++){
			i_[0] = offset_[0] + index[0];

			_sxy[ i_[1] ][ i_[0] ] = - _sxy[ i_[1] ][ offset_[0]-index[0]+1 ]; //odd extension
		}
	}

	} /* omp parallel */

  	return 0;
}

/*----------------------------------------------------------------------------*/
int esg_2d_4_v0( float ** _vx,  float ** _mv0, 
	         float ** _px,  float ** _sxy,
		 float * _epx,  float * _epy,
		 float ** _vx_x_I, float ** _vx_x_II, float ** _vx_x_III, float ** _vx_x_IV,
		 float ** _vx_y_I, float ** _vx_y_II, float ** _vx_y_III, float ** _vx_y_IV,
		 int pml_I_width_y, int pml_II_width_y, int pml_III_width_x, int pml_IV_width_x,
		 int * n_,  int * n_P0,  int * n_S0,
		 int * gs_, int * gs_P0, int * gs_S0,
		 float la_x, float la_y, float dt     ){
/*----------------------------------------------------------------------------*/
#ifdef MARIO_VERBOSE
	fprintf(stderr,"> esg_2d_4_v0\n");
#endif

	if ( n_[0] * n_[1] == 0 ) return 0;

	int tsz, tid;

	int width[] = {0,0};

	int offset_[]   = {0,0}; 
	int offset_P0[] = {0,0};
	int offset_S0[] = {0,0};

	int offset_III[]    = {0,0};
	int offset_P0_III[] = {0,0};
	int offset_S0_III[] = {0,0};

	int offset_IV[]    = {0,0};
	int offset_P0_IV[] = {0,0};
	int offset_S0_IV[] = {0,0};

	int index[] = {0,0};
	int i_[]    = {0,0};
	int i_P0[]  = {0,0};
	int i_S0[]  = {0,0};

	register float px3,  px2,  px1,  px0;
	register float sxy3, sxy2, sxy1, sxy0; 
	register float dfdx, dfdy, etaxdt, etaydt;
	register float aux;
	register float dt2 = dt / 2.0;
		

#pragma omp parallel private(tsz,tid,index,i_,i_P0,i_S0,offset_,offset_P0_III,offset_S0_III,offset_P0_IV,offset_S0_IV,v0,_mv0,_epx,_epy,_px,_sxy,px3,px2,px1,px0,sxy3,sxy2,sxy1,sxy0,_vx_x_I,_vx_y_I,_vx_x_II,_vx_y_II,_vx_x_III,_vx_y_III,_vx_x_IV,_vx_y_IV,dfdx,dfdy,etaxdt,etaydt)
	{
#ifdef _OPENMP
	tsz = omp_get_num_threads();
	tid = omp_get_thread_num();
#else
	tsz = 1;
	tid = 0;
#endif
	
#pragma omp single

#pragma omp barrier


#ifdef SUPER_MARIO_VERBOSE
	fprintf(stderr,"going into pml\n");
#endif

	//##################################################################################//
	//# PML + PHYSICAL REGIONS #########################################################//
	//##################################################################################//

	/* pml region I *********************************************************************/
	if ( pml_I_width_y>0 ){

		//computing width of loops and offsets
		width[0]  = n_[0];
		width[1]  = pml_I_width_y;
		
		offset_[0]   = gs_ [0];
		offset_P0[0] = gs_P0[0];
		offset_S0[0] = gs_S0[0];
		
		offset_[1]   = gs_ [1] + tid;
		offset_P0[1] = gs_P0[1] + tid;
		offset_S0[1] = gs_S0[1] + tid;

		//loop in y-direction
		for ( index[1] = tid; index[1] < width[1]; index[1]++ ){
			i_[1]   = index[1] + offset_[1];
			i_P0[1] = index[1] + offset_P0[1];
			i_S0[1] = index[1] + offset_S0[1];

			etaydt = _epy[i_[1]] * dt2;

			//loop in x-direction
			for ( index[0] = 0; index[0] < width[0]; index[0]++ ){
				i_[0]   = index[0] + offset_[0];
				i_P0[0] = index[0] + offset_P0[0];
				i_S0[0] = index[0] + offset_S0[0];

				etaxdt = _epx[i_[0]] * dt2;
				
				px0 = _px[ i_P0[1] ][ i_P0[0]-1 ];
				px1 = _px[ i_P0[1] ][ i_P0[0]   ];
				px2 = _px[ i_P0[1] ][ i_P0[0]+1 ];
				px3 = _px[ i_P0[1] ][ i_P0[0]+2 ];

				sxy0 = _sxy[ i_S0[1]-2 ][ i_S0[0] ];
				sxy1 = _sxy[ i_S0[1]-1 ][ i_S0[0] ];
				sxy2 = _sxy[ i_S0[1]   ][ i_S0[0] ];
				sxy3 = _sxy[ i_S0[1]+1 ][ i_S0[0] ];
				
				dfdx = ( (px0 - px3) + (px2 - px1) * 27.0 ) * la_x; 
				dfdy = ( (sxy0 - sxy3) + (sxy2 - sxy1) * 27.0 ) * la_y;

				_vx_x_I[ i_[1] ][ i_[0] ] = ( _vx_x_I[ i_[1] ][ i_[0] ]*(1.0f - etaxdt) 
                                                                                      + dfdx*_mv0[ i_[1] ][ i_[0] ] ) / (1.0f + etaxdt);
				_vx_y_I[ i_[1] ][ i_[0] ] = ( _vx_y_I[ i_[1] ][ i_[0] ]*(1.0f - etaydt) 
                                                                                      + dfdy*_mv0[ i_[1] ][ i_[0] ] ) / (1.0f + etaydt);

				_vx[ i_[1] ][ i_[0] ] = _vx_x_I[ i_[1] ][ i_[0] ] + _vx_y_I[ i_[1] ][ i_[0] ];
			}
		}
	}

#ifdef SUPER_MARIO_VERBOSE
	fprintf(stderr,"out of pml_I loop\n");
#endif

	/* pml region III, IV and physical region *******************************************/
	
	//computing width of loops and offsets
	width[0] = n_[0] - pml_III_width_x - pml_IV_width_x;
	width[1] = n_[1] - pml_I_width_y   - pml_II_width_y;

	offset_[0]       = pml_III_width_x + gs_[0];
	offset_P0[0]     = pml_III_width_x + gs_P0[0];
	offset_S0[0]     = pml_III_width_x + gs_S0[0];

	offset_[1]       = pml_I_width_y + gs_[1]   + tid;
	offset_P0[1]     = pml_I_width_y + gs_P0[1] + tid;
	offset_S0[1]     = pml_I_width_y + gs_S0[1] + tid;

	offset_III[0]    = gs_ [0];
	offset_P0_III[0] = gs_P0[0];
	offset_S0_III[0] = gs_S0[0];

	offset_IV[0]     = n_[0] - pml_IV_width_x + gs_[0];
	offset_P0_IV[0]  = n_[0] - pml_IV_width_x + gs_P0[0];
	offset_S0_IV[0]  = n_[0] - pml_IV_width_x + gs_S0[0];

	//loop in y-direction
	for ( index[1] = tid; index[1] < width[1]; index[1]++ ){
		i_[1]   = index[1] + offset_[1];
		i_P0[1] = index[1] + offset_P0[1];
		i_S0[1] = index[1] + offset_S0[1];

		etaydt = _epy[i_[1]] * dt2;

		/* pml region III */
		for ( index[0] = 0; index[0] < pml_III_width_x; index[0]++ ){
			i_[0]   = index[0] + offset_III[0];
			i_P0[0] = index[0] + offset_P0_III[0];
			i_S0[0] = index[0] + offset_S0_III[0];

			etaxdt = _epx[i_[0]] * dt2;
			
			px0 = _px[ i_P0[1] ][ i_P0[0]-1 ];
			px1 = _px[ i_P0[1] ][ i_P0[0]   ];
			px2 = _px[ i_P0[1] ][ i_P0[0]+1 ];
			px3 = _px[ i_P0[1] ][ i_P0[0]+2 ];

			sxy0 = _sxy[ i_S0[1]-2 ][ i_S0[0] ];
			sxy1 = _sxy[ i_S0[1]-1 ][ i_S0[0] ];
			sxy2 = _sxy[ i_S0[1]   ][ i_S0[0] ];
			sxy3 = _sxy[ i_S0[1]+1 ][ i_S0[0] ];

			dfdx = ( (px0 - px3) + (px2 - px1) * 27.0 ) * la_x; 
			dfdy = ( (sxy0 - sxy3) + (sxy2 - sxy1) * 27.0 ) * la_y;
			
			_vx_x_III[ i_[1] ][ i_[0] ] = ( _vx_x_III[ i_[1] ][ i_[0] ]*(1.0f - etaxdt) 
										  + dfdx*_mv0[ i_[1] ][ i_[0] ] ) / (1.0f + etaxdt);
			
			_vx_y_III[ i_[1] ][ i_[0] ] =   _vx_y_III[ i_[1] ][ i_[0] ] + dfdy*_mv0[ i_[1] ][ i_[0] ];
		
			_vx[ i_[1] ][ i_[0] ] = _vx_x_III[ i_[1] ][ i_[0] ] + _vx_y_III[ i_[1] ][ i_[0] ];
		}
	
		/* physical region */
		for ( index[0] = 0; index[0] < width[0]; index[0]++ ){
			i_[0]   = index[0] + offset_[0];
			i_P0[0] = index[0] + offset_P0[0];
			i_S0[0] = index[0] + offset_S0[0];

			px0 = _px[ i_P0[1] ][ i_P0[0]-1 ];
			px1 = _px[ i_P0[1] ][ i_P0[0]   ];
			px2 = _px[ i_P0[1] ][ i_P0[0]+1 ];
			px3 = _px[ i_P0[1] ][ i_P0[0]+2 ];

			sxy0 = _sxy[ i_S0[1]-2 ][ i_S0[0] ];
			sxy1 = _sxy[ i_S0[1]-1 ][ i_S0[0] ];
			sxy2 = _sxy[ i_S0[1]   ][ i_S0[0] ];
			sxy3 = _sxy[ i_S0[1]+1 ][ i_S0[0] ];
			
			dfdx = ( (px0 - px3) + (px2 - px1) * 27.0 ) * la_x; 
			dfdy = ( (sxy0 - sxy3) + (sxy2 - sxy1) * 27.0 ) * la_y;
			
			_vx[ i_[1] ][ i_[0] ] = _vx[ i_[1] ][ i_[0] ] + dfdx*_mv0[ i_[1] ][ i_[0] ] + dfdy*_mv0[ i_[1] ][ i_[0] ];
		}
	
		/* pml region IV */
		for ( index[0] = 0; index[0] < pml_IV_width_x; index[0]++ ){
			i_[0]   = index[0] + offset_IV[0];
			i_P0[0] = index[0] + offset_P0_IV[0];
			i_S0[0] = index[0] + offset_S0_IV[0];

			etaxdt = _epx[i_[0]] * dt2;

			px0 = _px[ i_P0[1] ][ i_P0[0]-1 ];
			px1 = _px[ i_P0[1] ][ i_P0[0]   ];
			px2 = _px[ i_P0[1] ][ i_P0[0]+1 ];
			px3 = _px[ i_P0[1] ][ i_P0[0]+2 ];

			sxy0 = _sxy[ i_S0[1]-2 ][ i_S0[0] ];
			sxy1 = _sxy[ i_S0[1]-1 ][ i_S0[0] ];
			sxy2 = _sxy[ i_S0[1]   ][ i_S0[0] ];
			sxy3 = _sxy[ i_S0[1]+1 ][ i_S0[0] ];

			dfdx = ( (px0 - px3) + (px2 - px1) * 27.0 ) * la_x; 
			dfdy = ( (sxy0 - sxy3) + (sxy2 - sxy1) * 27.0 ) * la_y;
			
			_vx_x_IV[ i_[1] ][ i_[0] ] = ( _vx_x_IV[ i_[1] ][ i_[0] ] * (1.0f - etaxdt) 
										+ dfdx*_mv0[ i_[1] ][ i_[0] ] ) / (1.0f + etaxdt);
			_vx_y_IV[ i_[1] ][ i_[0] ] =   _vx_y_IV[ i_[1] ][ i_[0] ] + dfdy*_mv0[ i_[1] ][ i_[0] ];
		
			_vx[ i_[1] ][ i_[0] ] = _vx_x_IV[ i_[1] ][ i_[0] ] + _vx_y_IV[ i_[1] ][ i_[0] ];
		}
	}
    	/*************************************************************************************/
#ifdef SUPER_MARIO_VERBOSE
	fprintf(stderr,"out of pml_III + pml_IV + physical loops\n");
#endif

    	/* pml region II ********************************************************************/
	if ( pml_II_width_y>0 ){

		//computing loop widths and offsets
		width[0] = n_[0];
		width[1] = pml_II_width_y;

		offset_[0]   = gs_[0];
		offset_P0[0] = gs_P0[0];
		offset_S0[0] = gs_S0[0];

		offset_[1]   = n_[1]   - pml_II_width_y + gs_[1]   + tid;
		offset_P0[1] = n_[1] - pml_II_width_y + gs_P0[1] + tid;
		offset_S0[1] = n_[1] - pml_II_width_y + gs_S0[1] + tid;

		//main loop in y-direction
		for ( index[1] = tid; index[1] < width[1]; index[1]++ ){
			i_[1]   = index[1] + offset_[1];
			i_P0[1] = index[1] + offset_P0[1];
			i_S0[1] = index[1] + offset_S0[1];

			etaydt = _epy[i_[1]] * dt2;
	
			//loop in x-direction
			for ( index[0] = 0; index[0] < width[0]; index[0]++ ){
				i_[0]   = index[0] + offset_[0];
				i_P0[0] = index[0] + offset_P0[0];
				i_S0[0] = index[0] + offset_S0[0];

				etaxdt = _epx[i_[0]] * dt2;
				
				px0 = _px[ i_P0[1] ][ i_P0[0]-1 ];
				px1 = _px[ i_P0[1] ][ i_P0[0]   ];
				px2 = _px[ i_P0[1] ][ i_P0[0]+1 ];
				px3 = _px[ i_P0[1] ][ i_P0[0]+2 ];

				sxy0 = _sxy[ i_S0[1]-2 ][ i_S0[0] ];
				sxy1 = _sxy[ i_S0[1]-1 ][ i_S0[0] ];
				sxy2 = _sxy[ i_S0[1]   ][ i_S0[0] ];
				sxy3 = _sxy[ i_S0[1]+1 ][ i_S0[0] ];
				
				dfdx = ( (px0 - px3) + (px2 - px1) * 27.0 ) * la_x; 
				dfdy = ( (sxy0 - sxy3) + (sxy2 - sxy1) * 27.0 ) * la_y;

				_vx_x_II[ i_[1] ][ i_[0] ] = ( _vx_x_II[ i_[1] ][ i_[0] ]*(1.0f - etaxdt) 
                                                                                      + dfdx*_mv0[ i_[1] ][ i_[0] ] ) / (1.0f + etaxdt);
				_vx_y_II[ i_[1] ][ i_[0] ] = ( _vx_y_II[ i_[1] ][ i_[0] ]*(1.0f - etaydt) 
                                                                                      + dfdy*_mv0[ i_[1] ][ i_[0] ] ) / (1.0f + etaydt);

				_vx[ i_[1] ][ i_[0] ] = _vx_x_II[ i_[1] ][ i_[0] ] + _vx_y_II[ i_[1] ][ i_[0] ];
			}
		}
	}
    	/*************************************************************************************/
#ifdef SUPER_MARIO_VERBOSE
	fprintf(stderr,"out of pml_II loops\n");
#endif

	} /* omp parallel */

  	return 0;
}

/*----------------------------------------------------------------------------*/
int esg_2d_4_v1( float ** _vy,  float ** _mv1, 
	         float ** _sxy, float ** _py,
		 float * _epx,  float * _epy,
		 float ** _vy_x_I, float ** _vy_x_II, float ** _vy_x_III, float ** _vy_x_IV,
		 float ** _vy_y_I, float ** _vy_y_II, float ** _vy_y_III, float ** _vy_y_IV,
		 int pml_I_width_y, int pml_II_width_y, int pml_III_width_x, int pml_IV_width_x,
		 int * n_,  int * n_S0,  int * n_P1,
		 int * gs_, int * gs_S0, int * gs_P1,
		 float la_x, float la_y, float dt     ){
/*----------------------------------------------------------------------------*/
#ifdef MARIO_VERBOSE
	fprintf(stderr,"> esg_2d_4_v1\n");
#endif

	int tsz, tid;

	int width[] = {0,0};

	int offset_[]   = {0,0}; 
	int offset_S0[] = {0,0};
	int offset_P1[] = {0,0};

	int offset_III[]    = {0,0};
	int offset_S0_III[] = {0,0};
	int offset_P1_III[] = {0,0};

	int offset_IV[]    = {0,0};
	int offset_S0_IV[] = {0,0};
	int offset_P1_IV[] = {0,0};

	int index[] = {0,0};
	int i_[]    = {0,0};
	int i_S0[]  = {0,0};
	int i_P1[]  = {0,0};

	register float sxy3, sxy2, sxy1, sxy0;
	register float py3,  py2,  py1,  py0; 
	register float dfdx, dfdy, etaxdt, etaydt;
	register float aux;
	register float dt2 = dt / 2.0;
	
	
#pragma omp parallel private(tsz,tid,index,i_,i_S0,i_P1,offset_,offset_S0_III,offset_P1,offset_III,offset_S0_IV,offset_P1_IV,_vy,_mv1,_esxy,_epy,_sxy,_py,sxy3,sxy2,sxy1,sxy0,py3,py2,py1,py0,_vy_x_I,_vy_y_I,_vy_x_II,_vy_y_II,_vy_x_III,_vy_y_III,_vy_x_IV,_vy_y_IV,dfdx,dfdy,etaxdt,etaydt)
	{
#ifdef _OPENMP
	tsz = omp_get_num_threads();
	tid = omp_get_thread_num();
#else
	tsz = 1;
	tid = 0;
#endif
	
#pragma omp single

#pragma omp barrier


#ifdef SUPER_MARIO_VERBOSE
	fprintf(stderr,"going into pml\n");
#endif

	//##################################################################################//
	//# PML + PHYSICAL REGIONS #########################################################//
	//##################################################################################//

	/* pml region I *********************************************************************/
	if ( pml_I_width_y>0 ){

		//computing width of loops and offsets
		width[0]  = n_[0];
		width[1]  = pml_I_width_y;
		
		offset_[0]   = gs_[0];
		offset_S0[0] = gs_S0[0];
		offset_P1[0] = gs_P1[0];
		
		offset_[1]   = gs_[1]   + tid;
		offset_S0[1] = gs_S0[1] + tid;
		offset_P1[1] = gs_P1[1] + tid;

		//loop in y-direction
		for ( index[1] = tid; index[1] < width[1]; index[1]++ ){
			i_[1]   = index[1] + offset_[1];
			i_S0[1] = index[1] + offset_S0[1];
			i_P1[1] = index[1] + offset_P1[1];

			etaydt = _epy[i_[1]] * dt2;

			//loop in x-direction
			for ( index[0] = 0; index[0] < width[0]; index[0]++ ){
				i_[0]   = index[0] + offset_ [0];
				i_S0[0] = index[0] + offset_S0[0];
				i_P1[0] = index[0] + offset_P1[0];

				etaxdt = _epx[i_[0]] * dt2;
				
				sxy0 = _sxy[ i_S0[1] ][ i_S0[0]-2 ];
				sxy1 = _sxy[ i_S0[1] ][ i_S0[0]-1 ];
				sxy2 = _sxy[ i_S0[1] ][ i_S0[0]   ];
				sxy3 = _sxy[ i_S0[1] ][ i_S0[0]+1 ];

				py0 = _py[ i_P1[1]-1 ][ i_P1[0] ];
				py1 = _py[ i_P1[1]   ][ i_P1[0] ];
				py2 = _py[ i_P1[1]+1 ][ i_P1[0] ];
				py3 = _py[ i_P1[1]+2 ][ i_P1[0] ];
				
				dfdx = ( (sxy0 - sxy3) + (sxy2 - sxy1) * 27.0 ) * la_x; 
				dfdy = ( (py0 - py3) + (py2 - py1) * 27.0 ) * la_y;

				_vy_x_I[ i_[1] ][ i_[0] ] = ( _vy_x_I[ i_[1] ][ i_[0] ]*(1.0f - etaxdt) 
                                                                                      + dfdx*_mv1[ i_[1] ][ i_[0] ] ) / (1.0f + etaxdt);
				_vy_y_I[ i_[1] ][ i_[0] ] = ( _vy_y_I[ i_[1] ][ i_[0] ]*(1.0f - etaydt) 
                                                                                      + dfdy*_mv1[ i_[1] ][ i_[0] ] ) / (1.0f + etaydt);

				_vy[ i_[1] ][ i_[0] ] = _vy_x_I[ i_[1] ][ i_[0] ] + _vy_y_I[ i_[1] ][ i_[0] ];
			}
		}
	}

#ifdef SUPER_MARIO_VERBOSE
	fprintf(stderr,"out of pml_I loop\n");
#endif

	/* pml region III, IV and physical region *******************************************/
	
	//computing width of loops and offsets
	width[0] = n_[0] - pml_III_width_x - pml_IV_width_x;
	width[1] = n_[1] - pml_I_width_y   - pml_II_width_y;

	offset_ [0]      = pml_III_width_x + gs_ [0];
	offset_S0[0]     = pml_III_width_x + gs_S0[0];
	offset_P1[0]     = pml_III_width_x + gs_P1[0];

	offset_ [1]      = pml_I_width_y + gs_ [1] + tid;
	offset_S0[1]     = pml_I_width_y + gs_S0[1] + tid;
	offset_P1[1]     = pml_I_width_y + gs_P1[1] + tid;

	offset_III [0]   = gs_ [0];
	offset_S0_III[0] = gs_S0[0];
	offset_P1_III[0] = gs_P1[0];

	offset_IV [0]    = n_[0] - pml_IV_width_x + gs_ [0];
	offset_S0_IV[0]  = n_[0] - pml_IV_width_x + gs_S0[0];
	offset_P1_IV[0]  = n_[0] - pml_IV_width_x + gs_P1[0];

	//loop in y-direction
	for ( index[1] = tid; index[1] < width[1]; index[1]++ ){
		i_ [1]  = index[1] + offset_ [1];
		i_S0[1] = index[1] + offset_S0[1];
		i_P1[1] = index[1] + offset_P1[1];

		etaydt = _epy[i_[1]] * dt2;

		/* pml region III */
		for ( index[0] = 0; index[0] < pml_III_width_x; index[0]++ ){
			i_ [0]  = index[0] + offset_III [0];
			i_S0[0] = index[0] + offset_S0_III[0];
			i_P1[0] = index[0] + offset_P1_III[0];

			etaxdt = _epx[i_[0]] * dt2;
			
			sxy0 = _sxy[ i_S0[1] ][ i_S0[0]-2 ];
			sxy1 = _sxy[ i_S0[1] ][ i_S0[0]-1 ];
			sxy2 = _sxy[ i_S0[1] ][ i_S0[0]   ];
			sxy3 = _sxy[ i_S0[1] ][ i_S0[0]+1 ];

			py0 = _py[ i_P1[1]-1 ][ i_P1[0] ];
			py1 = _py[ i_P1[1]   ][ i_P1[0] ];
			py2 = _py[ i_P1[1]+1 ][ i_P1[0] ];
			py3 = _py[ i_P1[1]+2 ][ i_P1[0] ];

			dfdx = ( (sxy0 - sxy3) + (sxy2 - sxy1) * 27.0 ) * la_x; 
			dfdy = ( (py0 - py3) + (py2 - py1) * 27.0 ) * la_y;
			
			_vy_x_III[ i_[1] ][ i_[0] ] = ( _vy_x_III[ i_[1] ][ i_[0] ]*(1.0f - etaxdt) 
										  + dfdx*_mv1[ i_[1] ][ i_[0] ] ) / (1.0f + etaxdt);
			
			_vy_y_III[ i_[1] ][ i_[0] ] =   _vy_y_III[ i_[1] ][ i_[0] ] + dfdy*_mv1[ i_[1] ][ i_[0] ];
		
			_vy[ i_[1] ][ i_[0] ] = _vy_x_III[ i_[1] ][ i_[0] ] + _vy_y_III[ i_[1] ][ i_[0] ];
		}
	
		/* physical region */
		for ( index[0] = 0; index[0] < width[0]; index[0]++ ){
			i_ [0]  = index[0] + offset_ [0];
			i_S0[0] = index[0] + offset_S0[0];
			i_P1[0] = index[0] + offset_P1[0];

			sxy0 = _sxy[ i_S0[1] ][ i_S0[0]-2 ];
			sxy1 = _sxy[ i_S0[1] ][ i_S0[0]-1 ];
			sxy2 = _sxy[ i_S0[1] ][ i_S0[0]   ];
			sxy3 = _sxy[ i_S0[1] ][ i_S0[0]+1 ];

			py0 = _py[ i_P1[1]-1 ][ i_P1[0] ];
			py1 = _py[ i_P1[1]   ][ i_P1[0] ];
			py2 = _py[ i_P1[1]+1 ][ i_P1[0] ];
			py3 = _py[ i_P1[1]+2 ][ i_P1[0] ];
			
			dfdx = ( (sxy0 - sxy3) + (sxy2 - sxy1) * 27.0 ) * la_x; 
			dfdy = ( (py0 - py3) + (py2 - py1) * 27.0 ) * la_y;
			
			_vy[ i_[1] ][ i_[0] ] = _vy[ i_[1] ][ i_[0] ] + dfdx*_mv1[ i_[1] ][ i_[0] ] + dfdy*_mv1[ i_[1] ][ i_[0] ];
		}
	
		/* pml region IV */
		for ( index[0] = 0; index[0] < pml_IV_width_x; index[0]++ ){
			i_ [0]  = index[0] + offset_IV [0];
			i_S0[0] = index[0] + offset_S0_IV[0];
			i_P1[0] = index[0] + offset_P1_IV[0];

			etaxdt = _epx[i_[0]] * dt2;

			sxy0 = _sxy[ i_S0[1] ][ i_S0[0]-2 ];
			sxy1 = _sxy[ i_S0[1] ][ i_S0[0]-1 ];
			sxy2 = _sxy[ i_S0[1] ][ i_S0[0]   ];
			sxy3 = _sxy[ i_S0[1] ][ i_S0[0]+1 ];

			py0 = _py[ i_P1[1]-1 ][ i_P1[0] ];
			py1 = _py[ i_P1[1]   ][ i_P1[0] ];
			py2 = _py[ i_P1[1]+1 ][ i_P1[0] ];
			py3 = _py[ i_P1[1]+2 ][ i_P1[0] ];

			dfdx = ( (sxy0 - sxy3) + (sxy2 - sxy1) * 27.0 ) * la_x; 
			dfdy = ( (py0 - py3) + (py2 - py1) * 27.0 ) * la_y;
			
			_vy_x_IV[ i_[1] ][ i_[0] ] = ( _vy_x_IV[ i_[1] ][ i_[0] ] * (1.0f - etaxdt) 
										+ dfdx*_mv1[ i_[1] ][ i_[0] ] ) / (1.0f + etaxdt);
			_vy_y_IV[ i_[1] ][ i_[0] ] =   _vy_y_IV[ i_[1] ][ i_[0] ] + dfdy*_mv1[ i_[1] ][ i_[0] ];
		
			_vy[ i_[1] ][ i_[0] ] = _vy_x_IV[ i_[1] ][ i_[0] ] + _vy_y_IV[ i_[1] ][ i_[0] ];
		}
	}
    	/*************************************************************************************/
#ifdef SUPER_MARIO_VERBOSE
	fprintf(stderr,"out of pml_III + pml_IV + physical loops\n");
#endif

    	/* pml region II ********************************************************************/
	if ( pml_II_width_y>0 ){

		//computing loop widths and offsets
		width[0] = n_[0];
		width[1] = pml_II_width_y;

		offset_[0]   = gs_[0];
		offset_S0[0] = gs_S0[0];
		offset_P1[0] = gs_P1[0];

		offset_[1]   = n_[1] - pml_II_width_y + gs_[1]   + tid;
		offset_S0[1] = n_[1] - pml_II_width_y + gs_S0[1] + tid;
		offset_P1[1] = n_[1] - pml_II_width_y + gs_P1[1] + tid;

		//main loop in y-direction
		for ( index[1] = tid; index[1] < width[1]; index[1]++ ){
			i_[1]   = index[1] + offset_[1];
			i_S0[1] = index[1] + offset_S0[1];
			i_P1[1] = index[1] + offset_P1[1];

			etaydt = _epy[i_[1]] * dt2;
	
			//loop in x-direction
			for ( index[0] = 0; index[0] < width[0]; index[0]++ ){
				i_[0]   = index[0] + offset_[0];
				i_S0[0] = index[0] + offset_S0[0];
				i_P1[0] = index[0] + offset_P1[0];

				etaxdt = _epx[i_[0]] * dt2;
				
				sxy0 = _sxy[ i_S0[1] ][ i_S0[0]-2 ];
				sxy1 = _sxy[ i_S0[1] ][ i_S0[0]-1 ];
				sxy2 = _sxy[ i_S0[1] ][ i_S0[0]   ];
				sxy3 = _sxy[ i_S0[1] ][ i_S0[0]+1 ];

				py0 = _py[ i_P1[1]-1 ][ i_P1[0] ];
				py1 = _py[ i_P1[1]   ][ i_P1[0] ];
				py2 = _py[ i_P1[1]+1 ][ i_P1[0] ];
				py3 = _py[ i_P1[1]+2 ][ i_P1[0] ];
				
				dfdx = ( (sxy0 - sxy3) + (sxy2 - sxy1) * 27.0 ) * la_x; 
				dfdy = ( (py0 - py3) + (py2 - py1) * 27.0 ) * la_y;

				_vy_x_II[ i_[1] ][ i_[0] ] = ( _vy_x_II[ i_[1] ][ i_[0] ]*(1.0f - etaxdt) 
                                                                                      + dfdx*_mv1[ i_[1] ][ i_[0] ] ) / (1.0f + etaxdt);
				_vy_y_II[ i_[1] ][ i_[0] ] = ( _vy_y_II[ i_[1] ][ i_[0] ]*(1.0f - etaydt) 
                                                                                      + dfdy*_mv1[ i_[1] ][ i_[0] ] ) / (1.0f + etaydt);

				_vy[ i_[1] ][ i_[0] ] = _vy_x_II[ i_[1] ][ i_[0] ] + _vy_y_II[ i_[1] ][ i_[0] ];
			}
		}
	}
    	/*************************************************************************************/
#ifdef SUPER_MARIO_VERBOSE
	fprintf(stderr,"out of pml_II loops\n");
#endif

	} /* omp parallel */

  	return 0;
}

